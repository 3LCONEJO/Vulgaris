---
title: "Teporingo"
date: 2026-07-27
draft: false
status: "in progress"
summary: "An academic AI demultiplexing tool leveraging Siamese networks and cosine similarity to map cells to donors."
tags: ["Python", "Complex Systems", "Simulation", "Deep Learning", "PyTorch"]
github: "https://github.com/3LCONEJO/Teporingo/"
stack: ["Python", "PyTorch", "NumPy", "SciPy"]
---

## The Problem

Standard demultiplexing tools (like `demuxlet`) are robust, but they suffer from two major drawbacks when applied to our datasets: they are incredibly slow, and they aren't explicitly designed for the extreme sparsity of scATAC-seq data. 

For my final AI course project, I decided to try to solve this my self. I built **Teporingo** — a deep learning module within the broader Tepolito pipeline (Just Planed, not even started (ง •̀_•́)ง). Teporingo uses a **Siamese Network** to project both the single-cell ATAC reads and the donor VCF genotypes into the same embedded latent space, allowing us to rapidly calculate the similarity (cosine distance) between a cell and a potential donor.

> **Try it out:** You can explore the interactive visualization of this latent space here: [Teporingo Web Demo](https://3lconejo.github.io/Teporingo/)

---

## Approach: Siamese Networks

Instead of calculating rigid statistical likelihoods for every single SNP across every single cell, a Siamese Network learns the fundamental *patterns* of genetic similarity. 

### Pipeline Flow

`[Cell ATAC Matrix]` -> `[Sub-network A]` v 
                                                
    `[Latent Space Embedding] -> Cosine Similarity -> [Match/No Match]`

`[Donor VCF Matrix]` -> `[Sub-network B]` ^

---

## Implementation Guide: The Training Pipeline

Below is the structured PyTorch training pipeline for Teporingo. It is designed to be highly modular, supporting YAML-driven configurations, Leave-One-Out (LOO) cross-validation for batch effects, and automated early stopping.

### 1. Configuration & Setup
The pipeline relies on YAML configuration to manage hyperparameters dynamically. This block merges the config file settings with command-line arguments to build a unified training configuration.

```python
#!/usr/bin/env python3
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from scipy.sparse import load_npz
from torch.optim import AdamW
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch.utils.data import DataLoader, random_split
from tqdm import tqdm

from buildnload import GenotypeMatchingDataset
from teporingo_demultiplexing.config import load_simple_yaml
from teporingo_demultiplexing.pipeline import normalize_pair_mode, run_pipeline
from teporingo_demultiplexing.models.siamese_network import SiameseNetwork, ShallowSiameseNetwork

# --- CONFIGURATION BUILDERS ---
def build_training_config(config):
    """Merge training settings from dedicated and legacy YAML sections."""
    pipeline_config = config.get('pipeline', {}) or {}
    training_config = config.get('training', {}) or config.get('model_training', {}) or {}

    return {
        'output_dir': training_config.get('output_dir', pipeline_config.get('output_dir', 'models/siamese')),
        'split_mode': training_config.get('split_mode', pipeline_config.get('split_mode', 'single-batch')),
        'batch_label': training_config.get('batch_label', pipeline_config.get('batch_label')),
        'epochs': training_config.get('epochs', 50),
        'batch_size': training_config.get('batch_size', 64),
        'learning_rate': training_config.get('learning_rate', 1e-3),
        'optimizer': training_config.get('optimizer', 'adamw'),
        'embedding_dim': training_config.get('embedding_dim', 512),
        'network': training_config.get('network', 'siamese'),
        'patience': training_config.get('patience', 10),
        'seed': training_config.get('seed', 42),
        'pair_mode': normalize_pair_mode(training_config.get('pair_mode', 'sampled')),
    }
```

### 2. Data Loading & Splits

Because scATAC-seq data is highly sparse, we load the BAM matrices using scipy.sparse. The pipeline supports random splitting for standard training, as well as leave-one-batch-out to ensure the model generalizes across different sequencing runs (e.g., Batches A, B, and C).

```python
def load_data(vcf_path, bam_path, metadata_path, assignments_path=None):
    """Load VCF tensors, sparse BAM matrices, and metadata."""
    vcf_data = torch.load(vcf_path, weights_only=False)
    X_VCF, donors = vcf_data['X_VCF'], vcf_data['donors']
    
    X_BAM = load_npz(bam_path)
    
    metadata = torch.load(metadata_path, weights_only=False)
    barcodes = metadata.get('barcodes', [])
    assignments = load_assignments_table(assignments_path) if assignments_path else metadata.get('assignments', {})

    return X_VCF, X_BAM, donors, barcodes, assignments

def create_random_dataloaders(X_VCF, X_BAM, donors, barcodes, assignments, batch_size=64, num_workers=4, seed=42, negative_ratio=1.0, pair_mode='sampled', pin_memory=False, hard_negatives_k=0):
    """Split the dataset into Train (70%), Val (15%), and Test (15%)."""
    full_dataset = GenotypeMatchingDataset(
        X_vcf=X_VCF, X_bam_sparse=X_BAM, cell_to_donor=assignments, barcodes=barcodes,
        donors=donors, negative_ratio=negative_ratio, pair_mode=pair_mode, hard_negatives_k=hard_negatives_k,
    )

    n_total = len(full_dataset)
    n_train, n_val = int(n_total * 0.7), int(n_total * 0.15)
    
    generator = torch.Generator().manual_seed(seed)
    train_ds, val_ds, test_ds = random_split(full_dataset, [n_train, n_val, n_total - n_train - n_val], generator=generator)

    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True, num_workers=num_workers, pin_memory=pin_memory)
    val_loader = DataLoader(val_ds, batch_size=batch_size * 2, shuffle=False, num_workers=num_workers, pin_memory=pin_memory)
    test_loader = DataLoader(test_ds, batch_size=batch_size * 2, shuffle=False, num_workers=num_workers, pin_memory=pin_memory)

    return train_loader, val_loader, test_loader
```

### 3. The Core Training Loop

The network trains using SmoothL1Loss to handle outliers in the similarity scores. During validation, we track both Area Under the ROC Curve (AUC) and Average Precision (AP) to ensure the model distinguishes true donor matches from hard negatives.

```python
def train_epoch(model, loader, criterion, optimizer, device, epoch):
    model.train()
    total_loss, correct, total = 0.0, 0, 0
    pbar = tqdm(loader, desc=f'Epoch {epoch} [Train]')

    for cell_batch, geno_batch, labels in pbar:
        cell_batch, geno_batch, labels = cell_batch.to(device), geno_batch.to(device), labels.to(device)

        predictions = model(cell_batch, geno_batch)
        loss = criterion(predictions, labels)

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()

        total_loss += loss.item() * len(labels)
        correct += ((predictions > 0.5).float() == labels).sum().item()
        total += len(labels)
        pbar.set_postfix({'loss': f'{loss.item():.4f}', 'acc': f'{100 * correct / total:.2f}%'})

    return total_loss / total, correct / total

@torch.no_grad()
def validate_epoch(model, loader, criterion, device, epoch, phase='Val'):
    from sklearn.metrics import average_precision_score, roc_auc_score
    model.eval()
    
    total_loss, correct, total = 0.0, 0, 0
    all_predictions, all_labels = [], []
    pbar = tqdm(loader, desc=f'Epoch {epoch} [{phase}]')

    for cell_batch, geno_batch, labels in pbar:
        cell_batch, geno_batch, labels = cell_batch.to(device), geno_batch.to(device), labels.to(device)

        predictions = model(cell_batch, geno_batch)
        loss = criterion(predictions, labels)

        total_loss += loss.item() * len(labels)
        correct += ((predictions > 0.9).float() == labels).sum().item()
        total += len(labels)

        all_predictions.extend(predictions.cpu().numpy())
        all_labels.extend(labels.cpu().numpy())

    auc = roc_auc_score(all_labels, all_predictions) if len(set(int(l) for l in all_labels)) >= 2 else float('nan')
    ap = average_precision_score(all_labels, all_predictions)

    return total_loss / total, correct / total, auc, ap
```

### 4. Orchestration & Checkpointing

We wrap the epochs in an orchestration function that manages the learning rate scheduler (ReduceLROnPlateau) and enforces early stopping if the validation loss stops improving, saving only the best_model.pt.

```python
def train_model(model, train_loader, val_loader, test_loader, criterion, optimizer, scheduler, device, epochs, output_dir, patience=10):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    best_val_loss = float('inf')
    epochs_no_improve = 0

    for epoch in range(1, epochs + 1):
        train_loss, train_acc = train_epoch(model, train_loader, criterion, optimizer, device, epoch)
        val_loss, val_acc, val_auc, val_ap = validate_epoch(model, val_loader, criterion, device, epoch)

        if isinstance(scheduler, ReduceLROnPlateau):
            scheduler.step(val_loss)
        else:
            scheduler.step()

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            epochs_no_improve = 0
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
            }, output_dir / 'best_model.pt')
            print(f'  ✓ New best model saved (val_loss: {val_loss:.4f})')
        else:
            epochs_no_improve += 1
            if epochs_no_improve >= patience:
                print(f'\nEarly stopping after {epoch} epochs')
                break

    # Final Evaluation
    checkpoint = torch.load(output_dir / 'best_model.pt', weights_only=False)
    model.load_state_dict(checkpoint['model_state_dict'])
    validate_epoch(model, test_loader, criterion, device, epoch='Final', phase='Test')
```