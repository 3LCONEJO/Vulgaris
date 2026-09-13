---
title: "JAGUAR-GRN"
aliases: ["/projects/jaguar-grn/"]
date: 2026-07-26
draft: false
status: "in_progress"
summary: "A single-cell ATAC-seq pipeline designed to build a Gene Regulation Network to unveil the regulatory systems on JAGUAR Project Data."
tags: ["scATAC-seq", "JAGUAR", "R", "Python", "Seurat"]
github: "https://github.com/JAGUAR-LATAM/CCAN-Visualizer/"
stack: ["Nextflow", "R", "Python", "Seurat", "LSF"]
mermaid: true
type: "guide"
---

## The problem

scATAC-seq tells you which parts of the genome are open in each cell. Getting from
there to **what regulates what** is a leap, and there is more than one way to take it.
This pipeline takes three at once and then compares them.

| Existing capabilities | New GRN capabilities |
|---|---|
| Filter scATAC-seq data by quality | Discover synchronously activated regions |
| Group cells by accessible regions | Map enhancer-to-target interactions |
| Perform differential accessibility | Identify regulatory topics and cellular archetypes |

## Three routes, one consensus

```mermaid
flowchart LR
    B["B · general"] --> A["A · pyCisTopic<br/>LDA topics"]
    B --> C["C · Cicero<br/>co-accessibility"]
    B --> G["G · scGLUE<br/>embeddings + RNA"]
    A --> D["D · integration<br/>three-way consensus"]
    C --> D
    G --> D
    D --> E["E · cisTarget"]
    E --> F["F · SCENIC+"]
    F --> NET["regulatory network"]

    classDef ok fill:#211d1c,stroke:#92c22c,color:#e6e8e1
    classDef partial fill:#211d1c,stroke:#cad54e,color:#e6e8e1
    classDef open fill:#12120f,stroke:#5f6b62,color:#8d9a8f
    class A,G ok
    class B,C,D partial
    class E,F,NET open
```

- **LDA topics** — sets of regions that open together.
- **Co-accessibility** — pairs of peaks that open at the same time, grouped into CCANs.
- **Embeddings with RNA** — the only route that sees expression, not just accessibility.

Only region-gene pairs backed by all three signals move on to the final stage. Pairs
where the signals disagree are **kept and labelled, not discarded** — disagreement
between independent methods is itself informative.

That the three agree at all was not a given. Two of them, built on completely
independent evidence, converge on the same regions **3.77× more often than chance**,
and — given a shared region — on the same gene **81.6% of the time** against ~0.1%
expected. The consensus is not noise.

## The full guide is in Spanish

The complete technical documentation — one section per branch, every Nextflow module,
and a glossary you can read by hovering over any note link — lives on the Spanish side
of this site:

**[→ Leer la guía completa](/Vulgaris/es/projects/jaguar/)**

It is written in Spanish on purpose: it started as my own lab notebook, and that is
where the explanations are clearest. Translating it is a separate job.
