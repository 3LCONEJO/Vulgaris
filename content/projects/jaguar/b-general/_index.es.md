---
title: "B — General"
weight: 20
eyebrow: "rama B // puente de entrada"
estado: "parcialmente migrada"
color: 1
---

No es una rama de análisis: es el puente entre el objeto Seurat integrado que llega y
todo lo demás. Exporta lo que las otras ramas necesitan como entrada.

## Los pasos

| Script | Qué hace | Estado |
|---|---|---|
| `export_seurat_data.r` | Exporta peaks, barcodes y metadata con rutas de fragmentos a partir del objeto Seurat integrado | Migrado a Nextflow — ver [[EXPORT SEURAT DATA]] |
| `Cytoscape.r` | Arma la red final para visualizar | Pendiente: espera su CSV en una ruta vieja, fuera de `2_data/` |

Aunque `export_seurat_data.r` no pertenezca a [A_pyCisTopic](../a-pycistopic/), es el
**primer paso real de todo el pipeline** — sin él no arranca ninguna de las tres rutas.

## Trabajo sin clasificar

Hay dos scripts en `2_python/` que no están conectados a ningún `run_pipeline.sh` ni
documentados en el plan: `compare_topics_vs_poissonvi.py` y
`train_poissonvi_topic_conditioned.py`. Parecen comparar los tópicos de Mallet contra
PoissonVI. Están anotados aquí para que no queden huérfanos.
