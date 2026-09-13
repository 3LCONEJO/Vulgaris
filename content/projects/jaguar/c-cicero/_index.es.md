---
title: "C — Cicero"
weight: 30
eyebrow: "rama C // co-accesibilidad"
estado: "corre en bash, sin migrar"
---

Segunda ruta hacia la red. En vez de tópicos, busca **co-accesibilidad**: qué pares de
picos se abren al mismo tiempo. Como los datos no tienen dimensión temporal, se calcula
un consenso sobre todas las células.

De ahí salen los **CCANs** (*Cis-Co-Accessibility Networks*): bloques de picos cercanos
que interactúan entre sí, no simples pares.

## Los pasos

| # | Script | Qué hace |
|---|---|---|
| 1 | `01_create_cicero_cds.r` | Seurat → CDS Monocle3 → CDS Cicero → co-accesibilidad (glasso por ventana deslizante) → CCANs (Louvain) → gene activity |
| 2 | `02_TestingCCANSandCONNS.R` | Track plots por CCAN y cromosoma |
| 3 | `03_assign_ccan_genes.r` | Genes por CCAN vía dominio regulatorio de GREAT (hasta 1Mb) |
| 4 | `04_annotate_ccans_go.r` | Enriquecimiento GO por CCAN (array de LSF) |
| 5 | `05_cluster_ccans_by_go.r` | CCAN × término-GO (TF-IDF) → PCA → UMAP → Leiden |
| 6 | `06_build_ccan_index_go.r` | Manifiesto JSON final por CCAN |

Aparte de la secuencia numerada hay scripts de anotación de motivos (HOMER),
RegulomeDB, y comparación de co-accesibilidad entre tipos celulares con test de
permutación por par de linajes.

## Dos cosas pendientes

**El peakset que usa no es el bueno.** `01_create_cicero_cds.r` no lee ningún BED:
toma los picos que ya trae el assay ATAC del objeto Seurat (210,045 picos de Signac).
[A_pyCisTopic](../a-pycistopic/) produce uno mejor — 487,005 picos, 2.3× más denso, con
resolución por tipo celular. Falta un script puente que recuantifique el assay contra
el consenso nuevo, usando los mismos fragmentos que ya se usan.

**Los scripts no son conscientes de las fechas.** A diferencia de los `.py` de
pycisTopic (que aceptan `--out-dir`), los ~8 `.r` de esta rama derivan su carpeta de
salida internamente y no se pueden redirigir. Se parcheó con symlinks de
compatibilidad, pero cualquier prefijo nuevo que se genere cae en una ruta plana hasta
que se les agregue un override de verdad.
