---
title: "A — pyCisTopic"
weight: 10
eyebrow: "rama A // tópicos LDA"
estado: "migrada a nextflow"
color: 0
---

Primera de las tres rutas hacia la red regulatoria. Modela la accesibilidad de la
cromatina como **tópicos LDA** (Mallet): cada tópico es un conjunto de regiones que
se abren juntas, y cada célula es una mezcla de tópicos.

Es la única rama migrada a Nextflow de punta a punta, así que es la que tiene más
notas. Dos ejes hermanos: [cómo se orquesta](nextflow/) y [qué hace cada script](codigo/).

## Los pasos

| # | Script | Qué hace |
|---|---|---|
| 00 | `export_seurat_data.r` | Vive en [B_general](../b-general/), pero es el primer paso real de todo |
| 01 | `01_pseudobulk_export.py` | BED/bigWig pseudobulk por grupo celular |
| 02 | `02_consensus_peak_calling.py` | MACS3 por grupo + consenso de picos (blacklist hg38) |
| 03 | `03_qc_generate_commands.py` | Genera los comandos de QC, uno por muestra |
| 04 | `04_qc_plots_thresholds.py` | Plots de QC + umbrales automáticos de barcode |
| 05 | `05_build_cistopic_objects.py` | Objetos cisTopic por muestra, Scrublet, modelos por muestra |
| 06 | `06_merge_and_model.py` | Fusiona todo y ajusta Mallet para los 11 valores de `n_topics` |
| 07 | `07_model_selection_umap.py` | Selección de modelo, Leiden, UMAP sobre θ |
| 08 | `08_binarization_topic_qc.py` | Binarización Otsu/Yen de β, métricas de tópico, anotación celular |
| 09 | `09_dars_gene_activity.py` | DARs + gene activity. **Roto**, ver abajo |

## Lo que cambió al migrar

Tres pasos dejaron de ser "un proceso que recorre todo por dentro" y pasaron a "un
proceso por unidad, en paralelo". No cambió la lógica de ningún script — solo cómo se
reparte el trabajo:

- **Por muestra** (`[2026-09-01]`): QC y `BUILD_CISTOPIC_OBJECTS`. Una corrida
  secuencial de ~20 muestras no cabía en ningún límite de tiempo de LSF, y no era
  resumible si el job moría a mitad.
- **Por valor de `n_topics`** (`[2026-09-07]`): el modelado se partió en
  `MERGE_ONLY` + `FIT_ONE_TOPIC` ×10 + `MERGE_TOPIC_STATES` + `COMBINE_MODELS`.
  Los 11 ajustes secuenciales hubieran tardado ~85h; en paralelo el total queda
  acotado por el más lento (`n_topics=50`, ~23h). Diferencia de días, no de minutos.

Cada fan-out necesita su fan-in: ver [[first()]] y [[collect()]] en el glosario, que
son las dos mitades del mismo abanico.

## La rama "solo anotadas"

`[2026-09-08]` — primera vez que el DAG **bifurca** de verdad. De las 177k células del
merge, ~63k no tienen anotación de tipo celular: pasaron el QC de ATAC de este pipeline
pero no el filtrado de label-transfer que produce esa etiqueta. La rama reajusta Mallet
solo sobre las 114k que sí la tienen.

No es "confiar a ciegas en que la calidad es buena": es una subselección concreta y
medida. Corre en paralelo a la original para poder compararlas, todavía no como su
reemplazo, y publica en una carpeta hermana — ver [[Publishdir Replace Not Merge]]
para el bug que casi hace que se pisaran.

## Lo que está roto

El paso **09 (`DARS_GENE_ACTIVITY`) lleva tres intentos, los tres muertos por
recursos, ninguno por lógica**: `TERM_RUNLIMIT` tras 12h (2026-08-20), `TERM_MEMLIMIT`
con 500GB pedidos y 501GB usados (2026-08-30), y `TERM_MEMLIMIT` con 64GB pedidos y
75GB usados bajo Nextflow (2026-09-09).

Importa más de lo que parece: ese paso produce el aporte de pycisTopic al consenso de
tres vías de [D_integrate](../d-integrate/), que **hoy sigue en cero**. No se le sube la
memoria otra vez a ciegas — ya se subió tres veces y el pico real solo crece, señal de
que el problema probablemente no es solo "más memoria". Ver [[DARS GENE ACTIVITY]].
