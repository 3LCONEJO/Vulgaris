---
title: "G — scGLUE"
weight: 70
eyebrow: "rama G // la que sí ve expresión"
estado: "implementada y probada"
color: 0
---

Tercera ruta, añadida en `[2026-08-24]`. Hasta entonces el consenso comparaba dos
vistas del **mismo** dato de ATAC: tópicos LDA contra co-accesibilidad. Ninguna de las
dos usa expresión.

scGLUE es la única señal realmente informada por RNA — RNA no pareada con la misma
ATAC, integrada vía embeddings, y de ahí su inferencia regulatoria.

La letra es G y no E porque E y F ya estaban reservadas para cisTarget y SCENIC+.
Conceptualmente es hermana de A y C, no algo que va después; se puso al final de la
lista para no romper la regla de nunca renumerar las letras.

## El hallazgo: la corrida existía, pero estaba a medias

Antes de escribir nada se encontró una corrida de scGLUE ya entrenada sobre las mismas
24 muestras y los mismos 7 tipos celulares de este proyecto (158,019 células de scATAC
+ 123,071 de scRNA no pareadas, precisión de transferencia de tipo celular 0.896).

Rastreando su script línea por línea se confirmó algo importante: el grafo que dejó
guardado es la salida del paso **previo** al entrenamiento — proximidad genómica pura,
cuerpo del gen más promotor de 2kb, sin ninguna señal de expresión. **Nunca se llamó al
paso de inferencia regulatoria**, el que sí usa el embedding entrenado. Simplemente no
se había hecho.

## La limitación que eso impone

Verificado contra el código fuente instalado, no contra el tutorial: el modelo indexa
cada nodo del grafo contra el vocabulario **fijo** con el que se entrenó, y falla si
falta alguno. No generaliza a nodos nuevos.

Consecuencia real: el universo de pares gen-pico que este modelo puede puntuar está
limitado a los ~18,176 nodos del grafo original (4,000 genes + ~14,176 picos), dentro
de una ventana de cuerpo del gen + promotor de 2kb. **No alcanza candidatos
genuinamente distales**, como sí hace el dominio GREAT de [Cicero](../c-cicero/) hasta
1Mb. Es una limitación de adoptar esa corrida tal cual, no un bug.

Se decidió adoptarla de todos modos en vez de reentrenar desde cero.

## Los pasos

| # | Script | Qué hace | Probado |
|---|---|---|---|
| 01 | `01_import_scglue_run.py` | Importa por symlink la corrida externa | Sí |
| 02 | `02_regulatory_inference.py` | El paso que faltaba: inferencia regulatoria sobre el mismo grafo | Sí — 14,454 pares candidatos, 11,066 pasan qval<0.05 |
| 03 | `03_reconcile_peakset.py` | Join por traslape contra el consenso de pycisTopic | Sí — el 99.7% de los picos coincide |
| 04 | `04_celltype_summary.py` | Especificidad de expresión por gen y tipo celular | **No** — cómputo real sobre 38GB, debe ir por LSF |
| 05 | `05_create_TF_GRN.py` | Red TF-gen vía la receta de scGLUE + poda con pySCENIC | Sí, de punta a punta |

## Un detalle para cuando se vuelva a tocar

Los nombres de gen tipo inmunoglobulina (`IGHV7-4-1`) terminan en el mismo patrón
`-<entero>-<entero>` que un ID de pico, y hay picos en scaffolds alternos que ni
siquiera empiezan con `chr`. Intentar distinguir gen de pico **por el nombre falla en
ambas direcciones**.

La solución correcta, verificada contra las 14,454 aristas reales: las aristas de tipo
`fwd` del grafo de guía son siempre gen→pico por construcción. La dirección de la
arista basta — no hay que adivinar.
