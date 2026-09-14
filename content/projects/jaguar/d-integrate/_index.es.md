---
title: "D — Integración"
weight: 40
eyebrow: "rama D // donde se cruzan"
estado: "consenso probado, visualizador abierto"
color: 3
---

Aquí es donde las tres rutas se comparan. Esta rama no tiene datos de entrada propios:
sus entradas son las salidas de [A](../a-pycistopic/), [C](../c-cicero/) y
[G](../g-scglue/).

## El consenso de tres vías

Un **join región-gen** sobre el espacio de regiones de pycisTopic. Cada método aporta
su lista de aristas región→gen:

- **pycisTopic** — pesos de gene activity (del paso 09, hoy en cero porque está roto)
- **Cicero** — CCANs + asignación de genes por dominio GREAT
- **scGLUE** — aristas gen-pico del embedding entrenado

Los pares donde las tres coinciden pasan a la siguiente fase. **Los que no coinciden se
conservan y se etiquetan, no se descartan** — el desacuerdo en sí es informativo.

Probado extremo a extremo (2026-08-24) con Cicero y scGLUE reales: 1,571,845 pares
totales, 8,849 con acuerdo de dos vías, 0 de tres vías — exactamente lo esperado sin el
aporte de pycisTopic.

## Que coincidan no era obvio

Contar pares en bruto no alcanza: los métodos difieren muchísimo en cuántos genes
asignan por región (Cicero, vía GREAT hasta 1Mb, tiene una mediana de 10 y llega a 157;
scGLUE, limitado a cuerpo del gen + promotor de 2kb, tiene mediana 1). Comparar al
nivel del par exacto subestima el acuerdo del método angosto.

Por eso se calculan dos pruebas separadas:

1. **A nivel región** (Fisher exacto sobre las 487,005 regiones de consenso): 10,764
   regiones en común contra 2,857 esperadas por azar — **3.77×**, p≈0.
2. **Acuerdo de gen condicionado** (nula por permutación, no Fisher): de esas 10,764
   regiones compartidas, coinciden también en el gen en el **81.6%** de los casos,
   contra ~0.1% esperado por azar.

Dos métodos completamente independientes — uno basado en co-accesibilidad de ATAC más
dominio genómico, el otro en un embedding entrenado con RNA no pareada — convergen en
las mismas regiones y, dado eso, en el mismo gen. El consenso no es ruido.

## Las pruebas por CCAN

| | Qué hace | Estado |
|---|---|---|
| C1 | Join por traslape genómico entre los dos peaksets | Implementado (innecesario si se aplica el puente de [C](../c-cicero/)) |
| C2 | Test continuo por permutación (β promedio del CCAN vs. nula) | Implementado |
| C3 | Fisher exacto sobre la binarización Otsu | Implementado |
| C4 | Test de proporciones para la etiqueta de tipo celular del tópico | Implementado |
| C5 | Corrección winner's curse (nula del máximo entre tópicos, no Bonferroni) | Implementado |
| C6 | Visualizador de cuatro vías | **No existe** |

## El hueco de verdad: C6

No hay ninguna tabla ni gráfico que junte las cuatro señales de tipo celular
(`predicted.l1` de Azimuth, `winners_celltype`, `topic_celltype_label` y
`scglue_celltype`). Falta:

1. La tabla unificada CCAN × las cuatro columnas + distancia de linaje.
2. Un diagrama alluvial/Sankey de cuatro ejes sobre esa tabla.
3. Una jerarquía tosca de linajes (linfoide / mieloide / otros) para distinguir el
   desacuerdo esperado (CD8T↔NK) del alarmante (B↔NK), en vez de dejarlo a ojo.

Ojo al leer desacuerdos: `scglue_celltype` es de naturaleza distinta a las otras tres.
Es una llamada de especificidad de expresión **por gen**, propagada al CCAN por voto
mayoritario — un paso más indirecta. Si es la única que discrepa, pesa menos que si dos
de las otras tres discrepan entre sí.
