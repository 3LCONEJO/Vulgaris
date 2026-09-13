---
title: "JAGUAR-GRN"
date: 2026-07-26
draft: false
status: "in_progress"
summary: "A single-cell ATAC-seq pipeline designed to build a Gene Regulation Network to unveil the regulatory systems on JAGUAR Project Data."
tags: ["scATAC-seq", "JAGUAR", "R", "Python", "Seurat"]
github: "https://github.com/JAGUAR-LATAM/CCAN-Visualizer/"
stack: ["Nextflow", "R", "Python", "Seurat", "LSF"]
mermaid: true
type: "guide"
cascade:
  type: "guide"
levels:
  - id: "1"
    name: "Qué y por qué"
    blurb: "La idea del proyecto y el flujo completo. Empieza aquí."
    href: ""
    hint: "esta página"
  - id: "2"
    name: "Cómo está armado"
    blurb: "Una sección por rama. Dentro de cada una, el pipeline y el código por separado."
    href: "a-pycistopic/"
    hint: "ramas A–G"
  - id: "3"
    name: "Qué significa esta palabra"
    blurb: "Glosario de Nextflow y Bash. Pasa el mouse sobre cualquier enlace de nota para leerlo sin salir de la página."
    href: "sintaxis/"
    hint: "glosario"
---

## La idea

De scATAC-seq a una **red de regulación génica**, por tres rutas independientes que se
cruzan al final.

Los datos dicen qué partes del genoma están abiertas en cada célula. De ahí a saber
*qué regula a qué* hay un salto, y hay más de una forma de darlo. Este pipeline da tres
a la vez:

- **Tópicos LDA** ([rama A](a-pycistopic/)) — conjuntos de regiones que se abren juntas.
- **Co-accesibilidad** ([rama C](c-cicero/)) — pares de picos que se abren al mismo tiempo, agrupados en bloques.
- **Embeddings con RNA** ([rama G](g-scglue/)) — la única ruta que ve expresión, no solo accesibilidad.

Las tres se comparan en un [paso de integración](d-integrate/). Solo los pares
región-gen respaldados por las tres señales pasan a la fase final. **Los pares donde
las señales no coinciden no se tiran: se conservan etiquetados**, porque el desacuerdo
entre métodos independientes también dice algo.

## El flujo

```mermaid
flowchart LR
    B["B · general<br/>exporta del Seurat"] --> A["A · pyCisTopic<br/>tópicos LDA"]
    B --> C["C · Cicero<br/>co-accesibilidad"]
    B --> G["G · scGLUE<br/>embeddings + RNA"]
    A -. "puente pendiente<br/>(peakset compartido)" .-> C
    A --> D["D · integración<br/>consenso de 3 vías"]
    C --> D
    G --> D
    D --> E["E · cisTarget"]
    E --> F["F · SCENIC+"]
    F --> RED["red regulatoria"]

    classDef ok fill:#211d1c,stroke:#92c22c,color:#e6e8e1
    classDef partial fill:#211d1c,stroke:#cad54e,color:#e6e8e1
    classDef open fill:#12120f,stroke:#5f6b62,color:#8d9a8f
    class A,G ok
    class B,C,D partial
    class E,F,RED open
```

## Dónde está cada cosa

| Rama | Qué hace | Estado |
|---|---|---|
| [B — general](b-general/) | Exporta peaks, barcodes y metadata del objeto Seurat integrado | Un paso migrado, uno pendiente |
| [A — pyCisTopic](a-pycistopic/) | Tópicos LDA sobre la accesibilidad | Migrada a Nextflow; el paso 09 está roto |
| [C — Cicero](c-cicero/) | CCANs por co-accesibilidad | Corre en bash; usa el peakset equivocado |
| [G — scGLUE](g-scglue/) | Inferencia regulatoria desde embeddings | Implementada y probada |
| [D — integración](d-integrate/) | Consenso de tres vías + comparación de señales | Consenso probado; falta el visualizador |
| [E — cisTarget](e-cistarget/) | Motivos de factores de transcripción | Sin empezar |
| [F — SCENIC+](f-scenicp/) | Red final | Sin empezar |

Transversal a todas: [la migración a Nextflow](nextflow/) y el [glosario](sintaxis/).

## Cómo leer esto

Tres niveles de zoom. Baja solo hasta donde te interese:

