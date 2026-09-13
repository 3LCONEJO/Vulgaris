---
title: "La migración a Nextflow"
weight: 5
eyebrow: "transversal"
estado: "en curso"
---

Notas sobre la migración en sí, no sobre ninguna rama en particular: qué es DSL2, cómo
está armado `nextflow.config`, cómo se validó todo sin gastar cómputo, y los dos bugs
reales que solo aparecieron al mandar trabajos de verdad a la farm.

**Hasta ahora solo terminó el piloto de [A_pyCisTopic](../a-pycistopic/)** (más
`export_seurat_data.r` de [B_general](../b-general/)). Cicero, integración y scGLUE
siguen corriendo con sus `run_pipeline.sh` de bash.

Empieza por [[OVERVIEW]]. Si una palabra no te suena, está en el [glosario](../sintaxis/).
