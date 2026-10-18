# Entendiendo los p-values: Estadística aplicada a bioinformática

Material de la sesión:

- qué es un p-value;
- pruebas de permutación;
- comparaciones múltiples: Bonferroni y Benjamini–Hochberg (FDR).

Todo usa datos de expresión génica de GEO.

## Contenido del repositorio

| Archivo | Descripción |
|---|---|
| `entendiendo_pvalues.R` | Script único de la clase.|
| `diapositivas.pdf` | Diapositivas de la sesión. |
| `temario.md` | Temario y guion de la sesión, bloque por bloque. |

## Datos (Drive)

Los datos no están en este repositorio. Descárgalos de Drive:

**➡ [Descargar `GSE9348_series_matrix.txt.gz` (Drive)](https://drive.google.com/file/d/17T8xLIXicxRN_Hj-V2P9zggUvNzItjFh/view?usp=share_link)**

| Archivo | Origen |
|---|---|
| `GSE9348_series_matrix.txt.gz` | [GEO GSE9348](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9348) |

GSE9348 contiene expresión génica de mucosa de colon:

- 12 controles sanos;
- 70 tumores de cáncer colorrectal en estadio temprano;
- plataforma Affymetrix HG-U133 Plus 2.0 (GPL570).

## Cómo ejecutarlo

1. Instala [R](https://cran.r-project.org/) y [RStudio](https://posit.co/download/rstudio-desktop/).
2. Descarga este repositorio: botón verde **Code → Download ZIP** (o `git clone` si usas línea de comandos).
3. Descarga `GSE9348_series_matrix.txt.gz` del Drive (enlace arriba) y guárdalo dentro de una carpeta `datos/`, junto al script.
4. Abre `entendiendo_pvalues.R` en RStudio, ve a **Session → Set Working Directory → To Source File Location**, y corre el script por bloques con `Ctrl/Cmd + Enter`.

## Qué hace el script

| Parte | Contenido |
|---|---|
| 0–3 | Lee la *series matrix* de GEO, identifica los grupos sano y enfermo y toma 8 muestras por grupo. |
| A | Calcula el p-value de un gen con una prueba de permutación (9 999 permutaciones) y lo compara con `t.test()`. |
| B | Hace una prueba t en cada una de las ~54 000 sondas. Repite el análisis con las etiquetas mezcladas para ver los falsos positivos. Aplica `p.adjust()` con Bonferroni y BH. |
| C | Muestra la tabla de resultados, sondas con p < 0.05 que no pasan FDR y un gráfico volcán. |
| Ejercicios | Propuestas para seguir explorando en casa. |

> **Nota:** el análisis es didáctico. Un análisis diferencial formal usaría herramientas como `limma`, control de calidad y covariables.
