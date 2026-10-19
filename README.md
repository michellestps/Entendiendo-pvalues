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

### 1. Instalar los programas necesarios

Antes de comenzar, instala los siguientes programas:

- **R:** descarga R desde [CRAN](https://cran.r-project.org/).
- **RStudio Desktop:** descarga [RStudio Desktop](https://posit.co/download/rstudio-desktop/).

Instala primero R y después RStudio. Ambos programas son gratuitos.

### 2. Descargar los archivos del taller

Descarga este repositorio desde GitHub mediante **Code → Download ZIP** y descomprime el archivo en tu computadora.

También puedes clonarlo con Git si ya lo tienes instalado:

```bash
git clone <URL-de-este-repo>
```

### 3. Organizar los archivos

Dentro de la carpeta del repositorio, crea una carpeta llamada `datos/` y coloca en ella el archivo de datos descargado desde Google Drive.

La estructura debe quedar así:

```
entendiendo_pvalues/
├── entendiendo_pvalues.R
└── datos/
    └── GSE9348_series_matrix.txt.gz
```

No descomprimas el archivo `.gz`, a menos que el código del taller lo indique.

### 4. Abrir el proyecto en RStudio

Abre el archivo `entendiendo_pvalues.R` desde RStudio.

En el menú, selecciona:

**Session → Set Working Directory → To Source File Location**

Esto permite que R encuentre los archivos utilizando las rutas relativas indicadas en el código.

### 5. Ejecutar el código

Ejecuta el script sección por sección para seguir las explicaciones del taller:

- **Windows/Linux:** `Ctrl + Enter`
- **macOS:** `Cmd + Enter`

Lee los comentarios del código y observa los resultados en la consola de RStudio antes de pasar a la siguiente sección.

## Qué hace el script

| Parte | Contenido |
|---|---|
| 0–3 | Lee la *series matrix* de GEO, identifica los grupos sano y enfermo y toma 8 muestras por grupo. |
| A | Calcula el p-value de un gen con una prueba de permutación (9 999 permutaciones) y lo compara con `t.test()`. |
| B | Hace una prueba t en cada una de las ~54 000 sondas. Repite el análisis con las etiquetas mezcladas para ver los falsos positivos. Aplica `p.adjust()` con Bonferroni y BH. |
| C | Muestra la tabla de resultados, sondas con p < 0.05 que no pasan FDR y un gráfico volcán. |
| Ejercicios | Propuestas para seguir explorando en casa. |

> **Nota:** el análisis es didáctico. Un análisis diferencial formal usaría herramientas como `limma`, control de calidad y covariables.
