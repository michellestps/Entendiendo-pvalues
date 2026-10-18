# =============================================================================
#  Entendiendo los p-values: Estadística Aplicada a Bioinformática
#  Script único de la sesión (~40 min)
#
#  Datos: GEO GSE9348 — mucosa de colon de controles sanos (n = 12) vs.
#         tumor de cáncer colorrectal en estadio temprano (n = 70).
#         Plataforma GPL570 (Affymetrix HG-U133 Plus 2.0).
#
#  Requisitos: solo R base. No hace falta instalar ningún paquete.
#
#  ANTES DE EMPEZAR
#  1. Descarga desde la carpeta de Drive del curso el archivo
#        GSE9348_series_matrix.txt.gz
#     (es el mismo que ofrece GEO en
#      https://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9348/matrix/ )
#  2. Guárdalo en una carpeta llamada "datos" junto a este script
#     (o cambia la ruta en la variable `archivo_geo` de abajo).
#  3. En RStudio: Session > Set Working Directory > To Source File Location.
#  4. Ejecuta el script sección por sección (Ctrl/Cmd + Enter).
# =============================================================================


# ---- 0. Configuración --------------------------------------------------------

archivo_geo   <- file.path("datos", "GSE9348_series_matrix.txt.gz")
n_por_grupo   <- 8       # muestras por grupo para el ejemplo (8 sanas, 8 enfermas)
B             <- 9999    # número de permutaciones
alfa          <- 0.05    # umbral convencional
semilla       <- 2026    # para que todos obtengamos los mismos resultados

# Sonda de ejemplo para la prueba de un solo gen.
# 204259_at corresponde a MMP7 en GPL570 (gen descrito como sobreexpresado en
# cáncer colorrectal). Si no está en el archivo, el script elige otra sonda.
sonda_ejemplo <- "204259_at"

set.seed(semilla)


# ---- 1. Leer la matriz de expresión de GEO -----------------------------------
# Una "series matrix" de GEO tiene dos partes:
#   * líneas de metadatos que empiezan con "!" (descripción de cada muestra)
#   * una tabla: filas = sondas (genes), columnas = muestras (GSM...)

if (!file.exists(archivo_geo)) {
  stop("No encuentro '", archivo_geo, "'. Descárgalo de la carpeta de Drive ",
       "y colócalo en la carpeta 'datos' junto a este script.")
}

lineas <- readLines(gzfile(archivo_geo))

# La tabla: read.delim ignora las líneas que empiezan con "!"
expr <- read.delim(archivo_geo, comment.char = "!", row.names = 1,
                   check.names = FALSE)
expr <- as.matrix(expr)

cat("Matriz de expresión:", nrow(expr), "sondas x", ncol(expr), "muestras\n")


# ---- 2. Identificar los grupos (sano / enfermo) -------------------------------
# Juntamos todos los metadatos de muestra (!Sample_...) y buscamos palabras clave.

metadatos <- grep("^!Sample_(title|source_name|characteristics|description)",
                  lineas, value = TRUE)
campos <- lapply(strsplit(metadatos, "\t"), function(x) gsub('"', "", x[-1]))
texto_muestra <- do.call(paste, campos)   # un texto por muestra

es_sano    <- grepl("healthy|normal|control", texto_muestra, ignore.case = TRUE)
es_enfermo <- grepl("cancer|tumou?r|crc|carcinoma", texto_muestra, ignore.case = TRUE)

grupo <- ifelse(es_sano & !es_enfermo, "Sano",
         ifelse(es_enfermo & !es_sano, "Enfermo", NA))

if (anyNA(grupo)) {
  print(texto_muestra[is.na(grupo)])
  stop("Hay muestras que no pude clasificar. Revisa los metadatos de arriba.")
}

grupo <- factor(grupo, levels = c("Sano", "Enfermo"))
print(table(grupo))     # para GSE9348 se espera: 12 sanas y 70 enfermas


# ---- 3. Preparar los datos ---------------------------------------------------
# Los microarreglos se analizan en escala log2. Si los valores no lo están,
# los transformamos (sumando 1 para evitar log de 0).

if (max(expr, na.rm = TRUE) > 50) {
  cat("Valores en escala lineal -> aplicamos log2(x + 1)\n")
  expr <- log2(pmax(expr, 0) + 1)
}

expr <- expr[complete.cases(expr), ]

# Para que el ejemplo coincida con la clase (8 vs 8), tomamos al azar
# 8 muestras de cada grupo. Cambia `n_por_grupo` para usar más.
idx_sano    <- sample(which(grupo == "Sano"),    n_por_grupo)
idx_enfermo <- sample(which(grupo == "Enfermo"), n_por_grupo)
idx         <- c(idx_sano, idx_enfermo)

x   <- expr[, idx]
grp <- droplevels(grupo[idx])
cat("Muestras usadas en el ejemplo:\n"); print(table(grp))


# =============================================================================
#  PARTE A (11–20 min): el p-value de UN gen mediante permutaciones
# =============================================================================

if (!sonda_ejemplo %in% rownames(x)) {
  # Si la sonda no está, usamos la de mayor variabilidad
  sonda_ejemplo <- names(which.max(apply(x, 1, var)))
}
y <- x[sonda_ejemplo, ]
cat("\nGen/sonda de ejemplo:", sonda_ejemplo, "\n")

# Estadístico: diferencia de medias (Enfermo - Sano)
dif_medias <- function(valores, etiquetas) {
  mean(valores[etiquetas == "Enfermo"]) - mean(valores[etiquetas == "Sano"])
}

d_obs <- dif_medias(y, grp)
cat("Diferencia observada (log2):", round(d_obs, 3), "\n")

# Distribución nula: mezclamos las etiquetas B veces y recalculamos
d_perm <- replicate(B, dif_medias(y, sample(grp)))

# p-value bilateral: ¿qué proporción de permutaciones es al menos tan extrema?
b       <- sum(abs(d_perm) >= abs(d_obs))
p_perm  <- (b + 1) / (B + 1)            # corrección Monte Carlo (b + 1)/(B + 1)
cat("Permutaciones al menos tan extremas:", b, "de", B, "\n")
cat("p-value por permutación:", signif(p_perm, 3), "\n")

hist(d_perm, breaks = 50, col = "grey85", border = "white",
     main = paste("Distribución nula por permutación -", sonda_ejemplo),
     xlab = "Diferencia de medias con etiquetas mezcladas (log2)",
     xlim = range(c(d_perm, d_obs, -d_obs)))
abline(v = c(d_obs, -d_obs), col = "firebrick", lwd = 2, lty = c(1, 2))
legend("topright", bty = "n", col = "firebrick", lwd = 2, lty = c(1, 2),
       legend = c("Diferencia observada", "Mismo valor, signo contrario"))

# Comparación con la prueba t de Welch clásica (otra forma de obtener un p-value)
print(t.test(y ~ grp))

# PREGUNTA: si cambias B a 99999, ¿el p-value necesariamente baja?


# =============================================================================
#  PARTE B (20–29 min): miles de pruebas y comparaciones múltiples
# =============================================================================

# Prueba t de Welch para TODAS las sondas a la vez (vectorizada con R base;
# da el mismo resultado que t.test() fila por fila, pero en segundos).
welch_por_fila <- function(m, etiquetas) {
  a <- m[, etiquetas == "Enfermo", drop = FALSE]
  s <- m[, etiquetas == "Sano",    drop = FALSE]
  na <- ncol(a); ns <- ncol(s)
  va <- apply(a, 1, var) / na
  vs <- apply(s, 1, var) / ns
  t  <- (rowMeans(a) - rowMeans(s)) / sqrt(va + vs)
  gl <- (va + vs)^2 / (va^2 / (na - 1) + vs^2 / (ns - 1))
  data.frame(log2FC = rowMeans(a) - rowMeans(s),
             t = t,
             p = 2 * pt(-abs(t), gl))
}

res <- welch_por_fila(x, grp)
res <- res[is.finite(res$p), ]           # quitamos sondas sin variabilidad

# Comprobación: coincide con t.test() para la sonda de ejemplo
cat("p (vectorizado) =", signif(res[sonda_ejemplo, "p"], 4),
    " | p (t.test) =", signif(t.test(y ~ grp)$p.value, 4), "\n")

m <- nrow(res)
cat("\nNúmero de pruebas (sondas):", m, "\n")
cat("Sondas con p < 0.05 sin ajustar:", sum(res$p < alfa), "\n")
cat("Esperadas por azar si TODAS las H0 fueran ciertas: ~", round(alfa * m), "\n")

# --- Experimento clave: etiquetas mezcladas (no hay diferencia real) ----------
grp_mezclado <- sample(grp)
res_nulo <- welch_por_fila(x, grp_mezclado)
res_nulo <- res_nulo[is.finite(res_nulo$p), ]
cat("\nCON ETIQUETAS MEZCLADAS (sin señal real):\n")
cat("  Sondas con p < 0.05:", sum(res_nulo$p < alfa), "<- ¡falsos positivos!\n")
cat("  Sondas con FDR (BH) < 0.05:",
    sum(p.adjust(res_nulo$p, "BH") < alfa), "\n")
# Nota: los genes están correlacionados entre sí, así que este número
# varía bastante de una mezcla a otra. Prueba a repetirlo.

op <- par(mfrow = c(1, 2))
hist(res$p, breaks = 40, col = "steelblue", border = "white",
     main = "Etiquetas reales", xlab = "p-value")
hist(res_nulo$p, breaks = 40, col = "grey70", border = "white",
     main = "Etiquetas mezcladas", xlab = "p-value")
par(op)
# Con etiquetas reales: pico cerca de 0 (señal) + fondo plano.
# Con etiquetas mezcladas: histograma aproximadamente plano (lo esperado bajo H0).

# --- Corrección por comparaciones múltiples (R base: p.adjust) ----------------
res$p_bonf <- p.adjust(res$p, method = "bonferroni")
res$FDR_BH <- p.adjust(res$p, method = "BH")

resumen <- data.frame(
  criterio = c("p < 0.05 sin ajustar", "Bonferroni < 0.05", "FDR (BH) < 0.05"),
  sondas   = c(sum(res$p < alfa), sum(res$p_bonf < alfa), sum(res$FDR_BH < alfa))
)
print(resumen)
cat("Umbral de Bonferroni equivalente: p <", signif(alfa / m, 2), "\n")


# =============================================================================
#  PARTE C (29–35 min): leer resultados sin sobreinterpretarlos
# =============================================================================

res <- res[order(res$p), ]
tabla <- data.frame(sonda   = rownames(res),
                    log2FC  = round(res$log2FC, 2),
                    p       = signif(res$p, 3),
                    FDR_BH  = signif(res$FDR_BH, 3),
                    p_bonf  = signif(res$p_bonf, 3))
print(head(tabla, 10), row.names = FALSE)

# Ejemplo del "gen B" de la clase: p < 0.05 pero FDR > 0.05
gen_B <- tabla[tabla$p < alfa & tabla$FDR_BH > alfa, ]
cat("\nSondas con p < 0.05 que NO pasan FDR < 0.05:", nrow(gen_B), "\n")
print(head(gen_B, 3), row.names = FALSE)

# Significancia vs. tamaño del efecto: gráfico volcán
plot(res$log2FC, -log10(res$p), pch = 16, cex = 0.4,
     col = ifelse(res$FDR_BH < alfa, "firebrick", "grey60"),
     xlab = "log2 fold change (Enfermo - Sano)", ylab = "-log10(p)",
     main = "Gráfico volcán (rojo: FDR < 0.05)")
abline(v = c(-1, 1), lty = 2, col = "grey40")

# Recordatorios para la discusión:
#  * FDR < 0.05 NO significa que cada gen tenga 5 % de probabilidad de ser falso.
#  * Un p pequeño no mide el tamaño del efecto: mira log2FC.
#  * Un gen significativo no es todavía un biomarcador clínico.
#  * Este análisis es didáctico: un análisis formal de microarreglos usaría,
#    por ejemplo, limma (Bioconductor), control de calidad y covariables.


# =============================================================================
#  EJERCICIOS PARA CASA
# =============================================================================
# 1. Cambia `n_por_grupo` a 4 y luego a 12. ¿Qué pasa con el número de sondas
#    con FDR < 0.05? ¿Por qué?
# 2. Cambia `semilla` y vuelve a ejecutar. ¿Cambia mucho el p-value de la
#    sonda de ejemplo? ¿Y el número de falsos positivos con etiquetas mezcladas?
# 3. Elige otra sonda (por ejemplo, la primera de `tabla`) y repite la Parte A.
# 4. Con etiquetas mezcladas, ¿cuántas sondas pasan Bonferroni? ¿Y BH?

sessionInfo()
