# =============================================================================
#  Entendiendo los p-values: Estadística Aplicada a Bioinformática
#
#  Datos: GEO GSE9348 — mucosa de colon de controles sanos (n = 12) vs.
#         tumor de cáncer colorrectal en estadio temprano (n = 70).
#         Plataforma GPL570 (Affymetrix HG-U133 Plus 2.0).
#
#  Requisitos: solo R base. No hace falta instalar ningún paquete.
#
#  Antes de empezar
#  1. Descarga desde la carpeta de Drive del curso el archivo
#        GSE9348_series_matrix.txt.gz
#     (es el mismo que ofrece GEO en
#      https://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9348/matrix/ )
#  2. Guárdalo en una carpeta llamada "datos" junto a este script
#     (o cambia la ruta en la variable `archivo_geo` de abajo).
#  3. En RStudio: Session > Set Working Directory > To Source File Location.
#     (si R no encuentra el archivo)
#  4. Ejecuta el script sección por sección (Ctrl/Cmd + Enter).
# =============================================================================


# ---- 0. Configuración --------------------------------------------------------

archivo_geo   <- file.path("datos", "GSE9348_series_matrix.txt.gz")
n_por_grupo   <- 8       # muestras por grupo para el ejemplo (8 sanas, 8 enfermas).
                          # Se usan pocas a propósito: con una muestra chica el
                          # resultado de un solo gen es más inestable, y eso es
                          # parte de lo que se quiere mostrar en la Parte A.
B             <- 9999    # número de permutaciones para la Parte A. 9999 es
                          # habitual para una prueba de un solo estadístico: es
                          # rápido (ver nota de tiempos abajo) y da resolución
                          # de p-value hasta 0.0001.
alfa          <- 0.05    # umbral convencional.
semilla       <- 2026    # fija el generador de números aleatorios para que,
                          # con la misma versión de R, todos obtengamos los
                          # mismos números.

# 9999 permutaciones de un solo gen toman bien por debajo de un
# segundo. Lo que de verdad tarda unos segundos es leer el archivo de GEO,
# una sola vez, al principio; eso no depende de B.

# Sonda de ejemplo para la prueba de un solo gen.
# 204259_at corresponde a MMP7 en GPL570 (gen descrito como sobreexpresado en
# cáncer colorrectal). Si no está en el archivo, el script elige otra sonda.
sonda_ejemplo <- "204259_at"

set.seed(semilla)


# ---- 1. Leer la matriz de expresión de GEO -----------------------------------
# Una "series matrix" de GEO tiene dos partes en el mismo archivo de texto:
#   * líneas de metadatos que empiezan con "!" (quién es cada muestra: sana,
#     tumor, edad, sexo, etc.)
#   * una tabla: filas = sondas, columnas = muestras (identificadas como GSM...)
# Se leen las dos partes por separado: `lineas` para los metadatos de texto,
# y `read.delim` para la tabla numérica (ignorando las líneas que empiezan
# con "!", que es justo lo que hace `comment.char = "!"`).

if (!file.exists(archivo_geo)) {
  stop("No encuentro '", archivo_geo, "'. Descárgalo de la carpeta de Drive ",
       "y colócalo en la carpeta 'datos' junto a este script.")
}

lineas <- readLines(gzfile(archivo_geo))

expr <- read.delim(archivo_geo, comment.char = "!", row.names = 1,
                   check.names = FALSE)
expr <- as.matrix(expr)

cat("Matriz de expresión:", nrow(expr), "sondas x", ncol(expr), "muestras\n")


# ---- 2. Identificar los grupos (sano / enfermo) -------------------------------
# GEO no trae una columna lista que diga "Sano" o "Enfermo": hay que
# construirla a partir del texto que escribió quien subió los datos.
# Por eso se juntan varios campos de metadatos (título, fuente del tejido,
# características, descripción) en un solo texto por muestra, y se buscan
# palabras clave en ese texto.

metadatos <- grep("^!Sample_(title|source_name|characteristics|description)",
                  lineas, value = TRUE)
campos <- lapply(strsplit(metadatos, "\t"), function(x) gsub('"', "", x[-1]))
texto_muestra <- do.call(paste, campos)   # un texto por muestra

es_sano    <- grepl("healthy|normal|control", texto_muestra, ignore.case = TRUE)
es_enfermo <- grepl("cancer|tumou?r|crc|carcinoma", texto_muestra, ignore.case = TRUE)

grupo <- ifelse(es_sano & !es_enfermo, "Sano",
         ifelse(es_enfermo & !es_sano, "Enfermo", NA))

# Si alguna muestra no cae claramente en ninguno de los dos grupos (o cae en
# los dos a la vez), es mejor detenerse y mostrar el texto antes que
# adivinar.
if (anyNA(grupo)) {
  print(texto_muestra[is.na(grupo)])
  stop("Hay muestras que no pude clasificar. Revisa los metadatos de arriba.")
}

grupo <- factor(grupo, levels = c("Sano", "Enfermo"))
print(table(grupo))     # para GSE9348 se espera: 12 sanas y 70 enfermas


# ---- 3. Preparar los datos ---------------------------------------------------

# Los microarreglos se suelen analizar en escala log2 porque la variabilidad de la
# expresión de genes crece junto con el nivel de expresión (un gen que se
# expresa mucho varía mucho más, en términos absolutos, que uno que apenas se
# expresa). El logaritmo aplana esa relación para que la varianza sea más
# parecida entre genes de expresión alta y baja.
if (max(expr, na.rm = TRUE) > 50) {
  cat("Valores en escala lineal -> aplicamos log2(x + 1)\n")
  expr <- log2(pmax(expr, 0) + 1)
}

# Sondas con NA en cualquier muestra se descartan: las pruebas de más abajo
# necesitan los 16 valores completos para poder comparar los dos grupos.
expr <- expr[complete.cases(expr), ]

# Para el ejemplo de la clase se usan pocas muestras (8 vs. 8),
# tomadas al azar de los 12 sanos y los 70 enfermos disponibles. Cambia
# `n_por_grupo` arriba para repetir el ejercicio con más muestras.
idx_sano    <- sample(which(grupo == "Sano"),    n_por_grupo)
idx_enfermo <- sample(which(grupo == "Enfermo"), n_por_grupo)
idx         <- c(idx_sano, idx_enfermo)

x   <- expr[, idx]
grp <- droplevels(grupo[idx])
cat("Muestras usadas en el ejemplo:\n"); print(table(grp))


# =============================================================================
#  PARTE A: el p-value de un gen mediante permutaciones
# =============================================================================

if (!sonda_ejemplo %in% rownames(x)) {
  # Si la sonda no está en este archivo, se usa la de mayor variabilidad
  # entre las 16 muestras, para no detener el script por una sonda ausente.
  sonda_ejemplo <- names(which.max(apply(x, 1, var)))
}
y <- x[sonda_ejemplo, ]
cat("\nGen/sonda de ejemplo:", sonda_ejemplo, "\n")

# --- Primero, ver los datos crudos----------
# Antes de hablar de distribuciones nulas conviene ver los 16 números tal
# cual son: ¿se ven separados los dos grupos a simple vista, o se traslapan?
# Nota: se usa method = "stack" , "stack"
# separa los puntos cercanos sin usar ningún número aleatorio.
stripchart(y ~ grp, vertical = TRUE, method = "stack", offset = 0.6,
           pch = 16, cex = 1.3,
           col = c("steelblue", "firebrick"),
           xlab = "", ylab = "Expresión (log2)",
           main = paste("Expresión de", sonda_ejemplo, "por paciente"))
# Una rayita horizontal en la media de cada grupo, para no tener que
# adivinarla a ojo entre los puntos.
segments(x0 = c(0.75, 1.75), x1 = c(1.25, 2.25),
         y0 = tapply(y, grp, mean), y1 = tapply(y, grp, mean),
         col = c("steelblue", "firebrick"), lwd = 3)
legend("topleft", bty = "n", lwd = 3, col = c("steelblue", "firebrick"),
       legend = c("Media Sano", "Media Enfermo"))

# Diferencia de medias (Enfermo - Sano). Se podría usar otra
# cosa (una mediana, un estadístico t), pero la diferencia de medias es la
# más fácil de explicar e interpretar: está en las mismas unidades (log2)
# que los datos originales.
dif_medias <- function(valores, etiquetas) {
  mean(valores[etiquetas == "Enfermo"]) - mean(valores[etiquetas == "Sano"])
}

d_obs <- dif_medias(y, grp)
cat("Diferencia observada (log2):", round(d_obs, 3), "\n")

# Distribución nula: si "Sano" y "Enfermo" fueran solo etiquetas sin ningún
# significado real, cualquier otra forma de repartir esas mismas 16
# etiquetas entre los mismos 16 valores sería igual de válida. Por eso se
# mezclan las etiquetas (sin tocar los valores) y se recalcula la diferencia
# de medias, B veces, para ver qué tan grande suele salir esa diferencia
# cuando no hay ninguna relación real entre el valor y la etiqueta.
d_perm <- replicate(B, dif_medias(y, sample(grp)))

# p-value bilateral: ¿qué proporción de esas permutaciones dio una diferencia
# al menos tan extrema (en cualquier dirección) como la que sí se observó?
b       <- sum(abs(d_perm) >= abs(d_obs))
# El "+1" en el numerador y en el denominador es
# importante: sin él, si ninguna permutación fuera tan extrema como la
# observada, se reportaría p = 0, lo cual sugeriría una certeza absoluta que
# en realidad no se tiene (solo se probaron B permutaciones, no todas las
# posibles). Con el "+1", el p-value más pequeño que se puede reportar es
# 1 / (B + 1), nunca exactamente cero.
p_perm  <- (b + 1) / (B + 1)
cat("Permutaciones al menos tan extremas:", b, "de", B, "\n")
cat("p-value por permutación:", signif(p_perm, 3), "\n")

hist(d_perm, breaks = 50, col = "grey85", border = "white",
     main = paste("Distribución nula por permutación -", sonda_ejemplo),
     sub  = paste0("B = ", B, " permutaciones  |  p = ", signif(p_perm, 3)),
     xlab = "Diferencia de medias con etiquetas mezcladas (log2)",
     xlim = range(c(d_perm, d_obs, -d_obs)))
abline(v = c(d_obs, -d_obs), col = "firebrick", lwd = 2, lty = c(1, 2))
legend("topright", bty = "n", col = "firebrick", lwd = 2, lty = c(1, 2),
       legend = c("Diferencia observada", "Mismo valor, signo contrario"))

# Comparación con la prueba t de Welch clásica: otra forma de llegar a un
# p-value para este mismo gen, sin permutar nada, asumiendo que los datos
# siguen una distribución normal. Suele dar un valor parecido al de
# permutación cuando ese supuesto es razonable.
print(t.test(y ~ grp))

# Pregunta: si cambias B a 99999, ¿el p-value necesariamente baja?


# =============================================================================
#  PARTE B: miles de pruebas y comparaciones múltiples
# =============================================================================

# Prueba t de Welch para todas las sondas a la vez. Está escrita "a mano"
# con operaciones de matrices (vectorizada) en vez de un bucle `for` con
# t.test() gen por gen: hace exactamente la misma cuenta, pero en segundos
# en vez de minutos, porque R hace la aritmética sobre toda la matriz de una
# sola vez.
welch_por_fila <- function(m, etiquetas) {
  a <- m[, etiquetas == "Enfermo", drop = FALSE]
  s <- m[, etiquetas == "Sano",    drop = FALSE]
  na <- ncol(a); ns <- ncol(s)
  va <- apply(a, 1, var) / na
  vs <- apply(s, 1, var) / ns
  t  <- (rowMeans(a) - rowMeans(s)) / sqrt(va + vs)
  # Grados de libertad de Welch-Satterthwaite: se usan en vez de los de la
  # t de Student clásica porque no se supone que los dos grupos tengan la
  # misma varianza.
  gl <- (va + vs)^2 / (va^2 / (na - 1) + vs^2 / (ns - 1))
  data.frame(log2FC = rowMeans(a) - rowMeans(s),
             t = t,
             p = 2 * pt(-abs(t), gl))
}

res <- welch_por_fila(x, grp)
res <- res[is.finite(res$p), ]           # quitamos sondas sin variabilidad
                                          # (varianza cero -> t indefinida)

# Comprobación: el resultado vectorizado debe coincidir con t.test() para la
# sonda de ejemplo
cat("p (vectorizado) =", signif(res[sonda_ejemplo, "p"], 4),
    " | p (t.test) =", signif(t.test(y ~ grp)$p.value, 4), "\n")

m <- nrow(res)
cat("\nNúmero de pruebas (sondas):", m, "\n")
cat("Sondas con p < 0.05 sin ajustar:", sum(res$p < alfa), "\n")
cat("Esperadas por azar si todas las H0 fueran ciertas: ~", round(alfa * m), "\n")

# --- Experimento clave: etiquetas mezcladas (no hay diferencia real) ----------
# Si se repite exactamente el mismo análisis pero con las etiquetas
# revueltas (así que "Sano" y "Enfermo" ya no significan nada), cualquier
# sonda que salga "significativa" aquí es, por construcción, un falso
# positivo: no puede haber una diferencia real entre dos grupos que se
# definieron al azar.
grp_mezclado <- sample(grp)
res_nulo <- welch_por_fila(x, grp_mezclado)
res_nulo <- res_nulo[is.finite(res_nulo$p), ]
cat("\ncon etiquetas mezcladas (sin señal real):\n")
cat("  Sondas con p < 0.05:", sum(res_nulo$p < alfa), "<- ¡falsos positivos!\n")
cat("  Sondas con FDR (BH) < 0.05:",
    sum(p.adjust(res_nulo$p, "BH") < alfa), "\n")

op <- par(mfrow = c(1, 2))
hist(res$p, breaks = 40, col = "steelblue", border = "white",
     main = "Etiquetas reales", xlab = "p-value")
hist(res_nulo$p, breaks = 40, col = "grey70", border = "white",
     main = "Etiquetas mezcladas", xlab = "p-value")
par(op)
# Cómo leer estos dos histogramas:
#  * Con etiquetas reales: un pico cerca de 0 (las sondas con señal real) y
#    un fondo que se esperaría plano (el resto, puro ruido de fondo).
#  * Con etiquetas mezcladas: histograma se esperaría plano, incluido cerca de 0; es lo que se espera si H0 es cierta en
#    todas las sondas.

# --- Corrección por comparaciones múltiples (R base: p.adjust) ----------------
# Bonferroni: multiplica cada p-value por el número de pruebas (en la
# práctica, equivale a exigir un umbral de alfa/m en vez de alfa). Controla
# la probabilidad de tener al menos un falso positivo en todo el
# experimento (FWER), por eso es tan estricto.
# Benjamini-Hochberg (FDR): en vez de eso, controla qué proporción de las
# sondas declaradas "significativas" se espera que sean falsos positivos.
# En ese sentido es menos estricto que Bonferroni y es el estándar en análisis
# exploratorios de expresión génica.
res$p_bonf <- p.adjust(res$p, method = "bonferroni")
res$FDR_BH <- p.adjust(res$p, method = "BH")

resumen <- data.frame(
  criterio = c("p < 0.05 sin ajustar", "Bonferroni < 0.05", "FDR (BH) < 0.05"),
  sondas   = c(sum(res$p < alfa), sum(res$p_bonf < alfa), sum(res$FDR_BH < alfa))
)
print(resumen)
cat("Umbral de Bonferroni equivalente: p <", signif(alfa / m, 2), "\n")

# El mismo resumen, pero como gráfica: ver cuánto se "encoge"
# la lista de sondas "significativas" según qué tan estricta sea la
# corrección
etiquetas_barras <- c("p < 0.05\nsin ajustar", "Bonferroni\n< 0.05", "FDR (BH)\n< 0.05")
barras <- barplot(resumen$sondas, names.arg = etiquetas_barras,
                   col = c("grey60", "firebrick", "steelblue"),
                   ylab = "Número de sondas",
                   main = paste0("De ", m, " sondas probadas, ¿cuántas 'sobreviven'?"))
text(barras, resumen$sondas, labels = resumen$sondas, pos = 3, xpd = TRUE)


# =============================================================================
#  PARTE C: leer resultados sin sobreinterpretarlos
# =============================================================================

res <- res[order(res$p), ]
tabla <- data.frame(sonda   = rownames(res),
                    log2FC  = round(res$log2FC, 2),
                    p       = signif(res$p, 3),
                    FDR_BH  = signif(res$FDR_BH, 3),
                    p_bonf  = signif(res$p_bonf, 3))
print(head(tabla, 10), row.names = FALSE)

# Ejemplo del "gen B" de la clase: pasa el umbral sin ajustar (p < 0.05) pero
# no sobrevive la corrección por FDR. Es justo el caso que muestra por qué
# no basta con mirar el p-value sin ajustar cuando se probaron miles de genes.
gen_B <- tabla[tabla$p < alfa & tabla$FDR_BH > alfa, ]
cat("\nSondas con p < 0.05 que NO pasan FDR < 0.05:", nrow(gen_B), "\n")
print(head(gen_B, 3), row.names = FALSE)

# Gráfico volcán: cada punto es una sonda. En el eje x, qué tan grande y en
# qué dirección es el cambio (log2FC); en el eje y, qué tan extremo es el
# p-value. Un punto arriba a la derecha o a la izquierda es a la vez grande
# Y estadísticamente firme; un punto arriba pero cerca del centro es
# estadísticamente firme pero de cambio chico.
plot(res$log2FC, -log10(res$p), pch = 16, cex = 0.5,
     col = ifelse(res$FDR_BH < alfa, "firebrick", "grey70"),
     xlab = "log2 fold change (Enfermo - Sano)", ylab = "-log10(p)",
     main = "Gráfico volcán (rojo: FDR < 0.05)")
abline(v = c(-1, 1), lty = 2, col = "grey40")
text(-1, par("usr")[4], "log2FC = -1", pos = 1, col = "grey40", cex = 0.8)
text( 1, par("usr")[4], "log2FC = 1",  pos = 1, col = "grey40", cex = 0.8)

# Se etiquetan con su nombre de sonda las 5 más significativas
top5 <- head(res[order(res$p), ], 5)
text(top5$log2FC, -log10(top5$p), labels = rownames(top5),
     pos = ifelse(top5$log2FC > 0, 4, 2), cex = 0.7, col = "black")

# Recordatorios para la discusión:
#  * FDR < 0.05 no significa que cada gen tenga 5 % de probabilidad de ser falso.
#  * Un p pequeño no mide el tamaño del efecto: mira log2FC.
#  * Un gen significativo no es todavía un biomarcador clínico.
#  * Este análisis es didáctico: un análisis diferencial formal usaría, por
#    ejemplo, limma (Bioconductor) que estabiliza la varianza de cada gen
#    compartiendo información entre todos los genes, en vez de estimarla
#    cada uno por su cuenta con pocas muestras, control de calidad de las
#    muestras, y covariables como edad o sexo que también vienen en los
#    metadatos de este mismo archivo.


# =============================================================================
#  Ejercicios para casa
# =============================================================================
# 1. Cambia `n_por_grupo` a 4 y luego a 12. ¿Qué pasa con el número de sondas
#    con FDR < 0.05? ¿Por qué?
# 2. Elige otra sonda (por ejemplo, la primera de `tabla`) y repite la Parte A.
# 3. Con etiquetas mezcladas, ¿cuántas sondas pasan Bonferroni? ¿Y BH?

sessionInfo()

