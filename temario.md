# Temario — Entendiendo los p-values

**Estadística aplicada a bioinformática**

### Objetivos de aprendizaje

Al finalizar el taller, las personas participantes podrán:

1. **Explicar qué representa un p-value** y cómo se interpreta en el contexto de una prueba de hipótesis.
2. **Describir intuitivamente una prueba de permutación** y comprender cómo permite evaluar si un resultado podría explicarse por el azar bajo una hipótesis nula.
3. **Explicar por qué las comparaciones múltiples aumentan el riesgo de obtener falsos positivos** cuando se realizan numerosas pruebas estadísticas.
4. **Distinguir entre p-values sin ajustar y ajustados**, e identificar el propósito de las correcciones de Bonferroni y Benjamini–Hochberg (FDR).

## Material

| Qué | Dónde |
|---|---|
| Script `entendiendo_pvalues.R` (solo R base) | GitHub |
| `diapositivas.pdf` (11 diapositivas) | GitHub |
| `GSE9348_series_matrix.txt.gz` | Carpeta de Drive del curso |

**Datos:** GEO [GSE9348](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9348). Mucosa de colon de 12 controles sanos y 70 tumores de cáncer colorrectal en estadio temprano, en la plataforma GPL570 (Affymetrix HG-U133 Plus 2.0). Son muestras independientes, así que la permutación libre de etiquetas es válida. En clase usamos 8 sanas y 8 enfermas tomadas al azar.

## Mapa del taller

| Bloque | Diapositivas | Script |
|---|---|---|
| Pregunta biológica | 1–2 | Secciones 0–3 (cargar datos) |
| Hipótesis y significado del p-value | 3–4 | — |
| Prueba de permutación | 5 | Parte A |
| Comparaciones múltiples | 6–7 | Parte B |
| Leer resultados sin sobreinterpretar | 8–9 | Parte C |
| Cierre y comprobación | 10–11 | Ejercicios |

---

## 1. Abrir con una pregunta biológica

*Diapositivas 1–2*

- La idea no es memorizar una definición, sino entender qué pregunta responde el p-value y qué errores se cometen al interpretarlo.
- **Caso:** comparamos mucosa de colon sana con tumor colorrectal (GSE9348). Un gen tiene una expresión media mayor en el tumor.
- ¿Basta esa diferencia? Todavía no: puede ser señal biológica, variabilidad natural entre muestras, efecto del tamaño de muestra o azar.
- Significancia estadística ≠ causalidad, ≠ biomarcador, ≠ importancia clínica.
- Presentar el dataset: qué es una *series matrix* de GEO (metadatos con `!` y una tabla de sondas × muestras), dónde está en Drive y que todo se hace con R base.
- **Hilo de la sesión:** hipótesis → p-value → permutación → ¿qué pasa con miles de pruebas?

> **Pregunta al grupo:** "Si analizamos 20 000 genes, ¿algunos podrían parecer significativos por azar aunque no haya diferencia real?" Acepta dos o tres respuestas y continúa.

**En R:** ejecutar las secciones 0–3. Mostrar `table(grupo)` (12 sanas y 70 enfermas) y explicar el submuestreo de 8 vs. 8.

## 2. Hipótesis y significado del p-value

*Diapositivas 3–4*

- **H₀:** no hay diferencia en la expresión media del gen entre sano y enfermo.
- **H₁:** sí hay diferencia.
- La prueba es bilateral si nos interesa cualquier dirección. La dirección se decide **antes** de ver los datos, nunca después para obtener un p más pequeño.
- **Estadístico:** resume la evidencia; aquí usamos la diferencia de medias.
- **Distribución nula:** los valores del estadístico que serían habituales si H₀ fuera cierta.
- **Definición:** probabilidad, *suponiendo H₀ y los supuestos del modelo*, de obtener un estadístico tan extremo como el observado o más.
- **Qué no es:**
  - la probabilidad de que H₀ sea verdadera;
  - la probabilidad de que el resultado "se deba al azar";
  - una medida del tamaño del efecto.
- p = 0.03 significa que, si H₀ fuera cierta, un resultado así de extremo aparecería alrededor del 3 % de las veces. No significa 97 % de probabilidad de que H₁ sea cierta.
- 0.05 es una convención, no una frontera mágica: p = 0.049 ≈ p = 0.051. Hay que reportar el valor exacto.
- Un p grande no demuestra ausencia de efecto: la ausencia de evidencia no es evidencia de ausencia.

> **Comprobación:** "¿p = 0.01 significa 99 % de probabilidad de que H₀ sea falsa?" → **No.** Es una probabilidad calculada suponiendo H₀.

## 3. Entender el p-value mediante permutaciones

*Diapositiva 5 · Parte A del script*

- Tenemos 8 sanas y 8 enfermas. Para una sonda calculamos la diferencia observada (Enfermo − Sano).
- **Idea:** si las etiquetas no significaran nada, ¿qué diferencias obtendríamos solo por cómo se reparten los valores?
- **Procedimiento:**
  1. calcular el estadístico observado;
  2. mezclar las etiquetas manteniendo fijos los valores;
  3. recalcular, repitiendo B = 9 999 veces;
  4. contar cuántas permutaciones son al menos tan extremas.
- **p Monte Carlo** = (b + 1) / (B + 1). Así nunca se reporta p = 0 solo porque ninguna permutación alcanzó el valor observado.
- **Supuesto:** las etiquetas deben ser intercambiables. Con datos pareados, bloques o medidas repetidas, las permutaciones tienen que respetar esa estructura. GSE9348 tiene muestras independientes, así que aquí se puede.
- **En R:** mostrar el histograma de la distribución nula con la línea de la diferencia observada. Comparar el p de permutación con el de `t.test()`; suelen ser parecidos.
- La sonda de ejemplo es `204259_at` (MMP7 en GPL570). Si no aparece, el script elige la sonda más variable.
- No hace falta saber programarlo desde cero: lo importante es poder explicar los 4 pasos.

> **Pregunta:** "Si aumentamos el número de permutaciones, ¿el p-value necesariamente baja?" → **No.** La estimación se vuelve más estable, pero no tiene por qué disminuir.

## 4. El problema de las comparaciones múltiples

*Diapositivas 6–7 · Parte B del script*

- En GPL570 hay unas 54 000 sondas, y probar cada una con α = 0.05 aumenta los falsos positivos.
- Si todas H₀ fueran ciertas, se esperarían unos 1 000 resultados con p < 0.05. Es un valor esperado, no exacto.
- **Demostración:**
  - con las etiquetas reales hay muchas sondas con p < 0.05;
  - con las etiquetas mezcladas, sin señal real, siguen saliendo cientos de p < 0.05;
  - con BH, casi ninguno sobrevive.
- Comparar los dos histogramas de p-values: con señal hay un pico cerca de 0; sin señal el histograma es aproximadamente plano.
- Como los genes están correlacionados, el número de falsos positivos varía de una mezcla a otra.
- **Bonferroni:**
  - controla la probabilidad de cometer al menos un falso positivo (FWER);
  - el umbral es α/m, por ejemplo 0.05/20 000 = 2.5 × 10⁻⁶;
  - es estricto y pierde efectos reales.
- **Benjamini–Hochberg:**
  - controla la FDR, la proporción esperada de falsos descubrimientos entre los resultados declarados significativos;
  - es lo habitual en análisis exploratorios de expresión de genes.
- FDR ≤ 0.05 no significa que cada gen tenga un 5 % de probabilidad de ser falso.
- En R basta con `p.adjust(p, method = "bonferroni")` y `p.adjust(p, method = "BH")`, sin paquetes.
- Ajustar no arregla un mal diseño, datos sesgados, confusión entre grupos ni supuestos incumplidos.

> **Pregunta:** "¿Cuál es más conservador, Bonferroni o BH?" → **Bonferroni**, porque controla el error familiar.

## 5. Leer resultados sin sobreinterpretarlos

*Diapositivas 8–9 · Parte C del script*

- **Tabla didáctica:**
  - Gen A: p = 0.001, FDR = 0.02 → pasa;
  - Gen B: p = 0.03, FDR = 0.18 → no pasa;
  - Gen C: p = 0.20, FDR = 0.40 → no pasa.
- El gen B muestra por qué no basta con el p sin ajustar. El script busca casos reales como este en GSE9348 ("sondas con p < 0.05 que no pasan FDR").
- Mostrar la tabla del script (log2FC, p, FDR, Bonferroni) y el gráfico volcán.
- Más allá de la significancia hay que revisar:
  - el tamaño y la dirección del cambio (log2FC);
  - la variabilidad y la calidad de las muestras;
  - las posibles variables de confusión;
  - la consistencia con otros estudios.
- Un efecto pequeño puede dar p muy bajo con muchas muestras; un efecto grande puede no ser significativo con pocas. **El p-value no mide el tamaño del efecto.**
- Pasar el umbral de FDR no convierte un gen en biomarcador clínico ni en diana terapéutica. Hacen falta validación independiente, experimentos funcionales y evaluación clínica.
- Este análisis es didáctico. Un análisis formal de microarreglos usaría, por ejemplo, limma, control de calidad y covariables.

> **Pregunta:** "Si un gen tiene FDR < 0.05, ¿ya es un buen biomarcador?" → **No.** Merece evaluación adicional; no es una validación.

## 6. Cierre

*Diapositivas 10–11*

1. El p-value se calcula suponiendo H₀ y unos supuestos.
2. No es la probabilidad de que H₀ sea cierta ni mide el tamaño del efecto.
3. Permutar etiquetas construye la distribución nula, siempre que el diseño lo permita.
4. Con miles de pruebas hay que corregir: Bonferroni (FWER) o BH (FDR), según el objetivo.
5. Significancia ≠ relevancia biológica, causalidad ni utilidad clínica.

**Pregunta para llevarse:** ¿qué hipótesis y supuestos hay detrás de este p-value, cuántas pruebas se hicieron y qué evidencia adicional necesito?

### Preguntas finales sugeridas

- ¿Qué significa un p-value de 0.03?
- ¿Por qué un p-value pequeño no demuestra que un gen cause una enfermedad?
- ¿Qué problema aparece al probar 20 000 genes?
- ¿Qué diferencia hay entre Bonferroni y BH/FDR?
- ¿Qué más revisarías antes de proponer un gen como biomarcador?

---

## Notas técnicas

- **Interpretación:** el p-value siempre es condicional a H₀ y a los supuestos de la prueba.
- **Permutación:** se usa una prueba bilateral con el valor absoluto de la diferencia de medias. El p se calcula como (b + 1)/(B + 1).
- **Intercambiabilidad:** GSE9348 tiene sanos y tumores de personas distintas. Con datos pareados (tumor y normal del mismo paciente) habría que permutar dentro de cada par.
- **BH:** controla la FDR bajo independencia y bajo ciertas formas de dependencia positiva, no bajo cualquier estructura de dependencia.
- **Escala:** el script detecta si los valores están en escala lineal y, en ese caso, aplica log2(x + 1).
- **Sondas, no genes:** en GPL570 varias sondas pueden corresponder al mismo gen. Para obtener símbolos de gen haría falta la anotación de la plataforma; se omite para no añadir paquetes.
- **Ejercicio 2 (cambiar la semilla):** al cambiar `semilla` se mueven dos cosas a la vez, no solo una: qué 8+8 pacientes se seleccionan y cómo se revuelven las etiquetas al permutar. Lo que se espera estable entre semillas es el fenómeno de fondo (con etiquetas mezcladas siempre aparecen varios cientos de sondas con p < 0.05 que casi todas desaparecen con FDR); el número exacto de falsos positivos sí varía un poco de una semilla a otra porque los genes están correlacionados. Lo que no se espera estable es el p-value de la sonda de ejemplo (MMP7): con una muestra de solo 8 contra 8, qué pacientes específicos tocaron sí importa, y el valor puede moverse bastante. La respuesta esperada del ejercicio no es "no cambió nada", sino notar que el patrón de falsos positivos por múltiples pruebas se repite de forma robusta, mientras que el resultado de un solo gen con pocas muestras puede ser inestable, lo cual conecta con por qué hacen falta más muestras o replicación para confiar en un hallazgo de un solo gen.
