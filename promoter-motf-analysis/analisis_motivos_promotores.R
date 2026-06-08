# =============================================================================
# ANÁLISIS DE MOTIVOS CIS-REGULADORES EN PROMOTORES
# Autor: Sergio Civera-Arroyo
# =============================================================================
#
# DESCRIPCIÓN GENERAL
# Este script realiza un análisis estadístico de motivos cis-reguladores
# (elementos reguladores en cis) en regiones promotoras de genes. El objetivo
# es determinar si ciertos motivos de secuencia aparecen con mayor frecuencia
# de la esperada por azar en los promotores analizados.
#
# El flujo general es:
#   1. Lectura y validación de datos de entrada (FASTA + definición de motivos)
#   2. Búsqueda de motivos en promotores reales
#   3. Simulación Monte Carlo para construir una distribución nula
#   4. Cálculo de enriquecimiento y p-valores estadísticamente correctos
#   5. Análisis posicional (absoluto y relativo)
#   6. Generación de gráficos y exportación de resultados
#
# =============================================================================
# ARCHIVOS NECESARIOS PARA EJECUTAR ESTE SCRIPT
# =============================================================================
#
# ---- ARCHIVO 1: FASTA DE PROMOTORES ----
#
# Formato: archivo FASTA estándar (.fa, .fasta, .txt)
# Contiene: una secuencia por promotor (región upstream de cada gen)
#
# Estructura esperada:
#   >NombreGen_o_ID_unico
#   ATCGATCGATCGATCG...
#   >OtroGen
#   GCTAGCTAGCTAGCTA...
#
# Requisitos:
#   - El encabezado (línea que empieza por ">") debe contener un identificador
#     único por promotor. El script tomará todo lo que hay tras ">" como nombre.
#   - Las secuencias pueden estar en una o varias líneas (el script las junta).
#   - Se aceptan bases A, T, C, G y también caracteres ambiguos IUPAC (N, R, Y,
#     etc.). Si hay caracteres ambiguos, quedan documentados en un warning.
#   - No se asume longitud fija: el script usa la longitud real de cada secuencia.
#
# ---- ARCHIVO 2: DEFINICIÓN DE MOTIVOS ----
#
# Formato: archivo CSV o TSV con cabecera
# Separador: coma (CSV) o tabulador (TSV). El script detecta automáticamente.
#
# Columnas mínimas obligatorias:
#   motif_id  → nombre o familia del motivo (p.ej. "ABRE", "DRE", "GBOX")
#   sequence  → secuencia concreta del motivo (p.ej. "ACGTG", "CACGTG")
#
# Columnas opcionales (se ignoran si no están presentes):
#   group     → agrupación funcional o familia más amplia
#   function  → función biológica descrita
#   source    → base de datos de origen (PLACE, JASPAR, etc.)
#
# Un mismo motif_id puede tener varias filas con distintas secuencias,
# representando variantes equivalentes del mismo motivo regulador.
#
# Ejemplo de contenido del archivo:
#   motif_id,sequence
#   ABRE,ACGTG
#   ABRE,CACGTG
#   DRE,GCCGAC
#   DRE,ACCGAC
#   GBOX,CACGTG
#
# Sobre caracteres IUPAC degenerados en motivos:
#   - Por defecto el script usa coincidencia exacta (A, T, C, G únicamente).
#   - Si se activa la opción usar_iupac = TRUE (ver sección de parámetros),
#     los códigos IUPAC en los motivos se expanden a expresiones regulares:
#       R=[AG], Y=[CT], S=[GC], W=[AT], K=[GT], M=[AC],
#       B=[CGT], D=[AGT], H=[ACT], V=[ACG], N=[ACGT]
#   - Esto permite buscar motivos degenerados como "RCCGAC" o "CACNNTG".
#
# ---- ARCHIVO 3 (OPCIONAL): CARPETA DE SALIDA ----
#
# Solo hay que definir la ruta en los parámetros de abajo.
# El script crea automáticamente la carpeta si no existe.
# En ella se guardarán todas las tablas CSV y gráficos PNG.
#
# =============================================================================
# PARÁMETROS CONFIGURABLES POR EL USUARIO
# (Solo es necesario modificar esta sección antes de ejecutar)
# =============================================================================

# ---- Rutas de entrada ----
ruta_fasta   <- "promotores.fa"       # Ruta al archivo FASTA con los promotores
ruta_motivos <- "motivos.tsv"         # Ruta al archivo CSV/TSV con los motivos
ruta_salida  <- "resultados_motivos"  # Carpeta donde se guardarán los resultados (se crea automáticamente)

# ---- Parámetros del análisis ----
n_sim         <- 2500   # Número de simulaciones Monte Carlo (mínimo 1000 recomendado)
semilla       <- 42      # Semilla aleatoria para reproducibilidad (cualquier entero)

# ---- Opciones de búsqueda ----
usar_iupac         <- TRUE  # TRUE: expande códigos IUPAC en motivos; FALSE: coincidencia exacta
buscar_rev_comp    <- TRUE   # TRUE: busca también en la hebra reverso-complementaria

# ---- Opciones de la simulación nula ---- 
# "mononucleotidica": usa frecuencias globales de A, T, C, G del FASTA real
# "dinucleotidica":   usa frecuencias de dinucleótidos para un modelo más realista
# "permutacion":      permuta aleatoriamente cada secuencia real (más conservador)
modelo_nulo <- "dinucleotidica"

# ---- Parámetros del análisis posicional ----
n_bins_relativo  <- 20   # Número de intervalos para el análisis posicional relativo (0-1)
tamano_bin_abs   <- 100  # Tamaño en pb de los bins para el análisis posicional absoluto

# ---- Ajuste de p-valores por múltiples comparaciones ----
metodo_ajuste <- "BH"   # "BH" = Benjamini-Hochberg (FDR); "bonferroni" también válido
umbral_padj   <- 0.05   # Umbral de significación tras ajuste

# =============================================================================
# LIBRERÍAS
# =============================================================================
# Comprobamos e instalamos automáticamente las librerías necesarias si faltan.

paquetes_necesarios <- c("ggplot2", "dplyr", "tidyr", "stringr", "readr",
                         "purrr", "scales", "ggridges")

for (pkg in paquetes_necesarios) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Instalando paquete:", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# =============================================================================
# INICIALIZACIÓN
# =============================================================================

# Fijamos la semilla aleatoria para garantizar reproducibilidad.
# Cualquier análisis con los mismos datos y la misma semilla producirá
# exactamente los mismos resultados.
set.seed(semilla)

# Registramos la fecha y hora de inicio del análisis
fecha_inicio <- Sys.time()
message("=== INICIO DEL ANÁLISIS: ", format(fecha_inicio, "%Y-%m-%d %H:%M:%S"), " ===")

# Creamos la carpeta de salida si no existe
if (!dir.exists(ruta_salida)) {
  dir.create(ruta_salida, recursive = TRUE)
  message("Carpeta de salida creada: ", ruta_salida)
}

# =============================================================================
# FUNCIONES AUXILIARES
# =============================================================================
# Definimos todas las funciones antes de ejecutar el pipeline principal.
# Esto facilita la lectura, el mantenimiento y la reutilización del código.

# -----------------------------------------------------------------------------
# FUNCIÓN: leer_fasta
# -----------------------------------------------------------------------------
# Qué hace: Lee un archivo FASTA y devuelve un data.frame con nombre,
#           secuencia y longitud real de cada promotor.
# Entrada:  ruta (character) → ruta al archivo FASTA
# Salida:   data.frame con columnas: gene, sequence, length
# Por qué:  El formato FASTA estándar requiere parseo propio porque las
#           secuencias pueden ocupar varias líneas y hay que unirlas.
# -----------------------------------------------------------------------------
leer_fasta <- function(ruta) {

  # Comprobamos que el archivo existe antes de intentar leerlo
  if (!file.exists(ruta)) {
    stop("ERROR: No se encuentra el archivo FASTA en la ruta: ", ruta,
         "\nRevisa la variable 'ruta_fasta' al inicio del script.")
  }

  lineas <- readLines(ruta, warn = FALSE)

  # Eliminamos líneas completamente vacías
  lineas <- lineas[nchar(trimws(lineas)) > 0]

  # Localizamos las líneas de encabezado (empiezan por ">")
  idx_header <- which(startsWith(lineas, ">"))

  if (length(idx_header) == 0) {
    stop("ERROR: El archivo FASTA no contiene encabezados válidos (líneas que empiezan por '>').")
  }

  # Extraemos nombres, secuencias y longitudes
  genes     <- character(length(idx_header))
  secuencias <- character(length(idx_header))

  for (i in seq_along(idx_header)) {
    # El nombre del gen es todo lo que hay tras ">"
    genes[i] <- trimws(sub("^>", "", lineas[idx_header[i]]))

    # Las líneas de secuencia van desde el siguiente encabezado hasta el anterior
    inicio_seq <- idx_header[i] + 1
    fin_seq    <- ifelse(i < length(idx_header), idx_header[i + 1] - 1, length(lineas))

    # Unimos las líneas de secuencia y pasamos a mayúsculas
    secuencias[i] <- paste(toupper(lineas[inicio_seq:fin_seq]), collapse = "")
  }

  df <- data.frame(
    gene     = genes,
    sequence = secuencias,
    length   = nchar(secuencias),
    stringsAsFactors = FALSE
  )

  return(df)
}

# -----------------------------------------------------------------------------
# FUNCIÓN: validar_promotores
# -----------------------------------------------------------------------------
# Qué hace: Comprueba la integridad del data.frame de promotores y avisa
#           de posibles problemas sin interrumpir el análisis.
# Entrada:  df_promotores (data.frame) con columnas gene, sequence, length
# Salida:   data.frame limpio (con los mismos datos, posiblemente filtrado)
# Por qué:  Es importante detectar problemas de calidad antes de análisis
#           para que los resultados sean confiables.
# -----------------------------------------------------------------------------
validar_promotores <- function(df_promotores) {

  n_total <- nrow(df_promotores)
  message("\n--- Validación de promotores ---")
  message("Promotores leídos: ", n_total)

  # Comprobamos nombres duplicados
  duplicados <- df_promotores$gene[duplicated(df_promotores$gene)]
  if (length(duplicados) > 0) {
    warning("Se encontraron nombres de promotor duplicados: ",
            paste(duplicados, collapse = ", "),
            "\nSe conserva solo la primera aparición de cada nombre.")
    df_promotores <- df_promotores[!duplicated(df_promotores$gene), ]
  }

  # Comprobamos secuencias vacías
  vacias <- df_promotores$gene[df_promotores$length == 0]
  if (length(vacias) > 0) {
    warning("Promotores con secuencia vacía eliminados: ",
            paste(vacias, collapse = ", "))
    df_promotores <- df_promotores[df_promotores$length > 0, ]
  }

  # Detectamos caracteres ambiguos (distintos de A, T, C, G)
  bases_validas <- "^[ATCGatcg]+$"
  tiene_ambiguos <- !grepl(bases_validas, df_promotores$sequence)
  if (any(tiene_ambiguos)) {
    warning("Los siguientes promotores contienen caracteres ambiguos (N, R, Y, etc.): ",
            paste(df_promotores$gene[tiene_ambiguos], collapse = ", "),
            "\nEstos caracteres se mantienen en la secuencia pero no coincidirán ",
            "con ningún motivo exacto.")
  }

  message("Promotores válidos tras validación: ", nrow(df_promotores))
  message("Longitud media: ", round(mean(df_promotores$length), 1), " pb")
  message("Longitud mínima: ", min(df_promotores$length), " pb")
  message("Longitud máxima: ", max(df_promotores$length), " pb")

  return(df_promotores)
}

# -----------------------------------------------------------------------------
# FUNCIÓN: leer_motivos
# -----------------------------------------------------------------------------
# Qué hace: Lee el archivo de definición de motivos (CSV o TSV), lo valida
#           y devuelve un data.frame limpio.
# Entrada:  ruta (character) → ruta al archivo de motivos
# Salida:   data.frame con columnas motif_id y sequence (al menos)
# Por qué:  Centralizar la lectura facilita cambiar el formato de entrada
#           sin tocar el resto del código.
# -----------------------------------------------------------------------------
leer_motivos <- function(ruta) {

  if (!file.exists(ruta)) {
    stop("ERROR: No se encuentra el archivo de motivos en: ", ruta,
         "\nRevisa la variable 'ruta_motivos' al inicio del script.")
  }

  # Detectamos si es CSV (coma) o TSV (tabulador) leyendo la primera línea
  primera_linea <- readLines(ruta, n = 1)
  sep_detectado <- ifelse(grepl("\t", primera_linea), "\t", ",")
  message("Separador detectado en archivo de motivos: '",
          ifelse(sep_detectado == "\t", "tabulador", "coma"), "'")

  df <- read.table(ruta, header = TRUE, sep = sep_detectado,
                   stringsAsFactors = FALSE, fill = TRUE, quote = "")

  # Comprobamos columnas obligatorias
  if (!"motif_id" %in% colnames(df)) {
    stop("ERROR: El archivo de motivos debe tener una columna llamada 'motif_id'.")
  }
  if (!"sequence" %in% colnames(df)) {
    stop("ERROR: El archivo de motivos debe tener una columna llamada 'sequence'.")
  }

  # Limpieza básica: mayúsculas, eliminar espacios, eliminar duplicados exactos
  df$motif_id  <- trimws(df$motif_id)
  df$sequence  <- toupper(trimws(df$sequence))
  df <- df[!duplicated(paste(df$motif_id, df$sequence)), ]
  df <- df[nchar(df$sequence) > 0, ]

  n_motivos  <- length(unique(df$motif_id))
  n_variantes <- nrow(df)
  message("Motivos únicos (familias): ", n_motivos)
  message("Variantes de secuencia totales: ", n_variantes)

  return(df)
}

# -----------------------------------------------------------------------------
# FUNCIÓN: iupac_a_regex
# -----------------------------------------------------------------------------
# Qué hace: Convierte una secuencia con códigos IUPAC degenerados en una
#           expresión regular de R que representa todas las bases posibles.
# Entrada:  seq (character) → secuencia con posibles códigos IUPAC
# Salida:   character → expresión regular equivalente
# Ejemplo:  "RCCGAC" → "[AG]CCGAC"
# Por qué:  Los motivos reguladores a veces se definen con degeneración.
#           Expandirlos a regex permite buscarlos sin enumerar todas las
#           variantes posibles a mano.
# -----------------------------------------------------------------------------
iupac_a_regex <- function(seq) {
  tabla_iupac <- c(
    "A" = "A", "T" = "T", "C" = "C", "G" = "G",
    "R" = "[AG]",  "Y" = "[CT]",  "S" = "[GC]",  "W" = "[AT]",
    "K" = "[GT]",  "M" = "[AC]",  "B" = "[CGT]", "D" = "[AGT]",
    "H" = "[ACT]", "V" = "[ACG]", "N" = "[ACGT]"
  )
  # Dividimos la secuencia en bases individuales y sustituimos cada una
  bases  <- strsplit(seq, "")[[1]]
  partes <- sapply(bases, function(b) {
    if (b %in% names(tabla_iupac)) tabla_iupac[b] else b
  })
  paste(partes, collapse = "")
}

# -----------------------------------------------------------------------------
# FUNCIÓN: rev_comp
# -----------------------------------------------------------------------------
# Qué hace: Calcula el reverso complementario de una secuencia de ADN.
# Entrada:  seq (character) → secuencia en mayúsculas
# Salida:   character → secuencia reverso-complementaria
# Por qué:  En un promotor de doble cadena, un motivo puede estar en
#           cualquiera de las dos hebras. Buscar en el reverso complementario
#           equivale a buscar en la hebra antisense.
# -----------------------------------------------------------------------------
rev_comp <- function(seq) {
  tabla_comp <- c("A"="T","T"="A","C"="G","G"="C","N"="N",
                  "R"="Y","Y"="R","S"="S","W"="W","K"="M","M"="K",
                  "B"="V","V"="B","D"="H","H"="D")
  bases <- strsplit(seq, "")[[1]]
  comp  <- sapply(bases, function(b) ifelse(b %in% names(tabla_comp), tabla_comp[b], b))
  paste(rev(comp), collapse = "")
}

# -----------------------------------------------------------------------------
# FUNCIÓN: buscar_motivo_en_secuencia
# -----------------------------------------------------------------------------
# Qué hace: Busca todas las ocurrencias de un motivo (con sus variantes) en
#           una secuencia promotora. Devuelve una tabla con posición y variante.
# Entrada:
#   promotor_name  → nombre del gen/promotor
#   seq_promotor   → secuencia del promotor en mayúsculas
#   motif_id       → nombre del motivo
#   variantes      → vector de secuencias variantes del motivo
#   usar_iupac     → (lógico) si TRUE, expande códigos IUPAC a regex
#   buscar_rev_comp → (lógico) si TRUE, busca también en la hebra complementaria
# Salida:  data.frame con columnas: gene, motif_id, variante, posicion_inicio,
#          posicion_fin, longitud_motivo, hebra
# Por qué: Esta función es el núcleo del análisis. Permite detectar y localizar
#          con precisión cada ocurrencia de cada motivo en cada promotor.
# -----------------------------------------------------------------------------
buscar_motivo_en_secuencia <- function(promotor_name, seq_promotor, motif_id,
                                        variantes, usar_iupac, buscar_rev_comp) {

  resultados <- list()

  for (var in variantes) {
    # Construimos el patrón de búsqueda: regex o coincidencia exacta
    patron <- if (usar_iupac) iupac_a_regex(var) else var

    # ---- Búsqueda en la hebra directa (sense) ----
    # gregexpr devuelve todas las posiciones de inicio de cada ocurrencia
    # Nota: usamos perl=TRUE para permitir lookahead y manejar solapamientos
    m <- gregexpr(patron, seq_promotor, perl = TRUE)[[1]]

    if (m[1] != -1) {
      # Para cada ocurrencia registramos su posición de inicio y fin
      for (pos in m) {
        resultados <- append(resultados, list(data.frame(
          gene            = promotor_name,
          motif_id        = motif_id,
          variante        = var,
          posicion_inicio = pos,
          posicion_fin    = pos + nchar(var) - 1,
          longitud_motivo = nchar(var),
          hebra           = "+",
          stringsAsFactors = FALSE
        )))
      }
    }

    # ---- Búsqueda en la hebra reverso-complementaria (antisense) ----
    if (buscar_rev_comp) {
      rc_seq <- rev_comp(seq_promotor)
      m_rc   <- gregexpr(patron, rc_seq, perl = TRUE)[[1]]

      if (m_rc[1] != -1) {
        len_prom <- nchar(seq_promotor)
        for (pos_rc in m_rc) {
          # Convertimos la posición en la RC a coordenada en la hebra directa
          pos_directa <- len_prom - pos_rc - nchar(var) + 2
          resultados <- append(resultados, list(data.frame(
            gene            = promotor_name,
            motif_id        = motif_id,
            variante        = var,
            posicion_inicio = pos_directa,
            posicion_fin    = pos_directa + nchar(var) - 1,
            longitud_motivo = nchar(var),
            hebra           = "-",
            stringsAsFactors = FALSE
          )))
        }
      }
    }
  }

  # Si no se encontró ninguna ocurrencia, devolvemos un data.frame vacío
  if (length(resultados) == 0) {
    return(data.frame(gene=character(), motif_id=character(), variante=character(),
                      posicion_inicio=integer(), posicion_fin=integer(),
                      longitud_motivo=integer(), hebra=character(),
                      stringsAsFactors=FALSE))
  }

  do.call(rbind, resultados)
}

# -----------------------------------------------------------------------------
# FUNCIÓN: buscar_todos_motivos
# -----------------------------------------------------------------------------
# Qué hace: Aplica la búsqueda de todos los motivos en todos los promotores
#           y devuelve la tabla completa de ocurrencias reales.
# Entrada:
#   df_promotores → data.frame con gene, sequence, length
#   df_motivos    → data.frame con motif_id, sequence
#   usar_iupac, buscar_rev_comp → parámetros heredados
# Salida:  data.frame con todas las ocurrencias encontradas
# Por qué: Centralizar la doble iteración (motivos × promotores) en una
#          función hace más fácil medir el tiempo y controlar el progreso.
# -----------------------------------------------------------------------------
buscar_todos_motivos <- function(df_promotores, df_motivos, usar_iupac, buscar_rev_comp) {

  # Obtenemos la lista de motif_ids únicos para iterar sobre ellos
  ids_motivos <- unique(df_motivos$motif_id)
  n_motivos   <- length(ids_motivos)
  n_prom      <- nrow(df_promotores)

  message("\n--- Búsqueda de motivos en promotores reales ---")
  message("Motivos a buscar: ", n_motivos, " | Promotores: ", n_prom)

  todas_ocurrencias <- list()

  for (i in seq_along(ids_motivos)) {
    mid <- ids_motivos[i]
    # Obtenemos todas las variantes de secuencia para este motif_id
    variantes <- df_motivos$sequence[df_motivos$motif_id == mid]

    message(sprintf("  [%d/%d] Buscando motivo: %s (%d variantes)",
                    i, n_motivos, mid, length(variantes)))

    # Comprobamos si alguna variante es más larga que algún promotor
    max_var <- max(nchar(variantes))
    prom_cortos <- df_promotores$gene[df_promotores$length < max_var]
    if (length(prom_cortos) > 0) {
      warning("El motivo '", mid, "' (longitud máxima: ", max_var, " pb) ",
              "es más largo que ", length(prom_cortos), " promotor(es). ",
              "Esos promotores se omiten para este motivo.")
    }

    # Buscamos en cada promotor
    for (j in seq_len(n_prom)) {
      prom <- df_promotores[j, ]

      # Omitimos promotores más cortos que la variante más larga
      if (prom$length < max_var) next

      ocurr <- buscar_motivo_en_secuencia(
        promotor_name   = prom$gene,
        seq_promotor    = prom$sequence,
        motif_id        = mid,
        variantes       = variantes,
        usar_iupac      = usar_iupac,
        buscar_rev_comp = buscar_rev_comp
      )

      if (nrow(ocurr) > 0) {
        todas_ocurrencias <- append(todas_ocurrencias, list(ocurr))
      }
    }
  }

  if (length(todas_ocurrencias) == 0) {
    warning("No se encontró ningún motivo en ningún promotor. ",
            "Revisa las rutas, los nombres y la configuración de búsqueda.")
    return(data.frame(gene=character(), motif_id=character(), variante=character(),
                      posicion_inicio=integer(), posicion_fin=integer(),
                      longitud_motivo=integer(), hebra=character()))
  }

  do.call(rbind, todas_ocurrencias)
}

# -----------------------------------------------------------------------------
# FUNCIÓN: calcular_metricas_observadas
# -----------------------------------------------------------------------------
# Qué hace: A partir de la tabla de ocurrencias reales, calcula para cada
#           motif_id el número de promotores positivos, el total de ocurrencias,
#           la frecuencia relativa y la densidad de ocurrencias.
# Entrada:
#   df_ocurrencias → data.frame de ocurrencias (salida de buscar_todos_motivos)
#   df_promotores  → data.frame con la información de todos los promotores
# Salida:  data.frame resumen con una fila por motif_id
# Por qué: Estas métricas observadas son el punto de comparación contra las
#          simulaciones. Sin ellas no podemos calcular enriquecimiento.
# -----------------------------------------------------------------------------
calcular_metricas_observadas <- function(df_ocurrencias, df_promotores) {

  n_total_prom   <- nrow(df_promotores)
  longitud_total <- sum(df_promotores$length)

  ids_motivos <- unique(df_ocurrencias$motif_id)

  resumen <- lapply(ids_motivos, function(mid) {
    sub <- df_ocurrencias[df_ocurrencias$motif_id == mid, ]
    prom_positivos   <- length(unique(sub$gene))
    total_ocurrencias <- nrow(sub)
    frecuencia_relativa <- prom_positivos / n_total_prom
    densidad_ocurrencias <- total_ocurrencias / longitud_total

    data.frame(
      motif_id              = mid,
      observed_promoters    = prom_positivos,
      observed_occurrences  = total_ocurrencias,
      freq_relativa         = frecuencia_relativa,
      densidad_ocurrencias  = densidad_ocurrencias,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, resumen)
}

# =============================================================================
# FUNCIONES PARA LA SIMULACIÓN MONTE CARLO
# =============================================================================

# -----------------------------------------------------------------------------
# FUNCIÓN: calcular_frecuencias_bases
# -----------------------------------------------------------------------------
# Qué hace: Calcula las frecuencias de bases (A, T, C, G) o dinucleótidos
#           en el conjunto de promotores reales para usar como modelo nulo.
# Entrada:
#   df_promotores → data.frame con columna sequence
#   modelo        → "mononucleotidica" o "dinucleotidica"
# Salida:  named vector de probabilidades
# Por qué: El modelo nulo debe ser fiel a la composición real del genoma.
#          Un promotor de AT-rich no debe compararse con una esperanza
#          calculada con 25% de cada base.
# -----------------------------------------------------------------------------
calcular_frecuencias_bases <- function(df_promotores, modelo = "mononucleotidica") {

  todas_seqs <- paste(df_promotores$sequence, collapse = "")

  if (modelo == "mononucleotidica") {
    bases <- c("A", "T", "C", "G")
    conteo <- sapply(bases, function(b) {
      nchar(todas_seqs) - nchar(gsub(b, "", todas_seqs, fixed = TRUE))
    })
    probs <- conteo / sum(conteo)
    names(probs) <- bases
    message("Frecuencias mononucleotídicas: ",
            paste(names(probs), round(probs, 3), sep = "=", collapse = " | "))
    return(probs)

  } else if (modelo == "dinucleotidica") {
    # Extraemos todos los dinucleótidos solapantes de todas las secuencias
    dinucs <- character(0)
    for (seq in df_promotores$sequence) {
      n <- nchar(seq)
      if (n >= 2) {
        dinucs <- c(dinucs, substring(seq, 1:(n-1), 2:n))
      }
    }
    # Nos quedamos solo con los dinucleótidos puros (sin ambiguos)
    dinucs <- dinucs[grepl("^[ATCG]{2}$", dinucs)]
    tabla  <- table(dinucs)
    probs  <- prop.table(tabla)
    message("Modelo dinucleotídico: ", length(probs), " dinucleótidos únicos calculados.")
    return(probs)

  } else {
    stop("modelo_nulo desconocido: use 'mononucleotidica', 'dinucleotidica' o 'permutacion'")
  }
}

# -----------------------------------------------------------------------------
# FUNCIÓN: generar_promotor_aleatorio
# -----------------------------------------------------------------------------
# Qué hace: Genera una secuencia de ADN aleatoria de la longitud indicada,
#           usando el modelo de composición de bases especificado.
# Entrada:
#   longitud → entero, longitud de la secuencia a generar
#   probs    → named vector de probabilidades de bases o dinucleótidos
#   modelo   → "mononucleotidica", "dinucleotidica" o "permutacion"
#   seq_real → secuencia original (solo necesaria si modelo == "permutacion")
# Salida:  character → secuencia aleatoria de la longitud indicada
# Por qué: Generamos promotores sintéticos que preservan la longitud real de
#          cada promotor observado. Esto construye una hipótesis nula más
#          realista que usar una longitud fija para todas las secuencias.
# -----------------------------------------------------------------------------
generar_promotor_aleatorio <- function(longitud, probs, modelo, seq_real = NULL) {

  if (modelo == "mononucleotidica") {
    # Muestreamos bases independientemente con las probabilidades globales
    bases <- sample(names(probs), size = longitud, replace = TRUE, prob = probs)
    return(paste(bases, collapse = ""))

  } else if (modelo == "dinucleotidica") {
    # Modelo de cadena de Markov de orden 1 basado en frecuencias de dinucleótidos
    # Calculamos las probabilidades de transición: P(base_j | base_i)
    bases  <- c("A", "T", "C", "G")
    trans  <- matrix(0, nrow = 4, ncol = 4, dimnames = list(bases, bases))

    for (b1 in bases) {
      for (b2 in bases) {
        dinuc <- paste0(b1, b2)
        if (dinuc %in% names(probs)) {
          trans[b1, b2] <- probs[dinuc]
        }
      }
      if (sum(trans[b1, ]) > 0) {
        trans[b1, ] <- trans[b1, ] / sum(trans[b1, ])
      } else {
        trans[b1, ] <- 0.25  # fallback uniforme
      }
    }

    # Generamos la secuencia base a base usando la cadena de Markov
    seq_gen <- character(longitud)
    # Base inicial aleatoria según la distribución marginal
    marginal <- colSums(trans) / sum(colSums(trans))
    seq_gen[1] <- sample(bases, 1, prob = marginal)
    for (k in 2:longitud) {
      seq_gen[k] <- sample(bases, 1, prob = trans[seq_gen[k-1], ])
    }
    return(paste(seq_gen, collapse = ""))

  } else if (modelo == "permutacion") {
    # Permutamos aleatoriamente las bases de la secuencia real.
    # Este es el modelo más conservador: preserva exactamente la composición
    # de cada promotor individual.
    if (is.null(seq_real)) stop("Para el modelo 'permutacion' se necesita seq_real.")
    bases <- strsplit(seq_real, "")[[1]]
    return(paste(sample(bases), collapse = ""))

  } else {
    stop("modelo_nulo desconocido.")
  }
}

# -----------------------------------------------------------------------------
# FUNCIÓN: simular_conjunto_completo
# -----------------------------------------------------------------------------
# Qué hace: Genera un conjunto completo de promotores simulados, uno por cada
#           promotor real, preservando exactamente las longitudes individuales.
# Entrada:
#   df_promotores → data.frame con gene, sequence, length
#   probs         → frecuencias de bases calculadas por calcular_frecuencias_bases
#   modelo        → string del modelo nulo
# Salida:  data.frame con gene, sequence, length (simulado)
# Por qué: Esta función garantiza que la hipótesis nula tenga la misma
#          estructura que los datos reales (mismo número de promotores,
#          mismas longitudes), haciendo la comparación estadísticamente válida.
# -----------------------------------------------------------------------------
simular_conjunto_completo <- function(df_promotores, probs, modelo) {
  seqs_sim <- mapply(function(longitud, seq_real) {
    generar_promotor_aleatorio(longitud, probs, modelo, seq_real)
  }, df_promotores$length, df_promotores$sequence)

  df_sim <- df_promotores
  df_sim$sequence <- seqs_sim
  return(df_sim)
}

# -----------------------------------------------------------------------------
# FUNCIÓN: ejecutar_monte_carlo_motivo
# -----------------------------------------------------------------------------
# Qué hace: Realiza n_sim simulaciones para UN motif_id concreto y devuelve
#           la distribución nula de promotores positivos y ocurrencias totales.
# Entrada:
#   mid           → character, nombre del motif_id
#   variantes     → vector de secuencias de variantes del motivo
#   df_promotores → data.frame con gene, sequence, length
#   probs         → frecuencias de bases
#   modelo        → modelo nulo
#   n_sim         → número de simulaciones
#   usar_iupac, buscar_rev_comp → opciones de búsqueda
# Salida:  list con vectores sim_promoters y sim_occurrences (longitud n_sim)
# Por qué: Repetir la búsqueda en conjuntos aleatorios nos da la distribución
#          de lo que esperaríamos por azar, que es nuestra hipótesis nula.
#          Cuanto mayor sea n_sim, más precisa es la estimación del p-valor.
# -----------------------------------------------------------------------------
ejecutar_monte_carlo_motivo <- function(mid, variantes, df_promotores, probs,
                                        modelo, n_sim, usar_iupac, buscar_rev_comp) {

  sim_promoters    <- integer(n_sim)
  sim_occurrences  <- integer(n_sim)
  max_var          <- max(nchar(variantes))

  for (s in seq_len(n_sim)) {
    # Generamos un conjunto completo de promotores aleatorios
    df_sim <- simular_conjunto_completo(df_promotores, probs, modelo)

    # Contamos promotores con al menos una ocurrencia y total de ocurrencias
    total_prom_pos <- 0
    total_ocurr    <- 0

    for (j in seq_len(nrow(df_sim))) {
      prom <- df_sim[j, ]
      if (prom$length < max_var) next

      ocurr <- buscar_motivo_en_secuencia(
        promotor_name   = prom$gene,
        seq_promotor    = prom$sequence,
        motif_id        = mid,
        variantes       = variantes,
        usar_iupac      = usar_iupac,
        buscar_rev_comp = buscar_rev_comp
      )

      if (nrow(ocurr) > 0) {
        total_prom_pos <- total_prom_pos + 1
        total_ocurr    <- total_ocurr + nrow(ocurr)
      }
    }

    sim_promoters[s]   <- total_prom_pos
    sim_occurrences[s] <- total_ocurr
  }

  return(list(
    sim_promoters   = sim_promoters,
    sim_occurrences = sim_occurrences
  ))
}

# -----------------------------------------------------------------------------
# FUNCIÓN: calcular_estadisticos_simulacion
# -----------------------------------------------------------------------------
# Qué hace: A partir de los vectores de simulación y los valores observados,
#           calcula todos los estadísticos de enriquecimiento y p-valores.
# Entrada:
#   sim_res       → list con sim_promoters y sim_occurrences
#   obs_promoters → entero, número observado de promotores positivos
#   obs_ocurr     → entero, número observado de ocurrencias totales
# Salida:  data.frame con una fila y todos los estadísticos calculados
# Por qué: El p-valor Monte Carlo se define como la proporción de simulaciones
#          en las que el valor simulado es mayor o igual que el observado.
#          Esto es equivalente a preguntar: ¿qué fracción de mundos aleatorios
#          produce un resultado tan extremo o más que el que observamos?
# -----------------------------------------------------------------------------
calcular_estadisticos_simulacion <- function(sim_res, obs_promoters, obs_ocurr, mid) {

  sp <- sim_res$sim_promoters
  so <- sim_res$sim_occurrences

  # Estadísticos de la distribución nula para promotores positivos
  exp_prom_mean <- mean(sp)
  exp_prom_med  <- median(sp)
  exp_prom_q025 <- quantile(sp, 0.025)
  exp_prom_q975 <- quantile(sp, 0.975)

  # Estadísticos de la distribución nula para ocurrencias totales
  exp_ocurr_mean <- mean(so)
  exp_ocurr_med  <- median(so)
  exp_ocurr_q025 <- quantile(so, 0.025)
  exp_ocurr_q975 <- quantile(so, 0.975)

  # Enriquecimiento:
  # > 1 indica más ocurrencias de las esperadas; < 1 indica depleción.
  # Evitamos dividir por cero sumando un pequeño epsilon.
  enrich_prom  <- obs_promoters / max(exp_prom_mean, 0.001)
  enrich_ocurr <- obs_ocurr     / max(exp_ocurr_mean, 0.001)

  # P-valor Monte Carlo unilateral superior (enriquecimiento):
  # Proporción de simulaciones donde el valor simulado >= valor observado.
  # Un p-valor pequeño indica que el motivo aparece más de lo esperado por azar.
  # Sumamos 1 en numerador y denominador (corrección de Davison-Hinkley)
  # para evitar p-valores exactamente 0, que serían demasiado optimistas.
  pval_prom  <- (sum(sp >= obs_promoters) + 1) / (length(sp) + 1)
  pval_ocurr <- (sum(so >= obs_ocurr) + 1)     / (length(so) + 1)

  data.frame(
    motif_id                = mid,
    expected_promoters_mean = exp_prom_mean,
    expected_promoters_med  = exp_prom_med,
    expected_promoters_q025 = as.numeric(exp_prom_q025),
    expected_promoters_q975 = as.numeric(exp_prom_q975),
    expected_occurrences_mean = exp_ocurr_mean,
    expected_occurrences_med  = exp_ocurr_med,
    expected_occurrences_q025 = as.numeric(exp_ocurr_q025),
    expected_occurrences_q975 = as.numeric(exp_ocurr_q975),
    enrichment_promoters    = enrich_prom,
    enrichment_occurrences  = enrich_ocurr,
    pvalue_promoters        = pval_prom,
    pvalue_occurrences      = pval_ocurr,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# FUNCIONES PARA EL ANÁLISIS POSICIONAL
# =============================================================================

# -----------------------------------------------------------------------------
# FUNCIÓN: calcular_posiciones_relativas
# -----------------------------------------------------------------------------
# Qué hace: Añade a la tabla de ocurrencias reales la posición relativa
#           (entre 0 y 1) de cada ocurrencia dentro de su promotor.
# Entrada:
#   df_ocurrencias → data.frame de ocurrencias reales
#   df_promotores  → data.frame con gene y length
# Salida:  df_ocurrencias con columna adicional pos_relativa
# Por qué: Como los promotores tienen distintas longitudes, comparar
#          posiciones absolutas sería engañoso. La posición relativa
#          permite visualizar si los motivos se concentran cerca del
#          inicio (0), del final (1), o en alguna zona particular.
# -----------------------------------------------------------------------------
calcular_posiciones_relativas <- function(df_ocurrencias, df_promotores) {

  # Añadimos la longitud de cada promotor al data.frame de ocurrencias
  df_merged <- merge(df_ocurrencias,
                     df_promotores[, c("gene", "length")],
                     by = "gene", all.x = TRUE)

  # Calculamos la posición relativa: inicio / longitud total del promotor
  # Da un valor de 0 (inicio) a ~1 (final del promotor)
  df_merged$pos_relativa <- df_merged$posicion_inicio / df_merged$length

  return(df_merged)
}

# -----------------------------------------------------------------------------
# FUNCIÓN: calcular_bins_posicionales
# -----------------------------------------------------------------------------
# Qué hace: Agrupa las ocurrencias en bins posicionales (absolutos o relativos)
#           y calcula el número y proporción de ocurrencias por bin.
# Entrada:
#   df_pos     → data.frame con ocurrencias y sus posiciones
#   tipo       → "relativo" o "absoluto"
#   n_bins_rel → número de bins para el análisis relativo
#   tam_bin_abs → tamaño en pb de cada bin absoluto
# Salida:  data.frame con bin, motif_id, n_ocurrencias, proporcion
# Por qué: Los bins resumen la distribución posicional de forma cuantitativa,
#          permitiendo detectar zonas de enriquecimiento posicional.
# -----------------------------------------------------------------------------
calcular_bins_posicionales <- function(df_pos, tipo = "relativo",
                                        n_bins_rel = 20, tam_bin_abs = 100) {

  ids_motivos <- unique(df_pos$motif_id)
  resultados  <- list()

  for (mid in ids_motivos) {
    sub <- df_pos[df_pos$motif_id == mid, ]
    if (nrow(sub) == 0) next

    if (tipo == "relativo") {
      # Dividimos el rango [0, 1] en n_bins_rel intervalos iguales
      breaks <- seq(0, 1, length.out = n_bins_rel + 1)
      sub$bin <- cut(sub$pos_relativa, breaks = breaks,
                     labels = paste0("[", round(breaks[-length(breaks)], 2),
                                     "-", round(breaks[-1], 2), ")"),
                     include.lowest = TRUE)
    } else {
      # Bins de tamaño fijo en pb (posición absoluta)
      max_pos <- max(sub$posicion_inicio)
      breaks  <- seq(0, max_pos + tam_bin_abs, by = tam_bin_abs)
      sub$bin <- cut(sub$posicion_inicio, breaks = breaks, include.lowest = TRUE)
    }

    tabla_bin <- as.data.frame(table(sub$bin))
    colnames(tabla_bin) <- c("bin", "n_ocurrencias")
    tabla_bin$motif_id  <- mid
    tabla_bin$proporcion <- tabla_bin$n_ocurrencias / sum(tabla_bin$n_ocurrencias)

    resultados <- append(resultados, list(tabla_bin))
  }

  if (length(resultados) == 0) return(data.frame())
  do.call(rbind, resultados)
}

# =============================================================================
# FUNCIONES DE VISUALIZACIÓN
# =============================================================================

# -----------------------------------------------------------------------------
# FUNCIÓN: guardar_grafico
# -----------------------------------------------------------------------------
# Pequeña función auxiliar para guardar gráficos con configuración consistente.
# -----------------------------------------------------------------------------
guardar_grafico <- function(p, nombre_archivo, ancho = 10, alto = 6) {
  ruta_completa <- file.path(ruta_salida, nombre_archivo)
  ggsave(ruta_completa, plot = p, width = ancho, height = alto, dpi = 150)
  message("Gráfico guardado: ", ruta_completa)
}

# -----------------------------------------------------------------------------
# FUNCIÓN: grafico_enriquecimiento
# -----------------------------------------------------------------------------
# Qué hace: Gráfico de barras del enriquecimiento de promotores para cada
#           motivo, coloreado por significación estadística (FDR).
# Entrada:  df_final → tabla resumen final con todos los estadísticos
# Salida:   objeto ggplot + archivo PNG guardado en ruta_salida
# Por qué:  Permite identificar de un vistazo qué motivos están sobre- o
#           sub-representados y cuáles son estadísticamente significativos.
# -----------------------------------------------------------------------------
grafico_enriquecimiento <- function(df_final) {

  df_plot <- df_final %>%
    arrange(desc(enrichment_promoters)) %>%
    mutate(
      motif_id    = factor(motif_id, levels = motif_id),
      significativo = padj_promoters < umbral_padj
    )

  p <- ggplot(df_plot, aes(x = motif_id, y = enrichment_promoters,
                            fill = significativo)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_fill_manual(values = c("TRUE" = "#E63946", "FALSE" = "#A8DADC"),
                      labels = c("TRUE" = paste0("padj < ", umbral_padj),
                                 "FALSE" = "No significativo")) +
    labs(
      title = "Enriquecimiento de motivos en promotores",
      subtitle = paste("Enriquecimiento = observado / esperado (simulación Monte Carlo, n =",
                       n_sim, "iteraciones)"),
      x = "Motivo",
      y = "Enriquecimiento (obs/esp)",
      fill = "Significación"
    ) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  guardar_grafico(p, "enrichment_motivos.png", ancho = 12, alto = 6)
  return(invisible(p))
}

# -----------------------------------------------------------------------------
# FUNCIÓN: grafico_ocurrencias_observadas
# -----------------------------------------------------------------------------
# Qué hace: Gráfico de barras con el número total de ocurrencias observadas
#           para cada motivo en los promotores reales.
# Entrada:  df_final → tabla resumen final
# Salida:   objeto ggplot + archivo PNG
# Por qué:  Da una perspectiva rápida de cuáles son los motivos más frecuentes
#           en términos absolutos en el conjunto de promotores.
# -----------------------------------------------------------------------------
grafico_ocurrencias_observadas <- function(df_final) {

  df_plot <- df_final %>%
    arrange(desc(observed_occurrences)) %>%
    mutate(motif_id = factor(motif_id, levels = motif_id))

  p <- ggplot(df_plot, aes(x = motif_id, y = observed_occurrences)) +
    geom_bar(stat = "identity", fill = "#457B9D", color = "black", linewidth = 0.3) +
    labs(
      title = "Ocurrencias totales de cada motivo en promotores reales",
      x = "Motivo",
      y = "Número total de ocurrencias"
    ) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  guardar_grafico(p, "ocurrencias_motivos.png", ancho = 12, alto = 6)
  return(invisible(p))
}

# -----------------------------------------------------------------------------
# FUNCIÓN: grafico_distribucion_nula
# -----------------------------------------------------------------------------
# Qué hace: Para cada motivo (o los más relevantes), genera un histograma
#           de la distribución nula simulada con una línea vertical marcando
#           el valor observado.
# Entrada:
#   mid         → motif_id del motivo a graficar
#   sim_res     → list con sim_promoters y sim_occurrences
#   obs_prom    → valor observado de promotores positivos
# Salida:   archivo PNG guardado
# Por qué:  Visualizar la distribución nula junto al valor observado es la
#           forma más intuitiva de entender el p-valor Monte Carlo.
#           Si el valor observado cae claramente en la cola derecha, el
#           motivo está significativamente enriquecido.
# -----------------------------------------------------------------------------
grafico_distribucion_nula <- function(mid, sim_res, obs_prom) {

  df_plot <- data.frame(sim_promoters = sim_res$sim_promoters)

  p <- ggplot(df_plot, aes(x = sim_promoters)) +
    geom_histogram(bins = 30, fill = "#A8DADC", color = "white", linewidth = 0.3) +
    geom_vline(xintercept = obs_prom, color = "#E63946", linewidth = 1.2,
               linetype = "solid") +
    annotate("text", x = obs_prom, y = Inf, label = paste("Observado:", obs_prom),
             color = "#E63946", hjust = -0.1, vjust = 1.5, size = 4) +
    labs(
      title = paste("Distribución nula simulada:", mid),
      subtitle = paste("Número de promotores con el motivo en cada simulación (n =", n_sim, ")"),
      x = "Promotores con el motivo (simulación)",
      y = "Frecuencia"
    ) +
    theme_bw(base_size = 13)

  guardar_grafico(p, paste0("distribucion_nula_", mid, ".png"), ancho = 8, alto = 5)
  return(invisible(p))
}

# -----------------------------------------------------------------------------
# FUNCIÓN: grafico_distribucion_posicional
# -----------------------------------------------------------------------------
# Qué hace: Histograma de la distribución de posiciones relativas de las
#           ocurrencias de un motivo a lo largo de los promotores.
# Entrada:
#   mid        → motif_id
#   df_pos     → data.frame con posiciones relativas (incluye todos los motivos)
# Salida:   archivo PNG guardado
# Por qué:  Si un motivo se concentra cerca de 0 (inicio del promotor, lejos
#           del TSS en promotores definidos upstream) o cerca de 1 (próximo
#           al gen), puede tener relevancia funcional diferente a si se
#           distribuye uniformemente.
# -----------------------------------------------------------------------------
grafico_distribucion_posicional <- function(mid, df_pos) {

  sub <- df_pos[df_pos$motif_id == mid, ]
  if (nrow(sub) == 0) {
    message("Sin ocurrencias para el motivo ", mid, ". Gráfico posicional omitido.")
    return(invisible(NULL))
  }

  p <- ggplot(sub, aes(x = pos_relativa)) +
    geom_histogram(aes(y = after_stat(density)), bins = n_bins_relativo,
                   fill = "#1D3557", color = "white", linewidth = 0.3) +
    geom_density(color = "#E63946", linewidth = 1) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_x_continuous(labels = scales::percent_format()) +
    labs(
      title = paste("Distribución posicional del motivo:", mid),
      subtitle = "Posición relativa dentro del promotor (0 = inicio, 1 = final)",
      x = "Posición relativa en el promotor",
      y = "Densidad de ocurrencias"
    ) +
    theme_bw(base_size = 13)

  guardar_grafico(p, paste0("distribucion_posicional_", mid, ".png"), ancho = 8, alto = 5)
  return(invisible(p))
}

# -----------------------------------------------------------------------------
# FUNCIÓN: grafico_mapa_posicional
# -----------------------------------------------------------------------------
# Qué hace: Gráfico de barras con el número de ocurrencias por bin posicional
#           relativo para todos los motivos significativos juntos o por motivo.
# Entrada:  df_bins → tabla de bins (salida de calcular_bins_posicionales)
# Salida:   archivo PNG guardado
# Por qué:  Permite detectar si alguna región del promotor concentra
#           consistentemente más motivos que otras zonas.
# -----------------------------------------------------------------------------
grafico_mapa_posicional <- function(df_bins) {

  if (nrow(df_bins) == 0) {
    message("Sin datos de bins posicionales. Gráfico omitido.")
    return(invisible(NULL))
  }

  p <- ggplot(df_bins, aes(x = bin, y = n_ocurrencias, fill = motif_id)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = "Distribución posicional de ocurrencias por motivo (bins relativos)",
      x = "Bin posicional (posición relativa en el promotor)",
      y = "Número de ocurrencias",
      fill = "Motivo"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))

  guardar_grafico(p, "mapa_posicional_relativo.png", ancho = 14, alto = 6)
  return(invisible(p))
}

# =============================================================================
# PIPELINE PRINCIPAL
# =============================================================================
# Esta función orquesta todo el análisis en orden. Cada paso está documentado
# con qué entra, qué hace y qué produce.
# -----------------------------------------------------------------------------
ejecutar_pipeline <- function() {

  message("\n========================================================")
  message("INICIO DEL PIPELINE DE ANÁLISIS DE MOTIVOS EN PROMOTORES")
  message("========================================================\n")

  # ------------------------------------------------------------------
  # PASO 1: CARGA DE DATOS
  # ------------------------------------------------------------------
  message(">>> PASO 1: Carga y validación de datos")

  # Leemos el archivo FASTA y lo convertimos en un data.frame limpio
  df_prom_raw <- leer_fasta(ruta_fasta)

  # Validamos y limpiamos el conjunto de promotores
  df_promotores <- validar_promotores(df_prom_raw)

  # Leemos el archivo de motivos
  df_motivos <- leer_motivos(ruta_motivos)

  # ------------------------------------------------------------------
  # PASO 2: BÚSQUEDA DE MOTIVOS EN PROMOTORES REALES
  # ------------------------------------------------------------------
  message("\n>>> PASO 2: Búsqueda de motivos en promotores reales")
  message("Opción IUPAC: ", usar_iupac, " | Reverso-complementario: ", buscar_rev_comp)

  df_ocurrencias <- buscar_todos_motivos(
    df_promotores   = df_promotores,
    df_motivos      = df_motivos,
    usar_iupac      = usar_iupac,
    buscar_rev_comp = buscar_rev_comp
  )

  message("Total de ocurrencias encontradas en promotores reales: ", nrow(df_ocurrencias))

  # Calculamos las métricas observadas por motivo
  df_obs <- calcular_metricas_observadas(df_ocurrencias, df_promotores)

  # ------------------------------------------------------------------
  # PASO 3: CÁLCULO DE COMPOSICIÓN DE BASES PARA EL MODELO NULO
  # ------------------------------------------------------------------
  message("\n>>> PASO 3: Cálculo de frecuencias de bases para simulación nula")
  message("Modelo nulo seleccionado: ", modelo_nulo)

  # Si el modelo es por permutación, las frecuencias no son necesarias
  # (se permutan directamente las secuencias). Para mononucleotídico y
  # dinucleotídico, calculamos las probabilidades del conjunto real.
  if (modelo_nulo != "permutacion") {
    probs_bases <- calcular_frecuencias_bases(df_promotores, modelo_nulo)
  } else {
    probs_bases <- NULL
    message("Modelo de permutación: no se precalculan frecuencias globales.")
  }

  # ------------------------------------------------------------------
  # PASO 4: SIMULACIÓN MONTE CARLO + ESTADÍSTICOS
  # ------------------------------------------------------------------
  message("\n>>> PASO 4: Simulación Monte Carlo (", n_sim, " iteraciones por motivo)")
  message("Este paso puede tardar varios minutos dependiendo del tamaño del dataset.")

  ids_motivos  <- unique(df_motivos$motif_id)
  lista_sim    <- list()     # Guardamos las distribuciones nulas
  lista_stats  <- list()     # Guardamos los estadísticos por motivo

  for (i in seq_along(ids_motivos)) {
    mid       <- ids_motivos[i]
    variantes <- df_motivos$sequence[df_motivos$motif_id == mid]

    message(sprintf("  Monte Carlo [%d/%d]: %s", i, length(ids_motivos), mid))

    # Valor observado para este motivo (puede ser 0 si no se encontró)
    obs_row   <- df_obs[df_obs$motif_id == mid, ]
    obs_prom  <- if (nrow(obs_row) > 0) obs_row$observed_promoters  else 0L
    obs_ocurr <- if (nrow(obs_row) > 0) obs_row$observed_occurrences else 0L

    # Ejecutamos las simulaciones Monte Carlo para este motivo
    sim_res <- ejecutar_monte_carlo_motivo(
      mid          = mid,
      variantes    = variantes,
      df_promotores = df_promotores,
      probs        = probs_bases,
      modelo       = modelo_nulo,
      n_sim        = n_sim,
      usar_iupac   = usar_iupac,
      buscar_rev_comp = buscar_rev_comp
    )

    lista_sim[[mid]] <- sim_res

    # Calculamos los estadísticos para este motivo
    stats_mid <- calcular_estadisticos_simulacion(sim_res, obs_prom, obs_ocurr, mid)
    lista_stats[[mid]] <- stats_mid

    # Generamos el gráfico de distribución nula para este motivo
    grafico_distribucion_nula(mid, sim_res, obs_prom)
  }

  # Combinamos todos los estadísticos en una tabla única
  df_stats <- do.call(rbind, lista_stats)

  # ------------------------------------------------------------------
  # PASO 5: CONSTRUCCIÓN DE LA TABLA RESUMEN FINAL
  # ------------------------------------------------------------------
  message("\n>>> PASO 5: Construcción de la tabla resumen final")

  # Contamos el número de variantes por motif_id
  n_variantes_por_motivo <- df_motivos %>%
    group_by(motif_id) %>%
    summarise(n_variantes = n(), .groups = "drop")

  # Unimos métricas observadas con estadísticos simulados
  # Primero aseguramos que todos los motivos aparecen (incluso con 0 ocurrencias)
  df_obs_completo <- data.frame(motif_id = ids_motivos) %>%
    left_join(df_obs, by = "motif_id") %>%
    mutate(
      observed_promoters   = ifelse(is.na(observed_promoters), 0L, observed_promoters),
      observed_occurrences = ifelse(is.na(observed_occurrences), 0L, observed_occurrences),
      freq_relativa        = ifelse(is.na(freq_relativa), 0, freq_relativa),
      densidad_ocurrencias = ifelse(is.na(densidad_ocurrencias), 0, densidad_ocurrencias)
    )

  df_final <- df_obs_completo %>%
    left_join(df_stats, by = "motif_id") %>%
    left_join(n_variantes_por_motivo, by = "motif_id")

  # Ajuste de p-valores por múltiples comparaciones (control del FDR)
  # Corrige el problema de que al hacer muchas pruebas simultáneas,
  # algunos resultados significativos pueden ser falsos positivos por azar.
  df_final$padj_promoters   <- p.adjust(df_final$pvalue_promoters, method = metodo_ajuste)
  df_final$padj_occurrences <- p.adjust(df_final$pvalue_occurrences, method = metodo_ajuste)

  # Ordenamos por significación en promotores (padj menor primero)
  df_final <- df_final %>%
    arrange(padj_promoters, desc(enrichment_promoters))

  # Seleccionamos y ordenamos las columnas de la tabla final
  df_final <- df_final %>%
    select(
      motif_id, n_variantes,
      observed_promoters, observed_occurrences, freq_relativa, densidad_ocurrencias,
      expected_promoters_mean, expected_promoters_q025, expected_promoters_q975,
      expected_occurrences_mean, expected_occurrences_q025, expected_occurrences_q975,
      enrichment_promoters, enrichment_occurrences,
      pvalue_promoters, pvalue_occurrences,
      padj_promoters, padj_occurrences
    )

  # ------------------------------------------------------------------
  # PASO 6: ANÁLISIS POSICIONAL
  # ------------------------------------------------------------------
  message("\n>>> PASO 6: Análisis posicional")

  if (nrow(df_ocurrencias) > 0) {
    # Calculamos posiciones relativas (0-1) para cada ocurrencia
    df_pos <- calcular_posiciones_relativas(df_ocurrencias, df_promotores)

    # Calculamos bins posicionales relativos
    df_bins_rel <- calcular_bins_posicionales(
      df_pos       = df_pos,
      tipo         = "relativo",
      n_bins_rel   = n_bins_relativo
    )

    # Calculamos bins posicionales absolutos (en pb)
    df_bins_abs <- calcular_bins_posicionales(
      df_pos      = df_pos,
      tipo        = "absoluto",
      tam_bin_abs = tamano_bin_abs
    )

    # Generamos gráfico de distribución posicional por motivo
    for (mid in ids_motivos) {
      grafico_distribucion_posicional(mid, df_pos)
    }

    # Gráfico resumen de bins posicionales relativos
    grafico_mapa_posicional(df_bins_rel)

  } else {
    df_pos      <- data.frame()
    df_bins_rel <- data.frame()
    df_bins_abs <- data.frame()
    message("Sin ocurrencias: análisis posicional omitido.")
  }

  # ------------------------------------------------------------------
  # PASO 7: VISUALIZACIONES GLOBALES
  # ------------------------------------------------------------------
  message("\n>>> PASO 7: Generando gráficos globales")

  grafico_enriquecimiento(df_final)
  grafico_ocurrencias_observadas(df_final)

  # ------------------------------------------------------------------
  # PASO 8: EXPORTACIÓN DE RESULTADOS
  # ------------------------------------------------------------------
  message("\n>>> PASO 8: Exportación de resultados a CSV")

  # Tabla de ocurrencias reales: una fila por ocurrencia
  write.csv(df_ocurrencias,
            file.path(ruta_salida, "tabla_ocurrencias_reales.csv"),
            row.names = FALSE)
  message("Exportado: tabla_ocurrencias_reales.csv")

  # Tabla resumen de motivos: una fila por motif_id
  write.csv(df_final,
            file.path(ruta_salida, "tabla_resumen_motivos.csv"),
            row.names = FALSE)
  message("Exportado: tabla_resumen_motivos.csv")

  # Tabla de posiciones relativas
  if (nrow(df_pos) > 0) {
    write.csv(df_pos,
              file.path(ruta_salida, "tabla_posiciones_relativas.csv"),
              row.names = FALSE)
    message("Exportado: tabla_posiciones_relativas.csv")
  }

  # Tabla de bins posicionales relativos
  if (nrow(df_bins_rel) > 0) {
    write.csv(df_bins_rel,
              file.path(ruta_salida, "tabla_bins_posicionales_relativos.csv"),
              row.names = FALSE)
    message("Exportado: tabla_bins_posicionales_relativos.csv")
  }

  # Tabla de bins posicionales absolutos
  if (nrow(df_bins_abs) > 0) {
    write.csv(df_bins_abs,
              file.path(ruta_salida, "tabla_bins_posicionales_absolutos.csv"),
              row.names = FALSE)
    message("Exportado: tabla_bins_posicionales_absolutos.csv")
  }

  # ------------------------------------------------------------------
  # PASO 9: RESUMEN TEXTUAL DEL ANÁLISIS
  # ------------------------------------------------------------------
  message("\n>>> PASO 9: Generando resumen textual")

  n_sig_prom  <- sum(df_final$padj_promoters < umbral_padj, na.rm = TRUE)
  n_sig_ocurr <- sum(df_final$padj_occurrences < umbral_padj, na.rm = TRUE)

  resumen_texto <- c(
    "=================================================",
    "RESUMEN DEL ANÁLISIS DE MOTIVOS EN PROMOTORES",
    "=================================================",
    paste("Fecha de análisis:      ", format(fecha_inicio, "%Y-%m-%d %H:%M:%S")),
    paste("Duración total:         ", round(as.numeric(Sys.time() - fecha_inicio, units = "mins"), 2), "minutos"),
    "",
    "--- DATOS DE ENTRADA ---",
    paste("Archivo FASTA:          ", ruta_fasta),
    paste("Archivo de motivos:     ", ruta_motivos),
    paste("Carpeta de salida:      ", ruta_salida),
    "",
    "--- PROMOTORES ANALIZADOS ---",
    paste("Número de promotores:   ", nrow(df_promotores)),
    paste("Longitud media (pb):    ", round(mean(df_promotores$length), 1)),
    paste("Longitud mínima (pb):   ", min(df_promotores$length)),
    paste("Longitud máxima (pb):   ", max(df_promotores$length)),
    "",
    "--- MOTIVOS ---",
    paste("Familias de motivos:    ", length(ids_motivos)),
    paste("Variantes totales:      ", nrow(df_motivos)),
    "",
    "--- PARÁMETROS DE BÚSQUEDA ---",
    paste("Uso de IUPAC:           ", usar_iupac),
    paste("Reverso-complementario: ", buscar_rev_comp),
    paste("Modelo nulo:            ", modelo_nulo),
    paste("Simulaciones (n_sim):   ", n_sim),
    paste("Semilla aleatoria:      ", semilla),
    paste("Ajuste p-valores:       ", metodo_ajuste),
    paste("Umbral padj:            ", umbral_padj),
    "",
    "--- RESULTADOS ---",
    paste("Motivos significativos (promotores, padj <", umbral_padj, "):", n_sig_prom),
    paste("Motivos significativos (ocurrencias, padj <", umbral_padj, "):", n_sig_ocurr),
    "",
    "--- MOTIVOS SIGNIFICATIVOS (ordenados por padj) ---"
  )

  sig_motivos <- df_final[df_final$padj_promoters < umbral_padj, ]
  if (nrow(sig_motivos) > 0) {
    for (k in seq_len(nrow(sig_motivos))) {
      r <- sig_motivos[k, ]
      resumen_texto <- c(resumen_texto,
        sprintf("  %s | obs_prom=%d | enrichment=%.2f | padj=%.4f",
                r$motif_id, r$observed_promoters,
                r$enrichment_promoters, r$padj_promoters))
    }
  } else {
    resumen_texto <- c(resumen_texto, "  (Ningún motivo superó el umbral de significación)")
  }

  resumen_texto <- c(resumen_texto, "=================================================")

  # Imprimimos el resumen en consola
  cat(paste(resumen_texto, collapse = "\n"), "\n")

  # Lo guardamos también como archivo de texto
  writeLines(resumen_texto,
             file.path(ruta_salida, "resumen_analisis.txt"))
  message("Exportado: resumen_analisis.txt")

  # ------------------------------------------------------------------
  # FIN DEL PIPELINE
  # ------------------------------------------------------------------
  message("\n=== ANÁLISIS COMPLETADO ===")
  message("Todos los resultados están en: ", normalizePath(ruta_salida))

  # Devolvemos una lista con todos los objetos principales por si el usuario
  # quiere inspeccionarlos directamente en la sesión de R
  return(invisible(list(
    promotores        = df_promotores,
    motivos           = df_motivos,
    ocurrencias       = df_ocurrencias,
    posiciones        = if (exists("df_pos")) df_pos else NULL,
    bins_relativos    = if (exists("df_bins_rel")) df_bins_rel else NULL,
    resumen_final     = df_final,
    simulaciones      = lista_sim
  )))
}

# =============================================================================
# EJECUCIÓN DEL PIPELINE
# =============================================================================
# Esta línea arranca todo el análisis. Es la única llamada que el usuario
# no necesita modificar.

resultados <- ejecutar_pipeline()

# El objeto 'resultados' contiene todos los data.frames principales para que
# el usuario pueda explorarlos en la consola de R si lo desea:
#   resultados$promotores     → datos de los promotores leídos
#   resultados$motivos        → definición de motivos cargada
#   resultados$ocurrencias    → tabla de ocurrencias reales
#   resultados$posiciones     → tabla de posiciones relativas
#   resultados$resumen_final  → tabla principal con estadísticos por motivo
#   resultados$simulaciones   → lista con distribuciones nulas por motivo
