
# =============================================================================
# EDTA TE Annotation Report Generator
# Genera un informe de anotación de elementos transponibles a partir de un
# archivo GFF3 producido por EDTA/panEDTA.
#
# Uso:
#   Rscript EDTA_TE_report.R input.gff3 genome_size_bp output_prefix
#
# Argumentos:
#   input.gff3       : ruta al archivo GFF3 de EDTA
#   genome_size_bp   : tamaño total del genoma en pb (p.ej. 91963648)
#   output_prefix    : prefijo para los archivos de salida (default: "TE_report")
#
# Nota metodológica:
#   El informe original de EDTA se genera desde el archivo .out de RepeatMasker,
#   que resuelve solapamientos entre anotaciones (mayor score gana). El GFF3
#   mantiene TODAS las anotaciones, incluidas las solapantes. Por ello, los
#   conteos de este script pueden ser ligeramente superiores (~5%) a los del
#   informe original. Las bp de los elementos homológicos son exactas; las de
#   los estructurales usan el body del LTRRT (sin TSD).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(beepr)
})

# ===========================================================================
# 0. Argumentos / configuración
# ===========================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  gff_file     <- "input.gff3"   # <-- cambia aquí
  genome_size  <- 1000000000    # <-- cambia aquí (tamaño del genoma en bp)
  out_prefix   <- "TE_report"
} else {
  gff_file    <- args[1]
  genome_size <- as.numeric(args[2])
  out_prefix  <- ifelse(length(args) >= 3, args[3], "TE_report")
}

cat("==========================================================\n")
cat("  EDTA TE Annotation Report\n")
cat("==========================================================\n")
cat("GFF3:", gff_file, "\n")
cat("Genome size:", format(genome_size, big.mark = ","), "bp\n")
cat("Output prefix:", out_prefix, "\n\n")

# ===========================================================================
# 1. Leer el GFF3
# ===========================================================================
cat("Leyendo GFF3...\n")

gff_raw <- fread(
  gff_file,
  header    = FALSE,
  comment.char = "#", # Ignora automáticamente las líneas con almohadillas
  sep       = "\t",
  skip      = "#",
  col.names = c("seqid", "source", "type", "start", "end",
                "score", "strand", "phase", "attributes"),
  quote     = ""
)

# Eliminar líneas de comentario que puedan haber quedado
gff_raw <- gff_raw[!grepl("^#", seqid)]

cat("  Features totales leídos:", nrow(gff_raw), "\n")

# ===========================================================================
# 2. Parsear atributos clave
# ===========================================================================
parse_attr <- function(attrs, key) {
  pattern <- paste0("(?:^|;)", key, "=([^;]+)")
  m <- regmatches(attrs, regexpr(pattern, attrs, perl = TRUE))
  ifelse(nchar(m) == 0, NA_character_,
         sub(paste0(".*", key, "=([^;]+).*"), "\\1", m))
}

gff <- gff_raw[, .(
  seqid, type, start, end,
  score,
  name           = parse_attr(attributes, "Name"),
  classification = parse_attr(attributes, "classification"),
  method         = parse_attr(attributes, "method"),
  has_parent     = grepl("Parent=", attributes),
  bp             = end - start + 1
)]

cat("  Atributos parseados OK\n\n")

# ===========================================================================
# 3. Selección de features para el conteo
# ===========================================================================
# Para la tabla Repeat Classes y Repeat Stats usamos:
# - Excluimos: long_terminal_repeat (hijos LTR estructurales)
#              target_site_duplication (hijos TSD estructurales)
#              repeat_region (padre estructural LTR, se cuenta via hijo LTRRT)
# - Incluimos: todos los demás (con o sin Parent=, siempre que no sean LTR/TSD padres)
#
# Para los LTRRT con Parent (hijos de repeat_region estructurales):
#   Tipo: Gypsy_LTR_retrotransposon, Copia_LTR_retrotransposon, LTR_retrotransposon
#   Estos representan el body del TE (sin TSD). Se incluyen.

exclude_types <- c("long_terminal_repeat", "target_site_duplication", "repeat_region")

te_features <- gff[!type %in% exclude_types]

cat("Features para análisis (excluidos LTR/TSD/repeat_region padres):",
    nrow(te_features), "\n\n")

# ===========================================================================
# 4. Mapeo de clasificaciones a clases del informe
# ===========================================================================
# La lógica usa el TIPO de feature (columna 3 del GFF) para asignar la clase.
# Esto es más robusto que usar el campo classification porque un mismo 
# Mutator_TIR_transposon puede tener classification=DNA/DTM o MITE/DTM.

map_type_to_class <- function(type_vec) {
  dplyr::case_when(
    type_vec == "L1_LINE_retrotransposon"         ~ "LINE/L1",
    type_vec == "RTE_LINE_retrotransposon"         ~ "LINE/RTE",
    type_vec %in% c("Gypsy_LTR_retrotransposon")  ~ "LTR/Gypsy",
    type_vec %in% c("Copia_LTR_retrotransposon")  ~ "LTR/Copia",
    type_vec == "LTR_retrotransposon"              ~ "LTR/unknown",
    type_vec == "tRNA_SINE_retrotransposon"        ~ "SINE/tRNA",
    type_vec == "CACTA_TIR_transposon"             ~ "TIR/CACTA",
    type_vec == "Mutator_TIR_transposon"           ~ "TIR/Mutator",
    type_vec == "PIF_Harbinger_TIR_transposon"     ~ "TIR/PIF_Harbinger",
    type_vec == "Tc1_Mariner_TIR_transposon"       ~ "TIR/Tc1_Mariner",
    type_vec == "hAT_TIR_transposon"               ~ "TIR/hAT",
    type_vec == "low_complexity"                   ~ "low_complexity",
    type_vec == "DIRS_YR_retrotransposon"          ~ "nonLTR/DIRS_YR",
    type_vec == "pararetrovirus"                   ~ "nonLTR/pararetrovirus",
    # GFF usa "pararetrovirus" (con una r), report usa "pararetrovirus"
    type_vec == "pararetrovirus"                   ~ "nonLTR/pararetrovirus",
    type_vec == "helitron"                         ~ "nonTIR/helitron",
    type_vec == "rRNA_gene"                        ~ "rDNA/45S",
    type_vec == "repeat_fragment"                  ~ "repeat_fragment",
    TRUE                                           ~ paste0("other/", type_vec)
  )
}

# Jerarquía de clases para agrupar en la tabla principal
map_class_to_superclass <- function(class_vec) {
  dplyr::case_when(
    class_vec %in% c("LINE/L1", "LINE/RTE")                        ~ "LINE",
    class_vec %in% c("LTR/Gypsy", "LTR/Copia", "LTR/unknown")     ~ "LTR",
    class_vec == "SINE/tRNA"                                        ~ "SINE",
    class_vec %in% c("TIR/CACTA", "TIR/Mutator", "TIR/PIF_Harbinger",
                     "TIR/Tc1_Mariner", "TIR/hAT")                 ~ "TIR",
    class_vec == "low_complexity"                                   ~ "low_complexity",
    class_vec %in% c("nonLTR/DIRS_YR", "nonLTR/pararetrovirus")    ~ "nonLTR",
    class_vec == "nonTIR/helitron"                                  ~ "nonTIR",
    class_vec == "rDNA/45S"                                         ~ "rDNA",
    class_vec == "repeat_fragment"                                  ~ "repeat_fragment",
    TRUE                                                            ~ "other"
  )
}

te_features[, te_class := map_type_to_class(type)]

# ===========================================================================
# 5. Tabla Repeat Classes
# ===========================================================================
cat("==========================================================\n")
cat("Repeat Classes\n")
cat("==============\n")
cat(sprintf("%-20s%10s%15s%10s\n", "Total Sequences:", length(unique(te_features$seqid)), "", ""))
cat(sprintf("%-20s%10s%15s%10s\n", "Total Length:",
            format(genome_size, big.mark = ","), "bp", ""))
cat(sprintf("\n%-30s%10s%15s%10s\n", "Class", "Count", "bpMasked", "%masked"))
cat(sprintf("%-30s%10s%15s%10s\n", "=====", "=====", "========", "======="))

# Calcular estadísticas por clase
class_stats <- te_features[, .(
  Count    = .N,
  bpMasked = sum(bp)
), by = te_class][order(te_class)]

class_stats[, pct := round(bpMasked / genome_size * 100, 2)]
class_stats[, superclass := map_class_to_superclass(te_class)]

# Orden de superclases como en el informe
superclass_order <- c("LINE", "LTR", "SINE", "TIR", "low_complexity",
                      "nonLTR", "nonTIR", "rDNA", "repeat_fragment", "other")

# Imprimir con estructura jerárquica
print_class_table <- function(df, genome_size) {
  lines <- list()
  
  for (sc in superclass_order) {
    sub <- df[superclass == sc]
    if (nrow(sub) == 0) next
    
    if (sc %in% c("LINE", "LTR", "SINE", "TIR", "nonLTR", "nonTIR", "rDNA")) {
      lines[[length(lines) + 1]] <- sprintf("%-30s%10s%15s%10s",
                                            sc, "--", "--", "--")
      for (i in seq_len(nrow(sub))) {
        cl <- sub[i, te_class]
        family <- sub("^[^/]+/", "", cl)
        lines[[length(lines) + 1]] <- sprintf("    %-26s%10s%15s%9s%%",
                                              family,
                                              format(sub[i, Count], big.mark = ","),
                                              format(sub[i, bpMasked], big.mark = ","),
                                              sub[i, pct])
      }
    } else {
      for (i in seq_len(nrow(sub))) {
        cl <- sub[i, te_class]
        lines[[length(lines) + 1]] <- sprintf("%-30s%10s%15s%9s%%",
                                              cl,
                                              format(sub[i, Count], big.mark = ","),
                                              format(sub[i, bpMasked], big.mark = ","),
                                              sub[i, pct])
      }
    }
  }
  
  # Total interspersed (excluye repeat_fragment y low_complexity)
  interspersed <- df[!te_class %in% c("repeat_fragment", "low_complexity")]
  total_int_count <- sum(interspersed$Count)
  total_int_bp    <- sum(interspersed$bpMasked)
  total_int_pct   <- round(total_int_bp / genome_size * 100, 2)
  
  lines[[length(lines) + 1]] <- sprintf("%45s", "---------------------------------")
  lines[[length(lines) + 1]] <- sprintf("    %-26s%10s%15s%9s%%",
                                        "total interspersed",
                                        format(total_int_count, big.mark = ","),
                                        format(total_int_bp, big.mark = ","),
                                        total_int_pct)
  lines[[length(lines) + 1]] <- ""
  lines[[length(lines) + 1]] <- strrep("-", 57)
  
  # Total global
  total_count <- sum(df$Count)
  total_bp    <- sum(df$bpMasked)
  total_pct   <- round(total_bp / genome_size * 100, 2)
  
  lines[[length(lines) + 1]] <- sprintf("%-30s%10s%15s%9s%%",
                                        "Total",
                                        format(total_count, big.mark = ","),
                                        format(total_bp, big.mark = ","),
                                        total_pct)
  
  invisible(lapply(lines, cat, "\n"))
  return(list(
    class_stats = df,
    total_count = total_count,
    total_bp    = total_bp,
    total_pct   = total_pct
  ))
}

result <- print_class_table(class_stats, genome_size)

# ===========================================================================
# 6. Tabla Repeat Stats (por Name/familia)
# ===========================================================================
cat("\n\nRepeat Stats\n")
cat("============\n")
cat(sprintf("%-20s%10s\n", "Total Sequences:", length(unique(te_features$seqid))))
cat(sprintf("%-20s%10s%5s\n", "Total Length:",
            format(genome_size, big.mark = ","), "bp"))
cat(sprintf("\n%-40s%10s%15s%10s\n",
            "ID", "Count", "bpMasked", "%masked"))
cat(sprintf("%-40s%10s%15s%10s\n",
            "================", "=====", "========", "======="))

# Para el Repeat Stats usamos el campo 'name' del GFF
# Separamos:
#  a) TEs homológicos y TIR/helitron/LINE/DIRS: directamente por Name
#  b) TEs estructurales LTR (LTRRT con Parent=): su Name se reporta como Name del repeat_region

# a) Features NO estructurales LTR (todos excepto LTRRT hijos de repeat_region)
#    Los LTRRT hijos tienen: has_parent=TRUE y type en c(Gypsy/Copia/LTR_retrotransposon)
structural_ltr_types <- c("Gypsy_LTR_retrotransposon",
                          "Copia_LTR_retrotransposon",
                          "LTR_retrotransposon")

# Features no estructurales (homología + TIR + LINE + helitron + etc)
non_struct <- te_features[!(type %in% structural_ltr_types & has_parent == TRUE)]

# Features estructurales LTR (LTRRT hijos): usamos sus bp
struct_ltr <- te_features[type %in% structural_ltr_types & has_parent == TRUE]

# Combinar: para los estructurales el Name viene del campo name (= nombre de la familia)
# Para los no-estructurales igual
all_for_stats <- rbind(non_struct, struct_ltr)

# Calcular estadísticas por Name
name_stats <- all_for_stats[!is.na(name), .(
  Count    = .N,
  bpMasked = sum(bp)
), by = name][order(name)]

name_stats[, pct := round(bpMasked / genome_size * 100, 2)]

# Imprimir tabla Repeat Stats
for (i in seq_len(nrow(name_stats))) {
  cat(sprintf("%40s%10s%15s%9s%%\n",
              name_stats[i, name],
              format(name_stats[i, Count], big.mark = ","),
              format(name_stats[i, bpMasked], big.mark = ","),
              name_stats[i, pct]))
}

cat(strrep("-", 76), "\n")
cat(sprintf("%40s%10s%15s%9s%%\n",
            "",
            format(sum(name_stats$Count), big.mark = ","),
            format(sum(name_stats$bpMasked), big.mark = ","),
            round(sum(name_stats$bpMasked) / genome_size * 100, 2)))

# ===========================================================================
# 7. Tabla "By Sequence" (por cromosoma/secuencia)
# ===========================================================================
cat("\n\nBy Sequence\n")
cat("===========\n")
cat(sprintf("%-30s%10s%15s\n", "Seq", "Count", "bpMasked"))
cat(sprintf("%-30s%10s%15s\n", "=====", "=====", "========"))

seq_stats <- all_for_stats[!is.na(name), .(
  Count    = .N,
  bpMasked = sum(bp)
), by = seqid][order(seqid)]

for (i in seq_len(nrow(seq_stats))) {
  cat(sprintf("%-30s%10s%15s\n",
              seq_stats[i, seqid],
              format(seq_stats[i, Count], big.mark = ","),
              format(seq_stats[i, bpMasked], big.mark = ",")))
}

# ===========================================================================
# 8. Guardar resultados en CSV
# ===========================================================================
class_out <- class_stats[, .(
  te_class, superclass, Count, bpMasked,
  pct_masked = pct
)]
fwrite(class_out, paste0(out_prefix, "_repeat_classes.csv"))

name_out <- name_stats
fwrite(name_out, paste0(out_prefix, "_repeat_stats.csv"))

seq_out <- seq_stats
fwrite(seq_out, paste0(out_prefix, "_by_sequence.csv"))

cat("\n\nArchivos generados:\n")
cat("  -", paste0(out_prefix, "_repeat_classes.csv"), "\n")
cat("  -", paste0(out_prefix, "_repeat_stats.csv"), "\n")
cat("  -", paste0(out_prefix, "_by_sequence.csv"), "\n")

# ===========================================================================
# 9. Resumen rápido en consola
# ===========================================================================
cat("\n==========================================================\n")
cat("RESUMEN\n")
cat("==========================================================\n")
cat(sprintf("  Secuencias analizadas : %d\n", length(unique(te_features$seqid))))
cat(sprintf("  Tamaño del genoma     : %s bp\n", format(genome_size, big.mark = ",")))
cat(sprintf("  Features TE totales   : %s\n",
            format(result$total_count, big.mark = ",")))
cat(sprintf("  bp enmascarados       : %s (%.2f%%)\n",
            format(result$total_bp, big.mark = ","), result$total_pct))
cat(sprintf("  Familias únicas       : %d\n", nrow(name_stats)))
cat("==========================================================\n")
beep(2)
