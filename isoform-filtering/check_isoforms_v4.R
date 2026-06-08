# =============================================================================
# check_isoforms_v4.R
# Limpieza de isoformas proteicas: selección de una proteína representativa
# por locus a partir de archivos GFF3 y FASTA.
#
# Criterio de selección:
#   1º Mayor longitud de proteína (aa) — proxy directo de CDS
#   2º Mayor longitud total de CDS    — desempate sumando todos los CDS
#   3º Menor número de isoforma       — reproducibilidad ante empates totales
#
# Compatible con GFFs NCBI con o sin columna locus_tag:
#   - Q. robur: tiene columna locus_tag rellena en genes
#   - Q. suber: no tiene locus_tag, usa columna gene (LOC112016752)
# =============================================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(rtracklayer)
  library(dplyr)
  library(cli)
  library(beepr)
})

# =============================================================================
# FUNCIÓN PRINCIPAL
# =============================================================================

check_isoforms <- function(gff_file, faa_file, output_dir) {
  
  # ---------------------------------------------------------------------------
  # 0. Validación de inputs
  # ---------------------------------------------------------------------------
  stopifnot(
    "gff_file no existe" = file.exists(gff_file),
    "faa_file no existe" = file.exists(faa_file)
  )
  
  cli_h1("Limpieza de isoformas")
  cli_alert_info("GFF  : {.path {gff_file}}")
  cli_alert_info("FASTA: {.path {faa_file}}")
  
  # ---------------------------------------------------------------------------
  # 1. Importar GFF
  # ---------------------------------------------------------------------------
  cli_progress_step("Importando GFF...")
  gff_df <- as.data.frame(import(gff_file))
  cli_alert_success("GFF importado: {nrow(gff_df)} features")
  
  # ---------------------------------------------------------------------------
  # 2a. Detectar qué columna usar como identificador de gen
  #
  #   Caso A — Q. robur (NCBI cromosómico): tiene columna locus_tag rellena
  #   Caso B — Q. suber (NCBI scaffold)  : sin locus_tag, usa columna gene
  #
  #   Jerarquía: locus_tag > gene > ID limpio (fallback)
  # ---------------------------------------------------------------------------
  genes_df <- gff_df[gff_df$type == "gene", ]
  
  has_locus_tag <- "locus_tag" %in% colnames(genes_df) &&
    any(!is.na(genes_df$locus_tag))
  
  has_gene_col  <- "gene" %in% colnames(genes_df) &&
    any(!is.na(genes_df$gene))
  
  id_col_used <- if (has_locus_tag) "locus_tag" else if (has_gene_col) "gene" else "ID (fallback)"
  cli_alert_info("Columna de ID de gen detectada: {.field {id_col_used}}")
  
  gene_tbl <- as_tibble(genes_df) |>
    transmute(
      gene_id   = as.character(ID),
      locus_tag = if (has_locus_tag) {
        if (has_gene_col) dplyr::coalesce(locus_tag, gene) else locus_tag
      } else if (has_gene_col) {
        gene
      } else {
        sub("^gene-", "", as.character(ID))
      }
    ) |>
    filter(!is.na(locus_tag))
  
  cli_alert_info("Genes con identificador de locus: {nrow(gene_tbl)}")
  
  # ---------------------------------------------------------------------------
  # 2b. Extraer tabla mRNA y propagar locus_tag desde el nivel gene
  # ---------------------------------------------------------------------------
  mrna_tbl <- as_tibble(gff_df) |>
    filter(type == "mRNA", !is.na(transcript_id)) |>
    transmute(
      mrna_id       = as.character(ID),
      parent_gene   = vapply(Parent, `[[`, character(1), 1),
      transcript_id,
      mrna_width    = width
    ) |>
    left_join(gene_tbl, by = c("parent_gene" = "gene_id")) |>
    filter(!is.na(locus_tag))
  
  n_mrna <- nrow(mrna_tbl)
  n_loci <- n_distinct(mrna_tbl$locus_tag)
  cli_alert_info("{n_mrna} mRNAs codificantes en {n_loci} loci")
  
  # ---------------------------------------------------------------------------
  # 3. Extraer tabla CDS: longitud total y protein_id por transcrito
  # ---------------------------------------------------------------------------
  cds_tbl <- as_tibble(gff_df) |>
    filter(type == "CDS", !is.na(protein_id)) |>
    transmute(
      parent_mrna = vapply(Parent, `[[`, character(1), 1),
      protein_id,
      width
    ) |>
    group_by(parent_mrna, protein_id) |>
    summarise(total_cds_length = sum(width), .groups = "drop")
  
  # ---------------------------------------------------------------------------
  # 4. Leer FASTA y calcular longitud de proteína en aa
  # ---------------------------------------------------------------------------
  cli_progress_step("Leyendo FASTA de proteínas...")
  fasta     <- readAAStringSet(faa_file)
  fasta_ids <- sub("^(\\S+).*", "\\1", names(fasta))
  
  protein_lengths <- tibble(
    protein_id     = fasta_ids,
    protein_length = width(fasta)
  )
  
  n_fasta <- length(fasta)
  cli_alert_info("{n_fasta} secuencias en el FASTA")
  
  # ---------------------------------------------------------------------------
  # 5. Unir toda la información
  # ---------------------------------------------------------------------------
  full_tbl <- mrna_tbl |>
    left_join(cds_tbl,         by = c("mrna_id" = "parent_mrna")) |>
    left_join(protein_lengths, by = "protein_id")
  
  n_missing <- sum(is.na(full_tbl$protein_id))
  if (n_missing > 0) {
    cli_alert_warning(
      "{n_missing} mRNA(s) sin protein_id tras el join. Revisa coherencia GFF-FASTA."
    )
  }
  
  # ---------------------------------------------------------------------------
  # 6. Seleccionar UNA proteína por locus
  # ---------------------------------------------------------------------------
  best_tbl <- full_tbl |>
    filter(!is.na(protein_id)) |>
    group_by(locus_tag) |>
    slice_max(
      order_by  = tibble(protein_length, total_cds_length, -mrna_width),
      n         = 1,
      with_ties = FALSE
    ) |>
    ungroup()
  
  n_selected <- nrow(best_tbl)
  cli_alert_success("{n_selected} loci con proteína representativa seleccionada")
  
  loci_lost <- setdiff(mrna_tbl$locus_tag, best_tbl$locus_tag)
  if (length(loci_lost) > 0) {
    cli_alert_warning(
      "{length(loci_lost)} loci sin proteína válida: {.val {head(loci_lost, 5)}}"
    )
  }
  
  # ---------------------------------------------------------------------------
  # 7. Filtrar FASTA
  # ---------------------------------------------------------------------------
  keep_ids       <- best_tbl$protein_id
  fasta_filtered <- fasta[fasta_ids %in% keep_ids]
  fasta_filtered <- fasta_filtered[match(keep_ids, fasta_ids[fasta_ids %in% keep_ids])]
  n_out          <- length(fasta_filtered)
  
  if (n_out == 0) {
    cli_abort(c(
      "x" = "El FASTA filtrado está vacío.",
      "i" = "Ejemplo GFF  : {.val {head(keep_ids, 3)}}",
      "i" = "Ejemplo FASTA: {.val {head(fasta_ids, 3)}}"
    ))
  }
  
  cli_alert_success("{n_out} proteínas seleccionadas")
  
  # ---------------------------------------------------------------------------
  # 8. Escribir FASTA de salida
  # ---------------------------------------------------------------------------
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  output_faa <- file.path(
    normalizePath(output_dir, mustWork = FALSE),
    sub("\\.faa$", "_cleaned_isoforms.faa", basename(faa_file))
  )
  
  writeXStringSet(fasta_filtered, filepath = output_faa, format = "fasta")
  cli_alert_success("FASTA guardado en: {.path {output_faa}}")
  
  # ---------------------------------------------------------------------------
  # 9. Estadísticas
  # ---------------------------------------------------------------------------
  stats <- tibble(
    total_mrna_features = n_mrna,
    total_loci          = n_loci,
    loci_with_protein   = n_selected,
    loci_lost           = length(loci_lost),
    fasta_input_seqs    = n_fasta,
    fasta_output_seqs   = n_out,
    isoforms_removed    = n_mrna - n_selected
  )
  
  cli_h2("Resumen")
  print(as.data.frame(t(stats)))
  
  invisible(list(
    table      = best_tbl,
    fasta      = fasta_filtered,
    output_faa = output_faa,
    stats      = stats
  ))
}


# =============================================================================
# EJECUCIÓN — edita las rutas y cambia FALSE → TRUE para ejecutar
# =============================================================================

# --- Especie 1 ---------------------------------------------------------------
if (FALSE) {
  result <- check_isoforms(
    gff_file   = "species1.gff",           # GFF/GFF3 annotation file
    faa_file   = "species1_protein.faa",   # Protein FASTA file
    output_dir = "output/species1"          # Output directory (created automatically)
  )
  result$stats
  beep(2)
}

# --- Especie 2 ---------------------------------------------------------------
if (FALSE) {
  result <- check_isoforms(
    gff_file   = "species2.gff",
    faa_file   = "species2_protein.faa",
    output_dir = "output/species2"
  )
  result$stats
  beep(2)
}