# =============================================================================
# extract_unique_blast_sequences.R
# Extract unique BLAST hit sequences from a multi-species reiterative BLAST
# results directory and merge them into a single annotated FASTA file.
#
# Expected directory structure:
#   results_dir/
#   ├── Species_A/        ← one folder per species
#   │   ├── blast_result1.xlsx
#   │   └── blast_result2.xlsx
#   ├── Species_B/
#   │   └── blast_result1.xlsx
#   └── ...
#
# Each .xlsx file must contain a 'Subject' column with protein IDs.
# =============================================================================

library(readxl)
library(writexl)
library(dplyr)
library(Biostrings)

# ── USER CONFIGURATION ────────────────────────────────────────────────────────

# Directory containing one subfolder per species with BLAST .xlsx results
results_dir <- "results/reiterative_blastp_results"

# Named list: species label → path to protein FASTA database
# Add or remove entries as needed
species_db <- list(
  "Species_A" = "data/raw/species_a_protein.faa",
  "Species_B" = "data/raw/species_b_protein.faa",
  "Species_C" = "data/raw/species_c_protein.faa"
)

# Optional: path to a FASTA with your own query sequences (set to NULL to skip)
query_fasta <- NULL   # e.g. "data/raw/my_query_proteins.fasta"
query_label <- NULL   # e.g. "(My species hap1)"

# Output file name
output_fasta <- "todas_las_secuencias.faa"

# ──────────────────────────────────────────────────────────────────────────────

# =============================================================================
# 1. Read BLAST results and deduplicate by Subject per species folder
# =============================================================================
setwd(results_dir)

direct_list <- trimws(list.files())

# Read all .xlsx files per species folder and keep unique Subject IDs
for (dir in direct_list) {

  ruta_dir <- file.path(".", dir)
  archivos <- list.files(path = ruta_dir, pattern = "\\.xlsx$", full.names = TRUE)

  # Skip Excel temp files and previously generated output
  archivos <- archivos[!grepl("~\\$", basename(archivos))]
  archivos <- archivos[basename(archivos) != "datos_unicos.xlsx"]

  if (length(archivos) == 0) next

  datos <- lapply(archivos, function(archivo) read_excel(archivo, col_types = "text"))
  datos_unidos <- bind_rows(datos)
  datos_unicos <- datos_unidos %>% distinct(Subject, .keep_all = TRUE)

  write_xlsx(datos_unicos, file.path(ruta_dir, "datos_unicos.xlsx"))
  message("Processed: ", dir, " — ", nrow(datos_unicos), " unique hits saved.")
}

# =============================================================================
# 2. Load deduplicated results into a named list
# =============================================================================
lista_datos <- lapply(direct_list, function(dir) {
  ruta <- file.path(".", dir, "datos_unicos.xlsx")
  if (file.exists(ruta)) read_excel(ruta) else NULL
})
names(lista_datos) <- direct_list
lista_datos <- Filter(Negate(is.null), lista_datos)

# =============================================================================
# 3. Load species FASTA databases and extract matching sequences
# =============================================================================
setwd("../../..")  # Return to project root

all_seqs_list <- list()

# Optional: add query sequences first
if (!is.null(query_fasta) && file.exists(query_fasta)) {
  seqs_query <- readAAStringSet(query_fasta)
  if (!is.null(query_label)) {
    names(seqs_query) <- paste0(names(seqs_query), " ", query_label)
  }
  all_seqs_list[["query"]] <- seqs_query
  message("Query sequences loaded: ", length(seqs_query))
}

# Extract sequences for each species
for (species_name in names(species_db)) {

  fasta_path <- species_db[[species_name]]

  if (!file.exists(fasta_path)) {
    warning("FASTA not found for ", species_name, ": ", fasta_path)
    next
  }

  if (!species_name %in% names(lista_datos)) {
    message("No BLAST results folder found for: ", species_name, " — skipping.")
    next
  }

  db <- readAAStringSet(fasta_path)
  ids_db <- sub(" .*", "", names(db))
  id_list <- unique(lista_datos[[species_name]]$Subject)

  seqs <- db[ids_db %in% id_list]
  names(seqs) <- paste0(names(seqs), " (", species_name, ")")

  all_seqs_list[[species_name]] <- seqs
  message(species_name, ": ", length(seqs), " sequences extracted.")
}

# =============================================================================
# 4. Merge all sequences and write output FASTA
# =============================================================================
todas_seqs <- do.call(c, all_seqs_list)

writeXStringSet(todas_seqs, output_fasta)
message("\nDone! ", length(todas_seqs), " sequences written to: ", output_fasta)
