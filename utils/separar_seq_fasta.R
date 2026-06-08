# =============================================================================
# separar_seq_fasta.R
# Splits a multi-sequence FASTA file into individual FASTA files,
# one per sequence, named after the sequence header.
# =============================================================================

library(Biostrings)

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
input_fasta <- "all_sequences.fasta"  # Input multi-sequence FASTA file
output_dir  <- "fasta_split"          # Output directory (created automatically)
# ──────────────────────────────────────────────────────────────────────────────

# Create output directory if it does not exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Read the multi-sequence FASTA file
archivo <- readAAStringSet(input_fasta)

message("Splitting ", length(archivo), " sequences into individual FASTA files...")

# Write one FASTA file per sequence, named after the sequence header
for (i in seq_along(archivo)) {
  nombre <- names(archivo)[i]
  filepath <- file.path(output_dir, paste0(nombre, ".fasta"))
  writeXStringSet(archivo[i], filepath = filepath)
}

message("Done! ", length(archivo), " files written to: ", output_dir)
