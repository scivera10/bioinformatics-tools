# =============================================================================
# extract_selected_proteins.R
# Extracts a subset of protein sequences from a FASTA database based on
# a list of accession IDs provided in an Excel file.
# =============================================================================

library(readxl)
library(Biostrings)

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
excel_file   <- "selected_proteins.xlsx"  # Excel file with protein IDs
id_column    <- "Accession"               # Column name containing the IDs
fasta_db     <- "proteome.faa"            # Input protein FASTA database
output_fasta <- "selected_proteins.faa"   # Output FASTA file
# ──────────────────────────────────────────────────────────────────────────────

# Read the Excel file
prot_id <- read_excel(excel_file)
ids <- unique(prot_id[[id_column]])
message("IDs to extract: ", length(ids))

# Load the protein database
db_raw <- readAAStringSet(fasta_db)
message("Sequences in database: ", length(db_raw))

# Extract matching sequences
db_ids   <- sub(" .*", "", names(db_raw))
prot_seq <- db_raw[db_ids %in% ids]
message("Sequences extracted: ", length(prot_seq))

# Write output FASTA
writeXStringSet(prot_seq, filepath = output_fasta)
message("Done! Written to: ", output_fasta)
