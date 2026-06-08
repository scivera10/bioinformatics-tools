


# Load necessary libraries
library(ape)        # For handling sequences
library(Biostrings) # For handling FASTA files
library(tidyverse)  # For manage data
library(openxlsx)   # For Excel export with multiple sheets

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
# Set your working directory to the folder containing your FASTA files before
# running this script (Session > Set Working Directory in RStudio, or use setwd())

query_fasta_file   <- "your_query.fasta"      # Query protein FASTA file
subject_fasta_file <- "your_subject.faa"      # Subject proteome FASTA file
blastp_executable  <- "blastp"                # Path to blastp binary (or just "blastp" if in PATH)
output_excel       <- "Blast_prot_local.xlsx" # Output Excel file name
# ──────────────────────────────────────────────────────────────────────────────

# 1. BLASTp
# ------------------------------------------------------------------------------

# Run the blastp command in R
blast_result <- system2(
  command = blastp_executable,
  args = c(
    "-query", query_fasta_file,
    "-subject", subject_fasta_file,
    "-outfmt", "6"
  ),
  stdout = TRUE
)
print(blast_result)
#writeLines(blast_result, "blast_result.txt")


# Convert blast result to a data frame
BLAST_AB <- read.table(textConnection(blast_result), header = FALSE, col.names = c(
  "Query", "Subject", "Identity", "AlignmentLength", "Mismatches", "GapOpenings",
  "QueryStart", "QueryEnd", "SubjectStart", "SubjectEnd", "EValue", "BitScore"
), sep = "\t")

print(BLAST_AB)


# Calculate the coverage for each alignment
# Read query sequences from the FASTA file
query_sequences <- readAAStringSet(query_fasta_file)

# Keep only the first word of the sequence name
## Transforms "chr5_hap1 - CSH1..." into "chr5_hap1"
names(query_sequences) <- sapply(strsplit(names(query_sequences), " "), `[`, 1)

# Create a named vector of query sequence lengths
query_lengths <- width(query_sequences)
names(query_lengths) <- names(query_sequences)

# Add a new column for the query length
BLAST_AB$QueryLength <- query_lengths[BLAST_AB$Query]

# Adjust the coverage calculation to account for gaps
BLAST_AB$EffectiveAlignmentLength <- (BLAST_AB$QueryEnd - BLAST_AB$QueryStart + 1)

# Calculate the coverage considering the gaps
BLAST_AB$Coverage <- (BLAST_AB$EffectiveAlignmentLength / BLAST_AB$QueryLength) * 100

head(BLAST_AB)

# 2. Filter
# ------------------------------------------------------------------------------
# Apply the three filters
BLAST_filtered <- BLAST_AB %>%
  filter(
    Identity >= 30,           # 1. Identidad >= 30%
    Coverage >= 21,           # 2. Cobertura >= 19%
    EValue < 1e-10            # 3. E-value < a 1e-10
  )

# Save the result in a CSV
# write.csv(BLAST_filtered, "blast_filered_final.csv", row.names = FALSE)

# 3. Save results
# ------------------------------------------------------------------------------
# Load the database
subject_db <- readAAStringSet(subject_fasta_file)

# Clean the df headers
names(subject_db) <- word(names(subject_db), 1)

# Create a list of unique IDs from the first BLAST
unique_subjects <- unique(BLAST_filtered$Subject)

# Create the Excel Workbook
wb <- createWorkbook()

# Add the first sheet to the Excel file
addWorksheet(wb, "First_Search")

# Paste the results into the first sheet
writeData(wb, "First_Search", BLAST_filtered)


# 4.Start Reiterative BLASTp
# ------------------------------------------------------------------------------
# Create a df with all filtered results
all_results_filtered <- BLAST_filtered[,c("Query","Subject","Identity","EValue","Coverage")]

# Create a list of the queries used
queries_used <- c()

# Start Loop
for (target_id in unique_subjects) {
  # Store the query used
  queries_used <- c(queries_used, target_id)
  # Extract the specific sequence for the new query
  new_query_seq <- subject_db[target_id]
  writeXStringSet(new_query_seq, "temp_reciprocal.fasta")
  
  # Run BLAST for this specific subject
  reciprocal_out <- system2(
    command = blastp_executable,
    args = c("-query", "temp_reciprocal.fasta", "-subject", subject_fasta_file, "-outfmt", "6"),
    stdout = TRUE
  )
  
  # Safety check: Only continue if BLAST returns a result
  if (length(reciprocal_out) > 0) {
    
    # Convert raw text to a df
    res_table <- read.table(textConnection(reciprocal_out), sep = "\t")
    
    # Set the column names
    colnames(res_table) <- colnames(BLAST_AB)[1:12]
    
    # Length of the query sequence used in this reciprocal BLAST search
    query_length <- width(new_query_seq)
    
    # Add the query length
    res_table$QueryLength <- query_length
    
    # Compute the effective alignment length
    res_table$EffectiveAlignmentLength <- res_table$QueryEnd - res_table$QueryStart + 1
    
    # Compute coverage
    res_table$Coverage <- (res_table$EffectiveAlignmentLength / res_table$QueryLength) * 100
    
    # Apply the filter
    res_table_filtered <- res_table %>%
      filter(
        Identity >= 30,
        Coverage >= 21, # Excepto en la primera todas a 26
        EValue < 1e-10
      )
    
    # Save only if results remain after filtering
    if (nrow(res_table_filtered) > 0) {
      
      # Trim the protein name to 31 characters (Excel max)
      sheet_name <- substr(target_id, 1, 31)
      
      # Create a new sheet named after that protein
      addWorksheet(wb, sheet_name)
      
      # Paste the filtered results into that sheet
      writeData(wb, sheet_name, res_table_filtered)
      
      # Add a column indicating the source search
      res_table_filtered$ReciprocalQuery <- target_id
      
      # Add the data to the combined df
      all_results_filtered <- rbind(all_results_filtered, res_table_filtered[,c("Query", "Subject","Identity","EValue","Coverage")])

      message(paste("Processed reciprocal BLAST for:", target_id))
    }
  }
}

# Keep only unique entries
all_results_filtered <- all_results_filtered[!duplicated(all_results_filtered$Subject), ]
queries_used <- unique(queries_used)

# Create a new sheet named "filtered_unique"
addWorksheet(wb, "filtered_unique")

# Paste the filtered results into that sheet
writeData(wb, "filtered_unique", all_results_filtered)

# SAVE THE RESULTS

saveWorkbook(wb, output_excel, overwrite = TRUE)

# Remove the temporary file
file.remove("temp_reciprocal.fasta") # Cleanup


# 5. Check which IDs have not been used yet
# ------------------------------------------------------------------------------

# Get the names of the sheets in the .xlsx file
SheetNames <- getSheetNames(output_excel)


file <- output_excel

# Sheets to analyze
sheets_to_check <- SheetNames[2:(length(SheetNames)-1)]

# All obtained queries (used and unused)
all_query_sheet <- c()

# Check whether a BLASTp has already been run for that ID
 for (sheet in sheets_to_check) {
   
   # Store the sheet results in the variable results_sheet
   results_sheet <- read.xlsx(file, sheet = sheet) 
   
   # Store the sheet queries in query_sheet
   query_sheet <- results_sheet[, 2]
   
   # Store all queries in the variable all_query_sheet
   all_query_sheet <- c(all_query_sheet,query_sheet)
 }

# Keep only the unique values
all_query_sheet <- unique(all_query_sheet)

# Create a vector to store unique IDs
unique_id <- c() 

# Check which IDs have not been searched and store them in the vector unique_id
for (id in all_query_sheet) {
  if (!(id %in% queries_used)){
    unique_id <- c(unique_id,id)
  }
}
# 6. Run a BLASTp for the remaining IDs
#-------------------------------------------------------------------------------
for (id in unique_id) {

    queries_used <- c(queries_used,id)
    
    # Extract the specific sequence for the new query
    new_query_seq <- subject_db[id]
    writeXStringSet(new_query_seq, "temp_reciprocal.fasta")
    
    # Run BLAST for this specific subject
    reciprocal_out <- system2(
      command = blastp_executable,
      args = c("-query", "temp_reciprocal.fasta", "-subject", subject_fasta_file, "-outfmt", "6"),
      stdout = TRUE
    )
    
    # Safety check: Only continue if BLAST returns a result
    if (length(reciprocal_out) > 0) {
      
      # Convert raw text to a DataFrame
      res_table <- read.table(textConnection(reciprocal_out), sep = "\t")
      
      # Name the colums
      colnames(res_table) <- colnames(BLAST_AB)[1:12]
      
      # Length of the query sequence used in this reciprocal BLAST search
      query_length <- width(new_query_seq)
      
      # Add the query length
      res_table$QueryLength <- query_length
      
      # Compute the effective alignment length
      res_table$EffectiveAlignmentLength <- res_table$QueryEnd - res_table$QueryStart + 1
      
      # Compute coverage
      res_table$Coverage <- (res_table$EffectiveAlignmentLength / res_table$QueryLength) * 100
      
      # Apply the same filter
      res_table_filtered <- res_table %>%
        filter(
          Identity >= 30,
          Coverage >= 21,
          EValue < 1e-10
        )
      
      # Save only if results remain after filtering
      if (nrow(res_table_filtered) > 0) {
        
        # Add a column indicating the source search
        res_table_filtered$ReciprocalQuery <- id
        
        # Trim the protein name to 31 characters (Excel max)
        sheet_name <- substr(id, 1, 31)
        
        # Create the sheet only if it doesn't exist
        if (!(sheet_name %in% names(wb))) {
          addWorksheet(wb, sheet_name)
          writeData(wb, sheet_name, res_table_filtered)
        } else {
          message(paste("La hoja ya existe:", sheet_name))
        }
        
        # Add the data to the unified df
        all_results_filtered <- rbind(
          all_results_filtered,
          res_table_filtered[, c("Query", "Subject", "Identity", "EValue", "Coverage")]
        )
        message(paste("Processed BLAST for:", id))
        
      }
    }
}  

 

# 7. SAVE FINAL EXCEL FILE
# ------------------------------------------------------------------------------
# Create a table of unique results
unique_results <- all_results_filtered %>%
  select(Query, Subject, Identity, EValue, Coverage) %>% 
  distinct(Subject, .keep_all = TRUE) %>%         
  arrange(desc(Identity))                         

# Create a new sheet in the Excel workbook
addWorksheet(wb, "Unique_results")

# Copy the results
writeData(wb, "Unique_results", unique_results)

saveWorkbook(
  wb,
  output_excel,
  overwrite = TRUE
)

# Remove the temporary file
file.remove("temp_reciprocal.fasta") # Cleanup

# Mensaje final
print("Analysis complete! Check Blast_prot_local.xlsx")


