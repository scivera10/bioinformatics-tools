# Generar secuencias promotoras de genes oleosina
# Incluyendo la región 5'UTR pero excluyendo el CDS

library(tidyr)
library(dplyr)
library(readxl)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(Biostrings)

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
gene_table_xlsx   <- "gene_table.xlsx"         # Excel with gene IDs ('Entrez ID' column required)
gff_file          <- "annotation.gff"          # GFF/GFF3 genome annotation file
gff_filtered_CG8  <- "gr_gff_filtrado_CG8.gff3"  # Curated GFF3 (after manually adding missing genes)
genome_fasta      <- "genome.fna"              # Reference genome FASTA (indexed with samtools faidx)
upstream_bp       <- 1500                      # bp upstream of TSS to include in promoter
# ──────────────────────────────────────────────────────────────────────────────

# --------------------------------------------------
# 1. Leer la tabla con los genes de interés
# --------------------------------------------------

Comparativa_chr4_ncbi_vs_cast1 <- read_excel(
  gene_table_xlsx
)

LOC_ID <- Comparativa_chr4_ncbi_vs_cast1$`Entrez ID`

if (any(is.na(LOC_ID))) {
  message("Hay valores NA en los IDs seleccionados. Serán eliminados.")
  LOC_ID <- LOC_ID[!is.na(LOC_ID)]
}

# --------------------------------------------------
# 2. Importar el GFF completo
# --------------------------------------------------

gr_gff <- import(
  gff_file
)

# --------------------------------------------------
# 3. Filtrar genes de oleosinas
# --------------------------------------------------

gr_genes <- gr_gff[gr_gff$type == "gene" & gr_gff$gene %in% LOC_ID]

# Exportar por si quieres revisar manualmente
export(gr_genes, "gr_gff_filtrado.gff3", format = "gff3")

message("Se ha creado el archivo gr_gff_filtrado.gff3.")
message("Si necesitas añadir manualmente el gen CG8, hazlo y guarda el archivo como gr_gff_filtrado_CG8.gff3.")

# Si has añadido CG8 manualmente, vuelve a importar el archivo corregido
gr_genes <- import(
  gff_filtered_CG8
)

# --------------------------------------------------
# 4. Generar promotores iniciales
#    (1500 pb upstream + 0 pb downstream del TSS)
# --------------------------------------------------

gr_promoters <- promoters(gr_genes, upstream = upstream_bp, downstream = 0)

# --------------------------------------------------
# 5. Filtrar CDS solo de esos genes
# --------------------------------------------------

gr_cds <- gr_gff[gr_gff$type == "CDS" & gr_gff$gene %in% gr_genes$gene]

# Comprobación básica
if (length(gr_cds) == 0) {
  stop("No se encontraron regiones CDS para los genes seleccionados.")
}

# --------------------------------------------------
# 6. Resumir coordenadas CDS por gen
#    cds_min = coordenada mínima del CDS
#    cds_max = coordenada máxima del CDS
# --------------------------------------------------

cds_by_gene <- split(gr_cds, gr_cds$gene)

cds_min <- sapply(cds_by_gene, function(x) min(start(x)))
cds_max <- sapply(cds_by_gene, function(x) max(end(x)))

# Alinear con el orden de los promotores
genes_prom <- gr_promoters$gene
cds_min_vec <- cds_min[genes_prom]
cds_max_vec <- cds_max[genes_prom]

# Comprobar genes sin CDS asociado
if (any(is.na(cds_min_vec)) | any(is.na(cds_max_vec))) {
  genes_sin_cds <- genes_prom[is.na(cds_min_vec) | is.na(cds_max_vec)]
  message(
    paste(
      "Hay genes sin CDS asociado en el GFF:",
      paste(unique(genes_sin_cds), collapse = ", ")
    )
  )
}

# --------------------------------------------------
# 7. Recortar promotores para incluir 5'UTR
#    pero excluir la región codificante
# --------------------------------------------------

gr_promoters_final <- gr_promoters

# Convertir strand a vector normal
strand_vec <- as.character(strand(gr_promoters_final))

idx_plus  <- strand_vec == "+"
idx_minus <- strand_vec == "-"

# Genes con CDS válido
idx_valid <- !is.na(cds_min_vec) & !is.na(cds_max_vec)

# Hebra +
end(gr_promoters_final)[idx_plus & idx_valid] <- 
  cds_min_vec[idx_plus & idx_valid] - 1

# Hebra -
start(gr_promoters_final)[idx_minus & idx_valid] <- 
  cds_max_vec[idx_minus & idx_valid] + 1

# --------------------------------------------------
# 8. Comprobaciones
# --------------------------------------------------

# Eliminar posibles rangos inválidos si aparecieran
rangos_validos <- start(gr_promoters_final) <= end(gr_promoters_final)
if (!all(rangos_validos)) {
  message("Se eliminarán promotores con coordenadas inválidas tras el recorte.")
  gr_promoters_final <- gr_promoters_final[rangos_validos]
}

# Verificar que ya no haya solapamiento con CDS
solapan_cds <- overlapsAny(gr_promoters_final, gr_cds)
if (any(solapan_cds)) {
  warning("Algunos promotores finales todavía solapan con CDS.")
} else {
  message("Los promotores finales no solapan con CDS.")
}

# Tabla resumen opcional
resumen_promotores <- data.frame(
  gen = gr_promoters_final$gene,
  strand = as.character(strand(gr_promoters_final)),
  start_promotor_original = start(gr_promoters)[match(gr_promoters_final$gene, gr_promoters$gene)],
  end_promotor_original = end(gr_promoters)[match(gr_promoters_final$gene, gr_promoters$gene)],
  start_promotor_final = start(gr_promoters_final),
  end_promotor_final = end(gr_promoters_final),
  width_final = width(gr_promoters_final)
)

print(resumen_promotores)

# --------------------------------------------------
# 9. Extraer secuencias del genoma
# --------------------------------------------------

fa <- FaFile(
  genome_fasta
)

open(fa)

seqs <- getSeq(fa, gr_promoters_final)

names(seqs) <- gr_promoters_final$ID

writeXStringSet(seqs, "secuencias_promotores_sin_CDS_con_UTR.fa")

close(fa)

message("Se ha creado el archivo: secuencias_promotores_sin_CDS_con_UTR.fa")

# --------------------------------------------------
# 10. Preparar archivos para visualizar en IGV
# --------------------------------------------------

# Asegurar que los promotores tengan un nombre legible
if (is.null(gr_promoters_final$gene)) {
  stop("gr_promoters_final no tiene columna 'gene'")
}

# Crear nombres para las regiones promotoras
gr_promoters_final$promoter_id <- paste0(gr_promoters_final$gene, "_promoter")

# --------------------------------------------------
# 10A. Exportar GFF3
# --------------------------------------------------

export(gr_genes,
       "IGV_genes_oleosinas.gff3",
       format = "gff3")

export(gr_promoters_final,
       "IGV_promotores_sin_CDS_con_5UTR.gff3",
       format = "gff3")

export(gr_cds,
       "IGV_CDS_oleosinas.gff3",
       format = "gff3")

# --------------------------------------------------
# 10B. Exportar BED de genes
# --------------------------------------------------

bed_genes <- data.frame(
  seqnames = as.character(seqnames(gr_genes)),
  start = start(gr_genes) - 1,
  end = end(gr_genes),
  name = as.character(gr_genes$gene),
  score = 0,
  strand = as.character(strand(gr_genes)),
  stringsAsFactors = FALSE
)

write.table(
  bed_genes,
  file = "IGV_genes_oleosinas.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# --------------------------------------------------
# 10C. Exportar BED de promotores finales
# --------------------------------------------------

bed_promoters <- data.frame(
  seqnames = as.character(seqnames(gr_promoters_final)),
  start = start(gr_promoters_final) - 1,   # BED = 0-based
  end = end(gr_promoters_final),
  name = gr_promoters_final$promoter_id,
  score = 0,
  strand = as.character(strand(gr_promoters_final))
)

write.table(
  bed_promoters,
  file = "IGV_promotores_sin_CDS_con_5UTR.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# --------------------------------------------------
# 10D. Exportar BED de CDS
# --------------------------------------------------

bed_cds <- data.frame(
  seqnames = as.character(seqnames(gr_cds)),
  start = start(gr_cds) - 1,   # BED = 0-based
  end = end(gr_cds),
  name = gr_cds$gene,
  score = 0,
  strand = as.character(strand(gr_cds))
)

write.table(
  bed_cds,
  file = "IGV_CDS_oleosinas.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

message("Archivos para IGV generados correctamente.")
