library(TCGAbiolinks)

# (Re-)create the query if needed
query <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Get file_id, file_name, md5sum
manifest_df <- getResults(
  query,
  cols = c("file_id", "file_name", "md5sum")
)

# Rename columns to GDC client format
colnames(manifest_df) <- c("id", "filename", "md5sum")

# Write out the manifest with tabs, no quotes
write.table(
  manifest_df,
  file      = "D:/MSC/gdc_manifest.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
