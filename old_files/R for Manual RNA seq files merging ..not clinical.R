# Load required libraries
library(tidyverse)

# Set the path to the folder where all the STAR .tsv files are
data_dir <- "D:/MSC/Mini Project/AllCounts"  # Change this if needed
setwd(data_dir)

# List all .tsv files
files <- list.files(pattern = "\\.tsv$", full.names = TRUE)

# Custom function to read STAR .tsv files and extract gene_id and unstranded counts
read_star_counts <- function(file) {
  sample_id <- tools::file_path_sans_ext(basename(file))
  df <- read.delim(file, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
  
  # Only keep gene_id and unstranded count
  df <- df[, c("gene_id", "unstranded")]
  colnames(df) <- c("GeneID", sample_id)
  return(df)
}

# Read all files safely (skip files that fail)
count_list <- lapply(files, function(f) {
  tryCatch(read_star_counts(f), error = function(e) NULL)
})

# Remove failed reads
count_list_clean <- count_list[!sapply(count_list, is.null)]

# Merge all into one data frame
merged_counts <- reduce(count_list_clean, full_join, by = "GeneID")

# Filter out special rows like "__no_feature", etc., if they exist (yours might not have this)
merged_counts <- merged_counts[!grepl("^(N_|__)", merged_counts$GeneID), ]


# Set gene IDs as rownames
rownames(merged_counts) <- merged_counts$GeneID
merged_counts <- merged_counts[, -1]

# Convert to numeric safely (retain rownames)
merged_counts[] <- lapply(merged_counts, as.numeric)


head(merged_counts[,1:10])



#tumor vs normal labels

library(TCGAbiolinks)
library(dplyr)

# 1. Query RNA-seq metadata from TCGA-PAAD
query <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Extract file metadata table
metadata <- getResults(query)
head(metadata[, c("file_id", "cases.submitter_id", "sample_type")])

# Extract UUIDs from column names
column_uuids <- gsub("\\.rna_seq\\.augmented_star_gene_counts", "", colnames(merged_counts))

# Create a named vector of sample types for each UUID
sample_types <- metadata$sample_type
names(sample_types) <- metadata$file_id

# Assign tumor (1) or normal (0) labels
labels <- sample_types[column_uuids]
labels_binary <- ifelse(labels == "Primary Tumor", 1,
                        ifelse(labels == "Solid Tissue Normal", 0, NA))

# Optional: view a summary
table(labels_binary, useNA = "ifany")




















-----------------------------------
head(colnames(merged_counts), 5)



# 📦 Load library
library(TCGAbiolinks)

# 🔁 Reuse the same query (if still in memory)
# If not, rerun this:
query <- GDCquery(
   project = "TCGA-PAAD",
   data.category = "Transcriptome Profiling",
   data.type = "Gene Expression Quantification",
   workflow.type = "STAR - Counts"
 )

# ✅ Prepare expression matrix + sample metadata
# Run this to download to a specific folder
GDCdownload(query, 
            files.per.chunk = 20, 
            method = "client",  # uses GDC Data Transfer Tool
            directory = "D:/MSC/Mini Project/TCGA")


# Prepare the TCGA object from downloaded files

# Define the download + storage directory
download_dir <- "D:/MSC/Mini Project/GDCdata"

# Redownload using GDC client method (no gzip errors)
GDCdownload(query, method = "client", directory = download_dir)

# Now prepare it using the same directory
tcga_data <- GDCprepare(query, directory = download_dir)


