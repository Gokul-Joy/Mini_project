#getting the TCGA data

# 🔧 STEP 1: Install and load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")  # Only once
library(TCGAbiolinks)


# 🎯 STEP 2: Query RNA-seq (HTSeq - Counts) for TCGA-PAAD
query_exp <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

# 💾 STEP 3: Download RNA-seq data
     # optional: reset download folder
GDCdownload(query_exp, method = "api", files.per.chunk = 20)


# 🧬 STEP 4: Prepare Expression Data
data_exp <- GDCprepare(query = query_exp)


if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
  BiocManager::install("SummarizedExperiment")

library(SummarizedExperiment)

# Now try
expression_matrix <- assay(data_exp)

sample_info <- colData(data_exp)

# Create labels (1 = tumor, 0 = normal)
labels <- ifelse(sample_info$shortLetterCode == "TP", 1, 0)



