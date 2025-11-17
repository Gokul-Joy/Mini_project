#==============TCGA version=====================================================
#------------------------[ 1. Load Libraries ]------------------------
# Install Bioconductor packages if not already installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")




#BiocManager::install(c("recount", "DESeq2", "biomaRt", "org.Hs.eg.db"))

library(TCGAbiolinks)
library(recount)
library(DESeq2)
library(biomaRt)
library(org.Hs.eg.db)

# Load TCGA data (original step)
data <- readRDS("D:/MSC/MiniProject/TCGA_PAAD_expr.rds")
expr_tcga <- assay(data)        # genes x samples
pheno_tcga <- colData(data)

View(expr_tcga)

# Filter TCGA for protein-coding genes
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl_ids <- gsub("\\..*","",rownames(expr_tcga))
gene_info <- getBM(attributes=c("ensembl_gene_id","gene_biotype","external_gene_name"),
                   filters="ensembl_gene_id",
                   values=ensembl_ids,
                   mart=ensembl)
protein_coding_ids <- gene_info$ensembl_gene_id[gene_info$gene_biotype=="protein_coding"]
keep_idx <- ensembl_ids %in% protein_coding_ids
expr_tcga <- expr_tcga[keep_idx, ]

# Map Ensembl IDs → gene symbols
gene_symbols <- gene_info$external_gene_name[match(ensembl_ids[keep_idx],
                                                   gene_info$ensembl_gene_id)]
valid_idx <- !is.na(gene_symbols)
expr_tcga <- expr_tcga[valid_idx, ]
rownames(expr_tcga) <- gene_symbols[valid_idx]

# Filter sample types in TCGA
sample_types <- substr(pheno_tcga$barcode, 14, 15)
labels <- ifelse(sample_types == "01", "Tumor", "Normal")
valid_samples <- sample_types %in% c("01", "11")
expr_tcga <- expr_tcga[, valid_samples]
labels <- labels[valid_samples]
pheno_tcga <- pheno_tcga[valid_samples, ]

View(expr_tcga)
dim(expr_tcga)

# tcga_qc_plots.R
library(DESeq2); library(ggplot2)
# vst normalization for visualization
dds_tmp <- DESeqDataSetFromMatrix(countData = expr_tcga, 
                                  colData = data.frame(condition=factor(labels)), design=~condition)
vst_tcga <- vst(dds_tmp, blind=TRUE)
mat_vst <- assay(vst_tcga)

# density plot
df <- data.frame(Expression = as.vector(mat_vst))
ggplot(df, aes(x = Expression)) + geom_density(alpha=0.5) +
  labs(title = "VST density of TCGA (processed)") + theme_minimal()






# ------------------------------
# STEP 1: Process GTEx data (your code)

library(recount3)

# Download GTEx pancreas data as RangedSummarizedExperiment
#rse_pancreas <- create_rse_manual(
# project = "PANCREAS",
#project_home = "data_sources/gtex",
#organism = "human",
#annotation = "gencode_v26",
#type = "gene"
#)

# Extract raw counts matrix

#expr_counts <- assay(rse_pancreas, "raw_counts")

# Optionally save for later use
#saveRDS(rse_pancreas, "gtex_pancreas_rse.rds")
#saveRDS(expr_counts, "D:/MSC/MiniProject/gtex_pancreas_raw_counts.rds")
expr_counts <- readRDS("D:/MSC/MiniProject/gtex_pancreas_raw_counts.rds")

gtex_ensembl_ids <- gsub("\\..*","",rownames(expr_counts))
# Get gene info from same biomaRt (make sure using same Ensembl version)
gene_info_gtex <- getBM(attributes=c("ensembl_gene_id","gene_biotype","external_gene_name"),
                        filters="ensembl_gene_id",
                        values=gtex_ensembl_ids,
                        mart=ensembl)
protein_coding_gtex <- gene_info_gtex$ensembl_gene_id[gene_info_gtex$gene_biotype=="protein_coding"]
keep_idx_gtex <- gtex_ensembl_ids %in% protein_coding_gtex

expr_gtex_filtered <- expr_counts[keep_idx_gtex, ]
gene_symbols_gtex <- gene_info_gtex$external_gene_name[match(gtex_ensembl_ids[keep_idx_gtex], gene_info_gtex$ensembl_gene_id)]
valid_gtex <- !is.na(gene_symbols_gtex)
expr_gtex_filtered <- expr_gtex_filtered[valid_gtex, ]
rownames(expr_gtex_filtered) <- gene_symbols_gtex[valid_gtex]
# Remove duplicate gene symbols if any
expr_gtex_filtered <- expr_gtex_filtered[!duplicated(rownames(expr_gtex_filtered)), ]



#===============================================================================

library(WGCNA)

# Your current TCGA expression matrix with duplicates
expr_tcga_matrix <- expr_tcga  # assuming your current expr_tcga object

# Clean up: Remove rows with blank or NA gene symbols (rownames) or Ensembl IDs
rowGroup <- rownames(expr_tcga_matrix)                   # gene symbols
rowID_raw <- rownames(expr_tcga_matrix)                  # original rownames with versions

# Extract Ensembl IDs by removing version suffix (e.g. ENSG000001.1 -> ENSG000001)
rowID <- gsub("\\..*", "", rowID_raw)

# Identify rows with blanks or NAs in gene symbols or Ensembl IDs
valid_rows <- !(rowGroup == "" | is.na(rowGroup) | rowID == "" | is.na(rowID))
expr_tcga_clean <- expr_tcga_matrix[valid_rows, ]
rowGroup_clean <- rownames(expr_tcga_clean)
rowID_clean_raw <- rownames(expr_tcga_clean)
rowID_clean <- gsub("\\..*", "", rowID_clean_raw)

# Ensure unique rowID by appending suffixes to duplicates if any remain
dup_flags <- duplicated(rowID_clean) | duplicated(rowID_clean, fromLast = TRUE)
rowID_unique <- rowID_clean
rowID_unique[dup_flags] <- make.unique(rowID_clean[dup_flags])

# Collapse duplicates by maxRowVariance using cleaned and unique IDs
collapsed <- collapseRows(
  datET = expr_tcga_clean,
  rowGroup = rowGroup_clean,
  rowID = rowID_unique,
  method = "maxRowVariance"
)

# Extract collapsed matrix with unique gene symbols
expr_tcga_collapsed <- collapsed$datETcollapsed

# Sanity checks
cat("Dimensions after collapsing:", dim(expr_tcga_collapsed), "\n")
cat("Duplicates after collapsing:", sum(duplicated(rownames(expr_tcga_collapsed))), "\n")

# Replace original expr_tcga with collapsed unique-gene matrix
expr_tcga <- expr_tcga_collapsed
# Recompute common genes after collapsing expr_tcga
common_genes <- intersect(rownames(expr_tcga_collapsed), rownames(expr_gtex_filtered))

# Subset both matrices to shared genes
expr_tcga_sub <- expr_tcga_collapsed[common_genes, , drop=FALSE]
expr_gtex_sub <- expr_gtex_filtered[common_genes, , drop=FALSE]

#===============================================================================




# Subset both to common genes
expr_tcga_sub <- expr_tcga[common_genes, ]
expr_gtex_sub <- expr_gtex_filtered[common_genes, ]
#=============================================================================
#=============================================================================
# ================= PCA-based GTEx sample selection (178 samples) =================
# ------------------------------
# PCA-based reduction of GTEx Normal samples to target N = 178
# Place this block immediately after:
# expr_gtex_sub <- expr_gtex_filtered[common_genes, , drop=FALSE]
# ------------------------------

target_n <- 174

if (ncol(expr_gtex_sub) > target_n) {
  set.seed(42)
  message("GTEx samples before reduction: ", ncol(expr_gtex_sub))
  
  libSizes <- colSums(expr_gtex_sub)
  libSizes[libSizes == 0] <- 1
  norm_counts <- t( t(expr_gtex_sub) / libSizes ) * median(libSizes)
  log2_norm <- log2(norm_counts + 1)
  
  # Remove samples with zero variance
  var_per_sample <- apply(log2_norm, 2, var)
  log2_norm <- log2_norm[, var_per_sample > 0, drop = FALSE]
  
  # ✅ Remove genes with zero variance
  var_per_gene <- apply(log2_norm, 1, var)
  log2_norm <- log2_norm[var_per_gene > 0, , drop = FALSE]
  
  # PCA
  pca <- prcomp(t(log2_norm), center = TRUE, scale. = TRUE)
  
  pcs12 <- pca$x[, 1:2, drop = FALSE]
  centroid <- colMeans(pcs12)
  dists_to_centroid <- sqrt(rowSums((pcs12 - matrix(centroid, nrow = nrow(pcs12), ncol = 2, byrow = TRUE))^2))
  
  sel_order <- order(dists_to_centroid)
  selected_samples <- rownames(pcs12)[ sel_order[1:target_n] ]
  
  expr_gtex_sub <- expr_gtex_sub[, selected_samples, drop = FALSE]
  message("GTEx samples after reduction: ", ncol(expr_gtex_sub))
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  plot(pca$x[,1], pca$x[,2], pch = 20, xlab = "PC1", ylab = "PC2",
       main = paste0("GTEx PCA (", ncol(expr_gtex_sub), " samples selected)"))
  points(pca$x[selected_samples, 1], pca$x[selected_samples, 2], pch = 20, cex = 1.2, col = "red")
}



#=============================================================================
#=============================================================================


# ------------------------------
# STEP 3: Combine the counts
expr_combined <- cbind(expr_tcga_sub, expr_gtex_sub)

# Create combined labels: TCGA tumor/normal and GTEx normals
tcga_labels <- factor(ifelse(labels=="Tumor", "Tumor", "Normal"), levels=c("Normal","Tumor"))
gtex_labels <- factor(rep("Normal", ncol(expr_gtex_sub)), levels=c("Normal","Tumor"))
combined_labels <- factor(c(as.character(tcga_labels), as.character(gtex_labels)),
                          levels=c("Normal","Tumor"))

coldata_combined <- data.frame(condition=combined_labels)
rownames(coldata_combined) <- colnames(expr_combined)


table(combined_labels)



# ------------------------------
# STEP 4: Run DESeq2 on the combined dataset
dds <- DESeqDataSetFromMatrix(countData=expr_combined,
                              colData=coldata_combined,
                              design=~condition)

dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","Tumor","Normal"))
res <- lfcShrink(dds, coef="condition_Tumor_vs_Normal", res=res, type="ashr")

# Save for downstream
deg_tcga <- as.data.frame(res[order(res$padj), ])

table(combined_labels)
#------------------------[ 5. Extract and Filter DEGs ]------------------------
deg_tcga_sig_1 <- deg_tcga[deg_tcga$padj < 0.05 & abs(deg_tcga$log2FoldChange) > 1, ]


# Convenience copies
deg_deseq2_1 <- deg_tcga_sig_1


cat("Significant DEGs (logFC>1):", nrow(deg_deseq2_1), "\n")


#------------------------[ 6. Map Gene Symbols (safe) ]------------------------
symbols_1 <- rownames(deg_tcga_sig_1)  # already gene symbols
# already gene symbols

genes_tcga_1 <- unname(symbols_1)


#===============================================================================
#===============================================================================




#install.packages("xml2")



# Load all libraries,almost ig
library(GEOquery); library(limma); library(DESeq2); library(TCGAbiolinks)
library(SummarizedExperiment); library(clusterProfiler); library(org.Hs.eg.db)
library(enrichplot); library(DOSE)



library(WGCNA)


# Load GSE62452
gse1 <- getGEO(filename = 'D:/MSC/MiniProject/Dataset/GSE62452_series_matrix.txt.gz')
expr_geo1 <- exprs(gse1)
pheno1 <- pData(gse1)
labels_geo1 <- ifelse(pheno1$`tissue:ch1` == "Pancreatic tumor", 1, 0)

# Load GSE183795
gse2 <- getGEO(filename = 'D:/MSC/MiniProject/Dataset/GSE183795_series_matrix.txt.gz')
expr_geo2 <- exprs(gse2)
pheno2 <- pData(gse2)

# Function to detect tumor/normal labels robustly
detect_label_column_and_assign <- function(pheno_df) {
  candidates <- c("tissue:ch1", "source_name_ch1", "characteristics_ch1", "tissue", "title")
  found <- intersect(candidates, colnames(pheno_df))
  if(length(found) > 0) {
    colvals <- as.character(pheno_df[[found[1]]])
  } else {
    idx <- grep("tissue|tumor|sample type|characteristics", tolower(colnames(pheno_df)))
    if(length(idx)>0) colvals <- as.character(pheno_df[[idx[1]]]) else colvals <- apply(pheno_df, 1, paste, collapse=" ")
  }
  is_tumor <- grepl("tumor|cancer|PDAC", tolower(colvals))
  is_normal <- grepl("normal|healthy|donor", tolower(colvals))
  labels <- ifelse(is_tumor, 1, ifelse(is_normal, 0, NA))
  return(labels)
}

labels_geo2 <- detect_label_column_and_assign(pheno2)

# Probe to gene mapping (platform assumed same for both)
gpl <- getGEO("GPL6244", AnnotGPL = TRUE)
gpl_table <- Table(gpl)
probe2gene <- gpl_table[, c("ID", "Gene symbol")]
colnames(probe2gene) <- c("PROBEID", "SYMBOL")
probe2gene$SYMBOL <- sapply(strsplit(probe2gene$SYMBOL, "///", fixed = TRUE), `[`, 1)
probe2gene <- probe2gene[!is.na(probe2gene$SYMBOL) & probe2gene$SYMBOL != "", ]

# Annotate and collapse: GSE62452
expr_geo1_df <- data.frame(PROBEID = rownames(expr_geo1), expr_geo1)
expr_geo1_annot <- merge(probe2gene, expr_geo1_df, by = "PROBEID")
expr_geo1_matrix <- as.matrix(expr_geo1_annot[, -(1:2)])
rowGroup1 <- expr_geo1_annot$SYMBOL
rowID1 <- expr_geo1_annot$PROBEID
collapsed1 <- collapseRows(datET = expr_geo1_matrix, rowGroup = rowGroup1, rowID = rowID1, method = "maxRowVariance")
expr_geo1_maxvar <- collapsed1$datETcollapsed

# Annotate and collapse: GSE183795
expr_geo2_df <- data.frame(PROBEID = rownames(expr_geo2), expr_geo2)
expr_geo2_annot <- merge(probe2gene, expr_geo2_df, by = "PROBEID")
expr_geo2_matrix <- as.matrix(expr_geo2_annot[, -(1:2)])
rowGroup2 <- expr_geo2_annot$SYMBOL
rowID2 <- expr_geo2_annot$PROBEID
collapsed2 <- collapseRows(datET = expr_geo2_matrix, rowGroup = rowGroup2, rowID = rowID2, method = "maxRowVariance")
expr_geo2_maxvar <- collapsed2$datETcollapsed

# Log2 transform if needed (heuristic check)
ensure_log2 <- function(mat) {
  q <- quantile(mat, probs = c(0.99))
  if(q[1] > 100) mat <- log2(mat + 1)
  return(mat)
}
expr_geo1_maxvar <- ensure_log2(expr_geo1_maxvar)
expr_geo2_maxvar <- ensure_log2(expr_geo2_maxvar)

# Keep only genes common to both
common_genes <- intersect(rownames(expr_geo1_maxvar), rownames(expr_geo2_maxvar))
e1 <- expr_geo1_maxvar[common_genes, , drop = FALSE]
e2 <- expr_geo2_maxvar[common_genes, , drop = FALSE]

# Combine expression and labels
expr_geo_maxvar <- cbind(e1, e2)
labels_geo <- c(labels_geo1, labels_geo2)
geo_batch <- c(rep("GSE62452", ncol(e1)), rep("GSE183795", ncol(e2)))

# Quantile normalization and batch correction
expr_geo_maxvar <- normalizeBetweenArrays(expr_geo_maxvar, method = "quantile")
modcombat <- model.matrix(~as.factor(labels_geo))
expr_geo_maxvar <- ComBat(dat = expr_geo_maxvar, batch = geo_batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)

# deg_geo is produced exactly as before (limma step using expr_geo_maxvar and labels_geo)
group_geo <- factor(labels_geo)
design_geo <- model.matrix(~group_geo)
fit_geo <- lmFit(expr_geo_maxvar, design_geo)
fit_geo <- eBayes(fit_geo)
deg_geo <- topTable(fit_geo, coef = 2, number = Inf, adjust.method = "fdr")
deg_geo_sig_1 <- deg_geo[deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 1, ]
gene_geo_1 <- rownames(deg_geo_sig_1)



#======================Block of common genes and GO nd KEGG=====================
# CHANGED: use gene lists created earlier (genes_tcga_1, genes_tcga_0.05) which we defined from deg_tcga
common_genes_1 <- intersect(gene_geo_1, genes_tcga_1)
length(gene_geo_1); length(genes_tcga_1);
common_genes_0.05 <- intersect(gene_geo_1, genes_tcga_1)
length(gene_geo_0.05); length(genes_tcga_1);

length(common_genes_1)






library(ggplot2)

# --- TCGA Volcano Plot ---
# Assumes deg_tcga is your DESeq2 result data.frame and rownames are gene symbols

deg_tcga$threshold <- as.factor(deg_tcga$padj < 0.05 & abs(deg_tcga$log2FoldChange) > 1)

ggplot(deg_tcga, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("gray", "orange")) +
  theme_minimal() +
  labs(title = "TCGA Volcano Plot (Tumor vs Normal)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")

# --- GEO Volcano Plot ---
# Assumes deg_geo is your limma topTable data.frame and rownames are gene symbols

deg_geo$threshold <- as.factor(deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 0.05)

ggplot(deg_geo, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("gray", "coral")) +
  theme_minimal() +
  labs(title = "GEO Volcano Plot (Tumor vs Normal)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")











































#|
#|
#|
#|
entrez_ids_1 <- bitr(common_genes_1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

#|
#|
# GO
ego_1 <- enrichGO(gene = entrez_ids_1$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                  ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
#|
#|
# KEGG
ekegg_1 <- enrichKEGG(gene = entrez_ids_1$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)

#|
#|
# Plots----------
#barplot(ego_0.05, showCategory = 15, title = "GO Enrichment")
#dotplot(ekegg_0.05, showCategory = 15, title = "KEGG Pathway Enrichment")
#===============================================================================








#------------------------[ Optional GSEA (TCGA only) ]------------------------
# CHANGED: ensure deg_tcga exists before using it
gene_stats_1<- deg_tcga$log2FoldChange   # CHANGED: adapted to deg_tcga frame (note: original code used deg_tcga$logFC — here use log2FoldChange)
names(gene_stats_1) <- rownames(deg_tcga)
#|

head(rownames(deg_tcga), 10)

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)  # optional, adds GSEA plotting functions

# Suppose your vector looks like this:
# gene_stats_1 <- setNames(deg_tcga$logFC, rownames(deg_tcga))

# 1️⃣ Convert SYMBOL → ENTREZ IDs
gene_symbols <- names(gene_stats_1)

# Map to Entrez IDs
symbol2entrez <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# 2️⃣ Merge back to retain stats for mapped genes
gene_stats_entrez <- gene_stats_1[symbol2entrez$SYMBOL]
names(gene_stats_entrez) <- symbol2entrez$ENTREZID

# Remove NAs
gene_stats_entrez <- gene_stats_entrez[!is.na(names(gene_stats_entrez))]

# 3️⃣ Run GSEA KEGG
gsea_kegg_1 <- gseKEGG(
  geneList = sort(gene_stats_entrez, decreasing = TRUE),
  organism = "hsa",
  pvalueCutoff = 0.2
)

# 4️⃣ View results
head(gsea_kegg_1@result)




dotplot(gsea_kegg_1, showCategory = 15, title = "KEGG GSEA logfc1 (TCGA-PAAD)")










#===============================================================================
# Clinical matching - TCGA/TCGAbiolinks
#===============================================================================



# Download TCGA clinical data for PAAD
clin_df <- GDCquery_clinic(project = "TCGA-PAAD", type = "clinical")

# Standardize sample barcodes (first 12 chars is participant)
clin_df$SampleID <- substr(clin_df$bcr_patient_barcode, 1, 12)

cat("Cleaned clinical SampleIDs:\n")
print(head(clin_df$SampleID))

# Extract base sample ids from expression matrix
expr_ids <- substr(colnames(expr_tcga), 1, 12)
cat("Cleaned expr IDs:\n")
print(head(expr_ids))

# Match
matched_ids <- intersect(expr_ids, clin_df$SampleID)
cat("✅ Matched samples: ", length(matched_ids), "\n")

# Subset and reorder
expr_tcga <- expr_tcga[, expr_ids %in% matched_ids]
clin_df <- clin_df[clin_df$SampleID %in% matched_ids, ]
clin_df <- clin_df[match(substr(colnames(expr_tcga), 1, 12), clin_df$SampleID), ]

# Coercion (update with correct TCGAbiolinks field names)
clin_df$vital_status <- as.numeric(clin_df$vital_status == "Dead")
clin_df$days_to_death <- as.numeric(clin_df$days_to_death)
clin_df$days_to_last_follow_up <- as.numeric(clin_df$days_to_last_follow_up)

# Survival objects
surv_time <- ifelse(is.na(clin_df$days_to_death), clin_df$days_to_last_follow_up, clin_df$days_to_death)
surv_status <- clin_df$vital_status

summary(surv_time)
table(surv_status)
dim(expr_tcga)







#===============================================================================
#LASSO + Cox Regression Modeling 
library(org.Hs.eg.db)

head(rownames(expr_tcga))
# Convert Entrez IDs in expr_tcga rownames to gene symbols
#symbols <- mapIds(org.Hs.eg.db,
#    keys = rownames(expr_tcga),   # Entrez IDs
#   column = "SYMBOL",
#  keytype = "ENTREZID",
#  multiVals = "first")

# Remove rows with NA symbols
#valid_idx <- !is.na(symbols)
#expr_tcga <- expr_tcga[valid_idx, ]
#symbols <- symbols[valid_idx]

# Assign gene symbols as rownames
#rownames(expr_tcga) <- symbols

# Remove duplicated gene symbols if any
#expr_tcga <- expr_tcga[!duplicated(rownames(expr_tcga)), ]

cat("New expression matrix rownames (symbols):\n")
print(head(rownames(expr_tcga)))

# Filter expression matrix to common genes (gene symbols)
expr_common_1 <- expr_tcga[rownames(expr_tcga) %in% common_genes_1, ]

cat("Expression matrix after gene filtering:\n")
print(dim(expr_common_1))

#===============================================================================
#|
#|
#|
#|
# Create survival objec=========================================================
library(survival)
library(glmnet)

#=============
#quick fix for surv time less than zero
#=============


# Remove samples with time <= 0 or NA
valid_samples <- !is.na(surv_time) & surv_time > 0 & !is.na(surv_status)

# Apply filter to everything
expr_common_1 <- expr_common_1[, valid_samples]
surv_time <- surv_time[valid_samples]
surv_status <- surv_status[valid_samples]


# Transpose expression (samples × genes)
x_1 <- t(expr_common_1)


# Create survival object
y <- Surv(time = surv_time, event = surv_status)

# Check dimensions again
cat("Final dimensions:\n")
print(dim(x_1))
#///////////////////////////////////////////////////////////////////////////////
n_runs <- 100
seed_base <- 123
alpha_val <- 1
nfolds_val <- 10
stability_cutoff <- 0.7  # keep genes appearing in >= 70% of runs

gene_list <- list()

# ===============================
# 100-iteration LASSO stability selection
# ===============================
for (i in 1:n_runs) {
  set.seed(seed_base + i)  # Different seed for variability
  
  fit <- cv.glmnet(
    x_1, y,
    family = "cox",
    alpha = alpha_val,
    nfolds = nfolds_val
  )
  
  coef_opt <- coef(fit, s = "lambda.min")
  selected_genes <- rownames(coef_opt)[which(coef_opt != 0)]
  selected_genes <- selected_genes[selected_genes != "(Intercept)"]  # remove intercept
  
  gene_list[[i]] <- selected_genes
}


#============================================================

gene_freq <- sort(table(unlist(gene_list)), decreasing = TRUE)
gene_stability <- as.numeric(gene_freq) / n_runs   # CHANGED: keep single canonical definition with names
names(gene_stability) <- names(gene_freq)

# View stable genes (≥ 70% of runs)
stable_genes <- gene_stability[gene_stability >= 0.7]
print(stable_genes)


library(tibble)
# or
library(dplyr)






# Plot stability
# #CHANGED: keep gene_stability as numeric vector with names so this works downstream
gene_df <- as.data.frame(gene_stability) %>%
  rownames_to_column("Gene") %>%
  rename(Stability = gene_stability)

ggplot(gene_df, aes(x = reorder(Gene, Stability), y = Stability)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = paste("Gene Stability over", n_runs, "LASSO CV runs"),
    x = "Gene",
    y = "Stability (Proportion of Runs)"
  ) +
  theme_minimal()

cat("✅ Stable genes selected (≥", stability_cutoff*100, "% runs):\n")
print(stable_genes)

# gene_stability <- gene_freq / n_runs   # CHANGED: removed duplicate re-assignment (was redundant)
#-----------------------------------------------------------------------


stable_gene_names <- names(stable_genes)
x_stable <- x_1[, stable_gene_names, drop = FALSE]
length(stable_genes) 
length(x_stable)
dim(x_stable) 









fit_final <- cv.glmnet(
  x_stable, y,
  family = "cox",
  alpha = alpha_val,
  nfolds = nfolds_val
)



lasso_coef <- coef(fit_final, s = fit_final$lambda.min)
lasso_genes <- rownames(lasso_coef)[lasso_coef[, 1] != 0]



coef_vector <- as.vector(lasso_coef[lasso_coef[, 1] != 0])
names(coef_vector) <- lasso_genes

# Risk scores
risk_score <- as.numeric(x_stable[, lasso_genes, drop = FALSE] %*% coef_vector)
median_cutoff <- median(risk_score)
risk_group <- ifelse(risk_score > median_cutoff, "High", "Low")
risk_group <- factor(risk_group, levels = c("Low", "High"))

# ===============================
# Kaplan-Meier Plot
# ===============================
library(survival)
library(survminer)

surv_obj <- Surv(surv_time, surv_status)
fit_km <- survfit(surv_obj ~ risk_group)

ggsurvplot(fit_km,
           data = data.frame(risk_group),
           pval = TRUE,
           risk.table = TRUE,
           title = "Survival Curve: High vs Low Risk (Stable LASSO Genes)",
           palette = c("blue", "red"))

# ===============================
# Time-dependent ROC
# ===============================
library(timeROC)

roc_obj <- timeROC(
  T = surv_time,
  delta = surv_status,
  marker = risk_score,
  cause = 1,
  times = c(365, 1095, 1825),
  iid = TRUE
)

auc_values <- roc_obj$AUC
names(auc_values) <- c("1yr", "3yr", "5yr")

plot(roc_obj, time = 365, col = "blue", title = TRUE)
plot(roc_obj, time = 1095, add = TRUE, col = "green")
plot(roc_obj, time = 1825, add = TRUE, col = "red")

legend("bottomright",
       legend = paste(names(auc_values), "AUC=", round(auc_values, 2)),
       col = c("blue", "green", "red"),
       lwd = 2)


cat("✅ Stable genes selected (≥", stability_cutoff*100, "% runs):\n")
print(stable_genes)
#===============================================================================
#===============================================================================
#===============================================================================




# slide1_datasets_counts.R
library(ggplot2)

df_counts <- data.frame(
  dataset = c("TCGA Tumor", "TCGA Normal", "GTEx Pancreas", "GEO GSE62452"),
  count = c(sum(combined_labels == "Tumor"),
            sum(combined_labels == "Normal") - ncol(expr_gtex_sub), # approximate TCGA normal count
            ncol(expr_gtex_sub),
            ncol(expr_geo_maxvar))
)

ggplot(df_counts, aes(x = dataset, y = count, fill = dataset)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = count), vjust = -0.25) +
  labs(title = "Samples by dataset", x = NULL, y = "Number of samples") +
  theme_minimal()




# slide2_tcga_density.R
library(ggplot2)

# use expr_tcga (collapsed, gene symbols). Log-transform for visualization.
vals <- as.vector(log2(expr_tcga + 1))
ggplot(data.frame(Expression = vals), aes(x = Expression)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of TCGA expression (log2(count+1))", x = "log2(count+1)", y = "Density") +
  theme_minimal()


# slide3_combined_boxplot.R

# -----------------------------
# Check total dimensions
# -----------------------------
cat("Combined expression matrix dimensions (Genes x Samples):\n")
dim(expr_combined)


# -----------------------------
# Define combined sample labels
# -----------------------------
samples_combined <- colnames(expr_combined)

# TCGA Tumor/Normal already defined in your 'labels' vector (same order as expr_tcga columns)
# Identify which columns belong to TCGA and which to GTEx
tcga_samples <- colnames(expr_tcga)
gtex_samples <- colnames(expr_gtex_sub)

# Match labels for TCGA part and add GTEx as "Normal"
labels_combined <- c(labels, rep("Normal", length(gtex_samples)))

# Double-check lengths match number of samples
cat("\nLength check (should match number of columns in expr_combined):\n")
length(labels_combined)

# -----------------------------
# Summarize Tumor vs Normal counts
# -----------------------------
cat("\nTumor vs Normal sample counts (after combining TCGA + GTEx):\n")
table(labels_combined)

# -----------------------------
# Optional: quick barplot visualization
# -----------------------------
barplot(table(labels_combined),
        col = c("tomato", "skyblue"),
        main = "Sample distribution after combining TCGA + GTEx",
        ylab = "Number of samples")

