#==============TCGA version=====================================================
#------------------------[ 1. Load Libraries ]------------------------
# Install Bioconductor packages if not already installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
<<<<<<< HEAD
 # install.packages("BiocManager")
=======
  install.packages("BiocManager")
>>>>>>> a7b5f766fde44ca6cf206a001507ceee8e2e5004

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
<<<<<<< HEAD
#expr_counts <- assay(rse_pancreas, "raw_counts")
=======
expr_counts <- assay(rse_pancreas, "raw_counts")
>>>>>>> a7b5f766fde44ca6cf206a001507ceee8e2e5004

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
deg_tcga_sig_1.5 <- deg_tcga[deg_tcga$padj < 0.05 & abs(deg_tcga$log2FoldChange) > 1, ]
deg_tcga_sig_0.5 <- deg_tcga[deg_tcga$padj < 0.05 & abs(deg_tcga$log2FoldChange) > 0.5, ]
deg_tcga_sig_0.05 <- deg_tcga[deg_tcga$padj < 0.05 & abs(deg_tcga$log2FoldChange) > 0.05, ]

# Convenience copies
deg_deseq2_1 <- deg_tcga_sig_1.5
deg_deseq2_0.5 <- deg_tcga_sig_0.5
deg_deseq2_0.05 <- deg_tcga_sig_0.05

cat("Significant DEGs (logFC>1):", nrow(deg_deseq2_1), "\n")
cat("Significant DEGs (logFC>0.5):", nrow(deg_deseq2_0.5), "\n")
cat("Significant DEGs (logFC>0.05):", nrow(deg_deseq2_0.05), "\n")

#------------------------[ 6. Map Gene Symbols (safe) ]------------------------
symbols_1.5 <- rownames(deg_tcga_sig_1.5)  # already gene symbols
symbols_0.05 <- rownames(deg_tcga_sig_0.05) # already gene symbols

genes_tcga_1.5 <- unname(symbols_1.5)
genes_tcga_0.05 <- unname(symbols_0.05)

#===============================================================================
#===============================================================================




#install.packages("xml2")



# Load all libraries,almost ig
library(GEOquery); library(limma); library(DESeq2); library(TCGAbiolinks)
library(SummarizedExperiment); library(clusterProfiler); library(org.Hs.eg.db)
library(enrichplot); library(DOSE)







#------------------------[ 2. GEO]------------------------
#------------------------[ 1. Load GEO Data ]------------------------
gse <- getGEO(filename = 'D:/MSC/MiniProject/Dataset/GSE62452_series_matrix.txt.gz')
expr_geo <- exprs(gse)
pheno <- pData(gse)

#------------------------[ 2. Sample Labels ]------------------------
labels_geo <- ifelse(pheno$`tissue:ch1` == "Pancreatic tumor", 1, 0)
stopifnot(length(labels_geo) == ncol(expr_geo))  # Check label/sample match

#------------------------[ 3. Probe to Gene Mapping ]------------------------
gpl <- getGEO("GPL6244", AnnotGPL = TRUE)
gpl_table <- Table(gpl)
probe2gene <- gpl_table[, c("ID", "Gene symbol")]
colnames(probe2gene) <- c("PROBEID", "SYMBOL")
probe2gene$SYMBOL <- sapply(strsplit(probe2gene$SYMBOL, "///", fixed = TRUE), `[`, 1)
probe2gene <- probe2gene[!is.na(probe2gene$SYMBOL) & probe2gene$SYMBOL != "", ]

expr_geo_df <- data.frame(PROBEID = rownames(expr_geo), expr_geo)
expr_geo_annot <- merge(probe2gene, expr_geo_df, by = "PROBEID")

#------------------------[ 4. Collapse Probes to Genes ]------------------------
library(WGCNA)
expr_geo_matrix <- as.matrix(expr_geo_annot[, -(1:2)])  # CHANGED: renamed so we don't overwrite TCGA expr
rowGroup <- expr_geo_annot$SYMBOL
rowID <- expr_geo_annot$PROBEID

collapsed <- collapseRows(datET = expr_geo_matrix,
                          rowGroup = rowGroup,
                          rowID = rowID,
                          method = "maxRowVariance")
expr_geo_maxvar <- collapsed$datETcollapsed

#------------------------[ 5. Differential Expression Analysis ]------------------------
library(limma)
group_geo <- factor(labels_geo)
design_geo <- model.matrix(~group_geo)
fit_geo <- lmFit(expr_geo_maxvar, design_geo)
fit_geo <- eBayes(fit_geo)
deg_geo <- topTable(fit_geo, coef = 2, number = Inf, adjust.method = "fdr")
dim(deg_geo)
#------------------------[ 6. Filter Significant DEGs ]------------------------
deg_geo_sig_1.5 <- deg_geo[deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 1, ]
deg_geo_sig_0.5 <- deg_geo[deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 0.5, ]
deg_geo_sig_0.05 <- deg_geo[deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 0.05, ]

cat("Number of significant DEGs 1.5:", nrow(deg_geo_sig_1.5), "\n")
cat("Number of significant DEGs 0.5:", nrow(deg_geo_sig_0.5), "\n")
cat("Number of significant DEGs 0.05:", nrow(deg_geo_sig_0.05), "\n")



# For TCGA (you used expr_tcga and labels)
table(combined_labels)          # tumor vs normal count

# For GEO (assuming labels_geo and expr_geo exist)
table(labels_geo)
head(rownames(deg_geo_sig_1.5))
head(probe2gene$PROBEID)
head(rownames(deg_geo))



# No re-mapping needed! Row names are already gene symbols.
gene_geo_1.5 <- rownames(deg_geo_sig_1.5)
gene_geo_0.5 <- rownames(deg_geo_sig_0.5)
gene_geo_0.05 <- rownames(deg_geo_sig_0.05)
cat("Number of significant DEGs 1.5:", length(gene_geo_1.5), "\n")
cat("Number of significant DEGs 0.5:", length(gene_geo_0.5), "\n")
cat("Number of significant DEGs 0.05:", length(gene_geo_0.05), "\n")


#------------------------[ 8. Plot Density of Collapsed Expression ]------------------------
library(ggplot2)
ggplot(data.frame(Expression = as.vector(expr_geo_maxvar)),
       aes(x = Expression)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Density of Collapsed Gene Expression",
       x = "Expression Value",
       y = "Density")





#==============================================================================|
# GEO extras (volcano etc.) - unchanged
#==============================================================================|

# Direction: Up and Down
table(sign(deg_geo_sig_0.05$logFC))  # -1 = down, +1 = up
summary(deg_geo_sig_0.05)

library(ggplot2)
deg_geo$threshold <- as.factor(deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 0.05)

ggplot(deg_geo, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of GEO DEGs",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")




#======================Block of common genes and GO nd KEGG=====================
# CHANGED: use gene lists created earlier (genes_tcga_1.5, genes_tcga_0.05) which we defined from deg_tcga
common_genes_1.5 <- intersect(gene_geo_1.5, genes_tcga_1.5)
length(gene_geo_1.5); length(genes_tcga_1.5);
common_genes_0.05 <- intersect(gene_geo_0.05, genes_tcga_0.05)
length(gene_geo_0.05); length(genes_tcga_0.05);

length(common_genes_1.5)
length(common_genes_0.05)
#|
#|
#|
entrez_ids_1.5 <- bitr(common_genes_1.5, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids_0.05 <- bitr(common_genes_0.05, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
#|
#|
# GO
ego_1.5 <- enrichGO(gene = entrez_ids_1.5$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                    ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
ego_0.05 <- enrichGO(gene = entrez_ids_0.05$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                     ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
#|
#|
#|
# KEGG
ekegg_1.5 <- enrichKEGG(gene = entrez_ids_1.5$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
ekegg_0.05 <- enrichKEGG(gene = entrez_ids_0.05$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
#|
#|
# Plots----------
#barplot(ego_0.05, showCategory = 15, title = "GO Enrichment")
#dotplot(ekegg_0.05, showCategory = 15, title = "KEGG Pathway Enrichment")
#===============================================================================








#------------------------[ Optional GSEA (TCGA only) ]------------------------
# CHANGED: ensure deg_tcga exists before using it
gene_stats_1.5<- deg_tcga$log2FoldChange   # CHANGED: adapted to deg_tcga frame (note: original code used deg_tcga$logFC — here use log2FoldChange)
names(gene_stats_1.5) <- rownames(deg_tcga)
#|

head(rownames(deg_tcga), 10)

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)  # optional, adds GSEA plotting functions

# Suppose your vector looks like this:
# gene_stats_1.5 <- setNames(deg_tcga$logFC, rownames(deg_tcga))

# 1️⃣ Convert SYMBOL → ENTREZ IDs
gene_symbols <- names(gene_stats_1.5)

# Map to Entrez IDs
symbol2entrez <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# 2️⃣ Merge back to retain stats for mapped genes
gene_stats_entrez <- gene_stats_1.5[symbol2entrez$SYMBOL]
names(gene_stats_entrez) <- symbol2entrez$ENTREZID

# Remove NAs
gene_stats_entrez <- gene_stats_entrez[!is.na(names(gene_stats_entrez))]

# 3️⃣ Run GSEA KEGG
gsea_kegg_1.5 <- gseKEGG(
  geneList = sort(gene_stats_entrez, decreasing = TRUE),
  organism = "hsa",
  pvalueCutoff = 0.2
)

# 4️⃣ View results
head(gsea_kegg_1.5@result)




dotplot(gsea_kegg_1.5, showCategory = 15, title = "KEGG GSEA logfc1.5 (TCGA-PAAD)")










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
expr_common_1.5 <- expr_tcga[rownames(expr_tcga) %in% common_genes_1.5, ]
expr_common_0.05 <- expr_tcga[rownames(expr_tcga) %in% common_genes_0.05, ]

cat("Expression matrix after gene filtering:\n")
print(dim(expr_common_1.5))
print(dim(expr_common_0.05))

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
expr_common_1.5 <- expr_common_1.5[, valid_samples]
expr_common_0.05 <- expr_common_0.05[, valid_samples]
surv_time <- surv_time[valid_samples]
surv_status <- surv_status[valid_samples]


# Transpose expression (samples × genes)
x_1.5 <- t(expr_common_1.5)
x_0.05 <- t(expr_common_0.05)

# Create survival object
y <- Surv(time = surv_time, event = surv_status)

# Check dimensions again
cat("Final dimensions:\n")
print(dim(x_0.05))
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
    x_0.05, y,
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
x_stable <- x_0.05[, stable_gene_names, drop = FALSE]
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
