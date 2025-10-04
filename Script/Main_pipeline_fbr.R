#------------------------[ 1. Download TCGA-PAAD RNA-seq Data with TCGAbiolinks ]------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)

# Query for gene expression (STAR counts) for both tumor and normal
query <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(query, method = "client", files.per.chunk = 10)
data <- GDCprepare(query)

#------------------------[ 2. Prepare Expression Matrix and Labels ]------------------------
expr_matrix <- assay(data)  # genes x samples
pheno <- colData(data)

# Extract sample types from barcode
sample_types <- substr(pheno$barcode, 14, 15)
labels <- ifelse(sample_types == "01", "Tumor", "Normal")

# Filter for tumor and normal samples
valid_types <- sample_types %in% c("01", "11")
expr_matrix <- expr_matrix[, valid_types]
labels <- labels[valid_types]
pheno <- pheno[valid_types, ]

#------------------------[ 3. Create DESeq2 Dataset ]------------------------
dds <- DESeqDataSetFromMatrix(
  countData = expr_matrix,
  colData = data.frame(condition = factor(labels, levels = c("Normal", "Tumor"))),
  design = ~ condition
)

#------------------------[ 4. Run DESeq2 Analysis ]------------------------
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
res <- lfcShrink(dds, coef = "condition_Tumor_vs_Normal", res = res, type = "ashr") # optional, for more accurate logFC

#------------------------[ 5. Extract and Filter DEGs ]------------------------
# Order by adjusted p-value
resOrdered <- res[order(res$padj), ]

# Filter for significant DEGs (adjust as needed)
deg_deseq2_1 <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)
deg_deseq2_0.5 <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 0.5)
deg_deseq2_0.05 <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 0.05)

cat("Number of significant DEGs (logFC>1):", nrow(deg_deseq2_1), "\n")
cat("Number of significant DEGs (logFC>0.5):", nrow(deg_deseq2_0.5), "\n")
cat("Number of significant DEGs (logFC>0.05):", nrow(deg_deseq2_0.05), "\n")

#==============================================================================================#|




library(org.Hs.eg.db)
symbols_1.5 <- mapIds(org.Hs.eg.db,
                      keys = rownames(deg_tcga_sig_1.5),
                      column = "SYMBOL",
                      keytype = "ENTREZID",   # since ID looks like ENTREZ (e.g., 100130426)
                      multiVals = "first")


symbols_0.05 <- mapIds(org.Hs.eg.db,
                       keys = rownames(deg_tcga_sig_0.05),
                       column = "SYMBOL",
                       keytype = "ENTREZID",   # since ID looks like ENTREZ (e.g., 100130426)
                       multiVals = "first")


# 6. Clean and save mapped symbols
symbols_1.5 <- na.omit(symbols_1.5)
symbols_0.05 <- na.omit(symbols_0.05)
length(symbols_0.05)
genes_tcga_1.5 <- unname(symbols_1.5)
genes_tcga_0.05<- unname(symbols_0.05)
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
expr_matrix <- as.matrix(expr_geo_annot[, -(1:2)])  # Remove PROBEID and SYMBOL columns
rowGroup <- expr_geo_annot$SYMBOL
rowID <- expr_geo_annot$PROBEID

collapsed <- collapseRows(datET = expr_matrix,
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

#------------------------[ 6. Filter Significant DEGs ]------------------------
deg_geo_sig_1.5 <- deg_geo[deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 1, ]
deg_geo_sig_0.5 <- deg_geo[deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 0.5, ]
deg_geo_sig_0.05 <- deg_geo[deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 0.05, ]

cat("Number of significant DEGs 1.5:", nrow(deg_geo_sig_1.5), "\n")
cat("Number of significant DEGs 0.5:", nrow(deg_geo_sig_0.5), "\n")
cat("Number of significant DEGs 0.05:", nrow(deg_geo_sig_0.05), "\n")

#------------------------[ 7. Map DEGs to Gene Symbols ]------------------------
deg_geo_sig_1.5$SYMBOL <- probe2gene$SYMBOL[match(rownames(deg_geo_sig_1.5), probe2gene$PROBEID)]
deg_geo_sig_1.5 <- deg_geo_sig_1.5[!is.na(deg_geo_sig_1.5$SYMBOL), ]
gene_geo_1.5 <- deg_geo_sig_1.5$SYMBOL

deg_geo_sig_0.05$SYMBOL <- probe2gene$SYMBOL[match(rownames(deg_geo_sig_0.05), probe2gene$PROBEID)]
deg_geo_sig_0.05 <- deg_geo_sig_0.05[!is.na(deg_geo_sig_0.05$SYMBOL), ]
gene_geo_0.05 <- deg_geo_sig_0.05$SYMBOL

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
#|
#|
#|
#==================Extra views ahn==============================================

# Direction: Up and Down
table(sign(deg_geo_sig_0.05$logFC))  # -1 = down, +1 = up
summary(deg_geo_sig_0.05)


library(ggplot2)# express change kanan vendi ahn this volvcano plot

deg_geo$threshold <- as.factor(deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 0.05)

ggplot(deg_geo, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of GEO DEGs",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")



library(pheatmap)



#==============================extra end========================================




#===============================================================================
#==================================GEO ending===================================
#===============================================================================







#======================Block of common genes and GO nd KEGG=====================
common_genes_1.5 <- intersect(gene_geo_1.5, genes_tcga_1.5)
length(gene_geo_1.5);length(genes_tcga_1.5);
common_genes_0.05 <- intersect(gene_geo_0.05, genes_tcga_0.05)
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
#barplot(ego_1.5, showCategory = 15, title = "GO Enrichment")!!!!!!!!
#dotplot(ekegg_1.5, showCategory = 15, title = "KEGG Pathway Enrichment")!!!!!!!

barplot(ego_0.05, showCategory = 15, title = "GO Enrichment")
dotplot(ekegg_0.05, showCategory = 15, title = "KEGG Pathway Enrichment")
#===============================================================================









#------------------------[ Optional GSEA (TCGA only) ]------------------------
gene_stats_1.5<- deg_tcga_sig_0.05$logFC
names(gene_stats_1.5) <- rownames(deg_tcga_sig_0.05)
#|

head(rownames(deg_tcga_sig_0.05), 10)


library(clusterProfiler)

gsea_kegg_1.5 <- gseKEGG(geneList = sort(gene_stats_1.5, decreasing = TRUE),
                         organism = "hsa",
                         pvalueCutoff = 0.2)



dotplot(gsea_kegg_1.5, showCategory = 15, title = "KEGG GSEA logfc1.5 (TCGA-PAAD)")


#===============================================================================



#------------------------[ 1. Download Clinical Data from TCGA ]------------------------
library(TCGAbiolinks)

# Download clinical data for TCGA-PAAD
clin_df <- GDCquery_clinic(project = "TCGA-PAAD", type = "clinical")
head(clin_df)

#------------------------[ 2. Match Clinical Data to Expression Data ]------------------------
# Extract patient barcodes (first 12 characters) from expression data
expr_ids <- substr(gsub("\\.", "", colnames(expr_matrix)), 1, 12)

# Match clinical data to expression data
matched_ids <- intersect(expr_ids, clin_df$submitter_id)
cat("✅ Matched samples: ", length(matched_ids), "\n")

# Subset and reorder expression and clinical data
expr_matrix <- expr_matrix[, expr_ids %in% matched_ids]
clin_df_matched <- clin_df[match(expr_ids[expr_ids %in% matched_ids], clin_df$submitter_id), ]

#------------------------[ 3. Extract Survival Information ]------------------------
# Convert survival columns to numeric and vital status to 0/1
clin_df_matched$vital_status <- ifelse(clin_df_matched$vital_status == "Alive", 0, 1)
clin_df_matched$days_to_death <- as.numeric(clin_df_matched$days_to_death)
clin_df_matched$days_to_last_follow_up <- as.numeric(clin_df_matched$days_to_last_follow_up)

# Calculate survival time and status
surv_time <- ifelse(is.na(clin_df_matched$days_to_death),
                    clin_df_matched$days_to_last_follow_up,
                    clin_df_matched$days_to_death)
surv_status <- clin_df_matched$vital_status

# Final check
summary(surv_time)
table(surv_status)
dim(expr_matrix)








#===============================================================================
#LASSO + Cox Regression Modeling 
#===============================================================================

length(common_genes_1.5)
length(common_genes_0.05)

head(rownames(expr_matrix), 10)

head(common_genes_1.5)
head(common_genes_0.05)








library(org.Hs.eg.db)

# Entrez 
symbols <- mapIds(org.Hs.eg.db,
                  keys = rownames(expr_matrix),
                  column = "SYMBOL",
                  keytype = "ENTREZID",
                  multiVals = "first")

# Remove rows with NA JUST FR SAFE
valid_idx <- !is.na(symbols)
expr_matrix <- expr_matrix[valid_idx, ]
symbols <- symbols[valid_idx]

# Assign gene symbols as rownames
rownames(expr_matrix) <- symbols

# Remove duplicated gene symbols (if any)
expr_matrix <- expr_matrix[!duplicated(rownames(expr_matrix)), ]

# Confirm new gene symbols
cat("New expression matrix rownames (symbols):\n")
print(head(rownames(expr_matrix)))

#///////////////////////////////////////////////////////////////////////////////
# Filter expression matrix to common genes======================================
expr_common_1.5 <- expr_matrix[rownames(expr_matrix) %in% common_genes_1.5, ]
expr_common_0.05 <- expr_matrix[rownames(expr_matrix) %in% common_genes_0.05, ]

# Double-check dimensions
cat("Expression matrix after gene filtering:\n")
dim(expr_common_1.5)
dim(expr_common_0.05)
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
gene_stability <- as.numeric(gene_freq) / n_runs
names(gene_stability) <- names(gene_freq)

# View stable genes (≥ 70% of runs)
stable_genes <- gene_stability[gene_stability >= 0.7]
print(stable_genes)



# Plot stability
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

gene_stability <- gene_freq / n_runs
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




#===============================================================================
#|                                                                            #|
#|                                                                            #|
#|                                                                            #|
#|                                                                            #|
#|                                                                            #|
#|                                                                            #|
#|                                                                            #|







