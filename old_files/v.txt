# ========================
# Trimmed & Cleaned Pipeline
# For Biomarker Discovery in Pancreatic Cancer
# GEO: GSE62452 | TCGA: PAAD (STAR counts)
# ========================

#------------------------[ 1. Setup ]------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "DESeq2", "TCGAbiolinks", 
                       "SummarizedExperiment", "clusterProfiler", 
                       "org.Hs.eg.db", "enrichplot", "DOSE"))

# Load all libraries
library(GEOquery); library(limma); library(DESeq2); library(TCGAbiolinks)
library(SummarizedExperiment); library(clusterProfiler); library(org.Hs.eg.db)
library(enrichplot); library(DOSE)

#------------------------[ 2. GEO: GSE62452 ]------------------------
gse <- getGEO("GSE62452", GSEMatrix = TRUE)[[1]]
expr_geo <- exprs(gse)
pheno <- pData(gse)
labels_geo <- ifelse(pheno$`tissue:ch1` == "Pancreatic tumor", 1, 0)

group_geo <- factor(labels_geo)
design_geo <- model.matrix(~group_geo)
fit_geo <- lmFit(expr_geo, design_geo)
fit_geo <- eBayes(fit_geo)
deg_geo <- topTable(fit_geo, coef = 2, number = Inf, adjust.method = "fdr")
deg_geo_sig <- deg_geo[deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 0.5, ]

# Probe -> Gene symbol (GPL6244)
gpl <- getGEO("GPL6244", AnnotGPL = TRUE)
gpl_table <- Table(gpl)
probe2gene <- gpl_table[, c("ID", "Gene symbol")]
colnames(probe2gene) <- c("PROBEID", "SYMBOL")
probe2gene$SYMBOL <- sapply(strsplit(probe2gene$SYMBOL, "///", fixed = TRUE), `[`, 1)
probe2gene <- probe2gene[!is.na(probe2gene$SYMBOL) & probe2gene$SYMBOL != "", ]
deg_geo_sig$SYMBOL <- probe2gene$SYMBOL[match(rownames(deg_geo_sig), probe2gene$PROBEID)]
deg_geo_sig <- deg_geo_sig[!is.na(deg_geo_sig$SYMBOL), ]
gene_geo <- deg_geo_sig$SYMBOL

#------------------------[ 3. TCGA: PAAD ]------------------------
query <- GDCquery(project = "TCGA-PAAD", data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  sample.type = c("Primary Tumor", "Solid Tissue Normal"))
GDCdownload(query)
data_exp <- GDCprepare(query)
expr_tcga <- assay(data_exp)
labels_tcga <- ifelse(colData(data_exp)$shortLetterCode == "TP", 1, 0)

# Filter low-expression genes
keep <- rowSums(expr_tcga > 5) >= 5
expr_tcga_filtered <- expr_tcga[keep, ]

# DESeq2
col_data <- data.frame(condition = factor(labels_tcga), row.names = colnames(expr_tcga_filtered))
dds <- DESeqDataSetFromMatrix(countData = expr_tcga_filtered, colData = col_data, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res_ordered <- res[order(res$padj), ]
deg_tcga_sig <- subset(res_ordered, padj < 0.05 & abs(log2FoldChange) > 0.5)

# ENSEMBL to SYMBOL
clean_ids <- gsub("\\..*", "", rownames(deg_tcga_sig))
tcga_symbols <- mapIds(org.Hs.eg.db, keys = clean_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
tcga_symbols <- na.omit(tcga_symbols)
genes_tcga <- unname(tcga_symbols)

#------------------------[ 4. Overlap & Enrichment ]------------------------
common_genes <- intersect(gene_geo, genes_tcga)
entrez_ids <- bitr(common_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO
ego <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

# KEGG
ekegg <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)

# Plots
barplot(ego, showCategory = 15, title = "GO Enrichment")
dotplot(ekegg, showCategory = 15, title = "KEGG Pathway Enrichment")

#------------------------[ 5. Optional GSEA (TCGA only) ]------------------------
gene_stats <- res_ordered$log2FoldChange
names(gene_stats) <- gsub("\\..*", "", rownames(res_ordered))
ensembl2entrez <- mapIds(org.Hs.eg.db, keys = names(gene_stats), column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
gene_stats <- gene_stats[!is.na(ensembl2entrez)]
names(gene_stats) <- ensembl2entrez[!is.na(ensembl2entrez)]
# Remove duplicate ENTREZ IDs
gene_stats <- gene_stats[!duplicated(names(gene_stats))]


gsea_kegg <- gseKEGG(geneList = gene_stats, organism = "hsa", pvalueCutoff = 0.2)
dotplot(gsea_kegg, showCategory = 15, title = "KEGG GSEA (TCGA-PAAD)")








# 1. Confirm which columns contain survival info
colnames(colData(data_exp))
# 2. Check which TCGA samples match gene expression matrix
dim(expr_tcga_filtered)
# 3. Check if you have enough survival follow-up samples (e.g., > 100)
table(colData(data_exp)$vital_status)







# Create survival time vector
surv_time <- ifelse(is.na(colData(data_exp)$days_to_death),
                    colData(data_exp)$days_to_last_follow_up,
                    colData(data_exp)$days_to_death)

# Status: 1 = dead, 0 = alive
surv_status <- ifelse(colData(data_exp)$vital_status == "Dead", 1, 0)

# Check
summary(surv_time)
table(surv_status)





#xtract exprs of common degs
# Make sure gene names are rownames of expr_tcga_filtered
rownames(expr_tcga_filtered) <- gsub("\\..*", "", rownames(expr_tcga_filtered))  # Clean if needed
# Map gene symbols (common_genes) to ENSEMBL IDs
symbol2ensembl <- mapIds(org.Hs.eg.db,
                         keys = common_genes,
                         column = "ENSEMBL",
                         keytype = "SYMBOL",
                         multiVals = "first")

# Clean NAs
symbol2ensembl <- na.omit(symbol2ensembl)

# Now intersect using ENSEMBL IDs
common_genes_tcga <- intersect(symbol2ensembl, rownames(expr_tcga_filtered))

# Extract and transpose
expr_risk <- t(expr_tcga_filtered[common_genes_tcga, ])


# Extract and transpose
expr_risk <- t(expr_tcga_filtered[common_genes_tcga, ])








install.packages("glmnet")
library(glmnet)
library(survival)
library(survminer)

# Prepare data
x <- as.matrix(expr_risk)
y <- Surv(time = surv_time, event = surv_status)



length(common_genes)  # How many common genes found between GEO and TCGA?
length(common_genes_tcga)  # How many are also in TCGA expression matrix?
dim(expr_tcga_filtered)  # Total genes in TCGA after filtering
dim(expr_risk)
dim(x)
# Run LASSO Cox model
set.seed(123)
cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)





summary(surv_time)
table(surv_time <= 0, useNA = "ifany")




# Step 1: Create a logical index of valid samples
valid_samples <- which(surv_time > 0 & !is.na(surv_time))

# Step 2: Subset survival and expression objects
x_valid <- x[valid_samples, ]
y_valid <- y[valid_samples]


dim(x_valid)  # e.g. [160, 68]
length(y_valid)  # Should match number of rows in x_valid




cvfit <- cv.glmnet(x_valid, y_valid, family = "cox", alpha = 1, nfolds = 10)





# Extract selected genes
coef_min <- coef(cvfit, s = "lambda.min")
selected_genes <- rownames(coef_min)[coef_min[, 1] != 0]
selected_genes <- selected_genes[selected_genes != "(Intercept)"]

selected_genes




plot(cvfit)
coef_min <- coef(cvfit, s = "lambda.min")
selected_genes <- rownames(coef_min)[coef_min[, 1] != 0]
selected_genes <- selected_genes[selected_genes != "(Intercept)"]
print(selected_genes)









# Extract the full coefficient vector at lambda.min
coef_vector <- coef(cvfit, s = "lambda.min")

# Keep only selected genes with non-zero coefficients
selected_genes <- rownames(coef_vector)[coef_vector[, 1] != 0]
selected_genes <- selected_genes[selected_genes != "(Intercept)"]

# Get gene expression data for valid samples and selected genes
expr_selected <- x_valid[, selected_genes]

# Extract matching non-zero coefficients
coefficients <- coef_vector[selected_genes, 1]

# Calculate risk score = sum(expr × coef)
risk_scores <- as.vector(expr_selected %*% coefficients)





#grp patients as high risk nd low risk
# Stratify by median risk score
risk_group <- ifelse(risk_scores > median(risk_scores), "High", "Low")

# Optional: make it a factor
risk_group <- factor(risk_group, levels = c("Low", "High"))

# View distribution
table(risk_group)








# Use only valid survival entries
surv_df <- data.frame(
  time = surv_time[valid_samples],
  status = surv_status[valid_samples],
  group = risk_group
)

# Fit KM model
fit_km <- survfit(Surv(time, status) ~ group, data = surv_df)

# Plot
ggsurvplot(fit_km, data = surv_df, pval = TRUE, risk.table = TRUE,
           palette = c("#00BFC4", "#F8766D"),
           title = "Kaplan-Meier Survival by Risk Group")







#ROC Curve for 1-, 3-, 5-Year Prediction
install.packages("survivalROC")
library(survivalROC)

# Define time points (in days)
time_points <- c(365, 1095, 1825)  # 1, 3, 5 years

# Plot ROC for each
par(mfrow = c(1, 3))  # Set layout

for (t in time_points) {
  roc_result <- survivalROC(Stime = surv_df$time,
                            status = surv_df$status,
                            marker = risk_scores,
                            predict.time = t,
                            method = "KM")
  
  plot(roc_result$FP, roc_result$TP, type = "l",
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = paste("ROC Curve -", t/365, "Year(s)"),
       col = "blue", lwd = 2)
  abline(0, 1, lty = 2)
  legend("bottomright", legend = paste("AUC =", round(roc_result$AUC, 3)))
}
