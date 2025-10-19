# Pancreatic Cancer Bioinformatics Pipeline
# Based on Lang et al., Frontiers in Surgery, 2022
# Adapted for pancreatic cancer using GSE62452 and TCGA-PAAD
# Fixed GDCdownload with explicit gdc-client path

# Step 1: Prerequisites and Setup
# Install required R packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GEOquery", "TCGAbiolinks", "limma", "ggplot2", "glmnet", "survival", "survminer", "pheatmap"))
install.packages(c("dplyr", "tidyr"))

# Load libraries
library(GEOquery)
library(TCGAbiolinks)
library(limma)
library(ggplot2)
library(glmnet)
library(survival)
library(survminer)
library(pheatmap)
library(dplyr)
library(tidyr)
library(SummarizedExperiment)



# Step 2: Data Preparation and Processing
# Downloading GEO dataset (GSE62452)
geo_data <- getGEO("GSE62452", GSEMatrix = TRUE, getGPL = TRUE)
gse62452 <- exprs(geo_data[[1]])  # Extract expression matrix
pheno_data <- pData(geo_data[[1]])  # Extract phenotype data (69 tumor, 61 normal)

# Downloading TCGA-PAAD RNA-seq data
query <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query, method = "client", directory = "D:/MSC/temp/new", files.per.chunk = 10)
tcga_data <- GDCprepare(query, directory = "D:/MSC/temp/new")
class(tcga_data)
# Filter samples with complete clinical information
tcga_clinical <- tcga_data$clinical
tcga_expr <- assay(tcga_data)
tcga_samples <- tcga_data$barcode
tcga_clinical <- tcga_clinical[!is.na(tcga_clinical$vital_status), ]  # Exclude incomplete clinical data
tcga_expr <- tcga_expr[, tcga_samples %in% rownames(tcga_clinical)]

# Step 3: Identification of Differentially Expressed Genes (DEGs)
# Normalize GEO data using limma
gse62452_norm <- normalizeBetweenArrays(gse62452)

# Identify DEGs in GEO dataset
design_geo <- model.matrix(~0 + pheno_data$source_name_ch1)  # Tumor vs Normal
colnames(design_geo) <- c("Normal", "Tumor")
unique(pheno_data$source_name_ch1)


pheno_data$source_name_ch1 <- ifelse(grepl("Non-tumor", pheno_data$source_name_ch1, ignore.case = TRUE), "Normal", "Tumor")
unique(pheno_data$source_name_ch1)  # Should show "Tumor" and "Normal"
design_geo <- model.matrix(~0 + pheno_data$source_name_ch1)
colnames(design_geo) <- c("Normal", "Tumor")


table(pheno_data$source_name_ch1)
ncol(gse62452)  # Should match
nrow(pheno_data)



fit_geo <- lmFit(gse62452_norm, design_geo)
contrast_matrix_geo <- makeContrasts(Tumor-Normal, levels=design_geo)
fit_geo <- contrasts.fit(fit_geo, contrast_matrix_geo)
fit_geo <- eBayes(fit_geo)
deg_geo <- topTable(fit_geo, number=Inf, p.value=0.05, lfc=1.5)







# Identify DEGs in TCGA dataset
design_tcga <- model.matrix(~0 + tcga_data$definition)  # Tumor vs Normal

unique(tcga_data$definition)


keep_samples <- tcga_data$definition %in% c("Primary solid Tumor", "Solid Tissue Normal")
tcga_data <- tcga_data[, keep_samples]  # Subset SummarizedExperiment
tcga_expr <- assay(tcga_data)  # Re-extract expression matrix


tcga_data$definition <- ifelse(tcga_data$definition == "Solid Tissue Normal", "Normal", "Tumor")


design_tcga <- model.matrix(~0 + tcga_data$definition)



colnames(design_tcga) <- c("Normal", "Tumor")
fit_tcga <- lmFit(tcga_expr, design_tcga)
contrast_matrix_tcga <- makeContrasts(Tumor-Normal, levels=design_tcga)
fit_tcga <- contrasts.fit(fit_tcga, contrast_matrix_tcga)
fit_tcga <- eBayes(fit_tcga)
deg_tcga <- topTable(fit_tcga, number=Inf, p.value=0.05, lfc=1.5)

# Find common DEGs
common_degs <- intersect(rownames(deg_geo), rownames(deg_tcga))

# Visualize DEGs (Volcano Plot and Heatmap)
# Volcano Plot for GEO
volcano_data <- deg_geo
volcano_data$Significant <- ifelse(volcano_data$P.Value < 0.05 & abs(volcano_data$logFC) > 1.5, "Significant", "Not Significant")
ggplot(volcano_data, aes(x=logFC, y=-log10(P.Value), color=Significant)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values=c("grey", "red")) +
  ggtitle("Volcano Plot of DEGs (GSE62452)")

# Heatmap for common DEGs
pheatmap(gse62452_norm[common_degs, ], annotation_col=pheno_data[, "source_name_ch1", drop=FALSE])


length(common_degs)
nrow(deg_geo)  # DEGs from GSE62452
nrow(deg_tcga)  # DEGs from TCGA-PAAD



gpl <- getGEO("GPL6480")
probe2gene <- Table(gpl)[, c("ID", "Gene Symbol")]

colnames(Table(gpl))


gpl <- getGEO("GPL6480")
probe2gene <- Table(gpl)[, c("ID", "GENE_SYMBOL")]


deg_geo$gene_symbol <- probe2gene[match(rownames(deg_geo), probe2gene$ID), "GENE_SYMBOL"]
deg_geo <- deg_geo[!is.na(deg_geo$gene_symbol), ]  # Remove unmapped probes

head(deg_geo$gene_symbol)
head(rownames(deg_geo))
head(probe2gene$ID)


deg_geo_unfiltered <- topTable(fit_geo, number=Inf)  # No p-value or logFC filter
nrow(deg_geo_unfiltered)
dim(gse62452_norm)
dim(design_geo)




deg_geo <- topTable(fit_geo, number=Inf, p.value=0.1, lfc=1)  # Looser: p < 0.1, |logFC| > 1
nrow(deg_geo)  # Check how many DEGs now
head(rownames(deg_geo))  # Should show probe IDs like "A_23_P..."



gpl <- getGEO("GPL6480")
probe2gene <- Table(gpl)[, c("ID", "GENE_SYMBOL")]
deg_geo$gene_symbol <- probe2gene[match(rownames(deg_geo), probe2gene$ID), "GENE_SYMBOL"]
deg_geo <- deg_geo[!is.na(deg_geo$gene_symbol), ]  # Remove unmapped probes
head(deg_geo$gene_symbol)  # Should show gene symbols now


deg_geo <- topTable(fit_geo, number=Inf, p.value=0.1, lfc=1)
nrow(deg_geo)
head(rownames(deg_geo))
sum(is.na(probe2gene$GENE_SYMBOL))
sum(is.na(deg_geo$gene_symbol))
head(rownames(deg_tcga))

-------------------
  
  
  
gpl <- getGEO("GPL6480")
probe2gene <- Table(gpl)[, c("ID", "GENE_SYMBOL")]
deg_geo$gene_symbol <- probe2gene[match(rownames(deg_geo), probe2gene$ID), "GENE_SYMBOL"]
head(deg_geo$gene_symbol)  # Should now show gene symbols

------------------------------
  
head(rownames(deg_geo), 10)
head(probe2gene$ID, 10)




head(rownames(deg_geo), 10)
sum(rownames(deg_geo) %in% probe2gene$ID)


geo_data <- getGEO("GSE62452", GSEMatrix = TRUE, getGPL = TRUE)
gse62452 <- exprs(geo_data[[1]])
fdata <- fData(geo_data[[1]])

head(rownames(gse62452))
head(fdata[, 1:5])

-------------------------------
gpl <- getGEO("GPL6480", AnnotGPL = TRUE)
gpl_table <- Table(gpl)

# Check available columns
colnames(gpl_table)




# Check possible numeric match
head(gpl_table$SPOT_ID)
head(gpl_table$ID)
head(gpl_table$GENE_SYMBOL)

# Attempt matching rownames(gse62452) to SPOT_ID
probe2gene <- gpl_table[, c("SPOT_ID", "GENE_SYMBOL")]
colnames(probe2gene) <- c("ID", "GENE_SYMBOL")  # rename for consistency

# Map gene symbols
deg_geo$gene_symbol <- probe2gene[match(rownames(deg_geo), probe2gene$ID), "GENE_SYMBOL"]
head(deg_geo$gene_symbol)

-----------------------------------

geo_data <- getGEO("GSE62452", GSEMatrix = TRUE, getGPL = TRUE)
gse62452 <- geo_data[[1]]
expr_matrix <- exprs(gse62452)
fdata <- fData(gse62452)

head(rownames(expr_matrix))         # should be like "A_23_P..."
head(fdata[, c("ID", "GENE_SYMBOL")], 10)


colnames(fdata)


# Extract gene symbols from gene_assignment
fdata$gene_symbol <- sapply(strsplit(as.character(fdata$gene_assignment), " /// "), function(x) {
  first_entry <- x[1]
  gene <- unlist(strsplit(first_entry, " "))[2]
  return(gene)
})

head(fdata$gene_symbol)
head(fdata$gene_assignment, 10)



gpl <- getGEO("GPL6480")
gpl_table <- Table(gpl)
colnames(gpl_table)

#-----------------
  # Force reload and get correct GSE62452 data with probe IDs
# ⚠️ Set tempdir to avoid Windows issues (optional but helpful)
options(timeout = 300)
# Force R to use a new, simpler temp directory
tempdir <- "D:/MSC/temp"
dir.create(tempdir, showWarnings = FALSE)
Sys.setenv(TMPDIR = tempdir)
options(timeout = 600)  # increase timeout in case of slow internet

# 🔁 Re-download full GEO dataset (GSEMatrix = TRUE + getGPL = TRUE)
geo_data <- getGEO("GSE62452", GSEMatrix = TRUE, getGPL = TRUE)

# Extract expression and pheno
gse62452 <- geo_data[[1]]
expr_data <- exprs(gse62452)
pheno_data <- pData(gse62452)


# Manually parse the series matrix
gse_file <- GEOquery::getGEO("GSE62452", GSEMatrix = TRUE, getGPL = TRUE)[[1]]
expr_data <- exprs(gse_file)
pheno_data <- pData(gse_file)
head(rownames(expr_data), 5)









# Step 4: Functional Enrichment Analysis
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)

# Convert gene symbols to Entrez IDs
gene_ids <- bitr(common_degs, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# GO Enrichment
go_enrich <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05)
# KEGG Enrichment
kegg_enrich <- enrichKEGG(gene = gene_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)

# Visualize enrichment results
dotplot(go_enrich, showCategory=10)
dotplot(kegg_enrich, showCategory=10)

# Step 5: Cox and Lasso Regression Analysis
# Prepare survival data
surv_data <- tcga_clinical[, c("vital_status", "days_to_last_follow_up", "days_to_death")]
surv_data$time <- ifelse(surv_data$vital_status == "Dead", surv_data$days_to_death, surv_data$days_to_last_follow_up)
surv_data$status <- ifelse(surv_data$vital_status == "Dead", 1, 0)

# Univariate Cox regression
cox_results <- lapply(common_degs, function(gene) {
  cox_model <- coxph(Surv(time, status) ~ tcga_expr[gene, ], data=surv_data)
  summary(cox_model)$coefficients
})
cox_significant <- names(cox_results)[sapply(cox_results, function(x) x[5] < 0.05)]

# Lasso regression
x <- t(tcga_expr[cox_significant, ])
y <- Surv(surv_data$time, surv_data$status)
lasso_model <- glmnet(x, y, family="cox")
lasso_coef <- coef(lasso_model, s=lambda.min)
selected_genes <- rownames(lasso_coef)[lasso_coef != 0]

# Multivariate Cox regression
final_model <- coxph(Surv(time, status) ~ ., data=data.frame(t(tcga_expr[selected_genes, ])))
risk_scores <- predict(final_model, type="risk")
median_risk <- median(risk_scores)
risk_groups <- ifelse(risk_scores > median_risk, "High", "Low")

# Step 6: Survival Analysis
# Kaplan-Meier survival curves
surv_object <- Surv(time=surv_data$time, event=surv_data$status)
fit_surv <- survfit(surv_object ~ risk_groups)
ggsurvplot(fit_surv, data=surv_data, pval=TRUE, risk.table=TRUE)

# ROC curves
library(survivalROC)
roc_1yr <- survivalROC(Stime=surv_data$time, status=surv_data$status, marker=risk_scores, predict.time=365, method="KM")
roc_3yr <- survivalROC(Stime=surv_data$time, status=surv_data$status, marker=risk_scores, predict.time=365*3, method="KM")
roc_5yr <- survivalROC(Stime=surv_data$time, status=surv_data$status, marker=risk_scores, predict.time=365*5, method="KM")
plot(roc_1yr$FP, roc_1yr$TP, type="l", col="blue", xlab="False Positive Rate", ylab="True Positive Rate", main="ROC Curves")
lines(roc_3yr$FP, roc_3yr$TP, col="red")
lines(roc_5yr$FP, roc_5yr$TP, col="green")
legend("bottomright", c("1-year", "3-year", "5-year"), col=c("blue", "red", "green"), lty=1)

# Step 7: Immune Infiltration Analysis (CIBERSORT)
# Note: Requires CIBERSORT R script and LM22 signature matrix from CIBERSORT website
# For simplicity, assume CIBERSORT results are available as 'cibersort_results'
# wilcox.test(cibersort_results$cell_type ~ risk_groups)

# Step 8: Gene Set Enrichment Analysis (GSEA)
BiocManager::install("fgsea")
library(fgsea)
# Prepare ranked gene list based on risk scores
gene_rank <- rankGenes(tcga_expr[common_degs, ])
pathways <- gmtPathways("c2.cp.kegg.v7.4.symbols.gmt")  # Download from MSigDB
fgsea_results <- fgsea(pathways, gene_rank, nperm=1000)
plotEnrichment(pathways[[1]], gene_rank)

# Save results
write.csv(common_degs, "common_degs_pancreatic.csv")
write.csv(risk_scores, "risk_scores_pancreatic.csv")