
raw <- read.delim('paad.txt', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Extract real sample IDs from first row (excluding gene_id column)
sample_ids <- as.character(unlist(raw[1, -1]))  # remove first column

# Data starts from 3rd row onward — remove 1st (sample names) and 2nd row (repeated "normalized_counts")
expr_data <- raw[-1, ]

View(expr_data)

colnames(expr_data)[1]<-'gene_id'

# Remove repeated "gene_id" in expression data rows
expr_data$gene_id <- sub(".*\\|", "", expr_data$gene_id)

# Set gene_id as rownames and drop that column
rownames(expr_data) <- expr_data$gene_id
expr_data <- expr_data[, -1]

View(expr_data)



# Convert expression values to numeric
expr_matrix <- as.data.frame(lapply(expr_data, as.numeric))
expr_matrix <- as.matrix(expr_matrix)
rownames(expr_matrix) <- rownames(expr_data)

View(expr_matrix)


sample_ids <- colnames(expr_matrix)  # Save column names first

# Extract the 14th-15th characters (sample type codes)
sample_types <- substr(sample_ids, 14, 15)

# Double check lengths match
length(sample_types) == ncol(expr_matrix)  # should return TRUE



valid_types <- sample_types %in% c("01", "11")

# Now safely subset the matrix and other variables
expr_matrix <- expr_matrix[, valid_types]
sample_types <- sample_types[valid_types]
labels_tcga <- ifelse(sample_types == "01", 1, 0)

sample_ids <- colnames(expr_matrix)
View(expr_matrix)

View(sample_types)
table(labels)





# Load limma package
library(limma)

# Create design matrix: tumor (1) vs normal (0)
group <- factor(labels, levels = c(0, 1))  # 0 = normal, 1 = tumor
design <- model.matrix(~ group)

# Apply limma
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)

# Extract DEGs
deg_tcga <- topTable(fit, coef = 2, number = Inf, adjust = "fdr")

# Filter significant DEGs
deg_tcga_sig <- subset(deg_tcga, adj.P.Val < 0.05 & abs(logFC) > 1)

# Summary
cat("Total DEGs:", nrow(deg_tcga_sig), "\n")
head(deg_tcga_sig)

View(expr_matrix)


library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db,
                  keys = rownames(deg_tcga_sig),
                  column = "SYMBOL",
                  keytype = "ENTREZID",   # since ID looks like ENTREZ (e.g., 100130426)
                  multiVals = "first")

# 6. Clean and save mapped symbols
symbols <- na.omit(symbols)
genes_tcga <- unname(symbols)
#===============================================================================
#===============================================================================








# Load all libraries,almost ig
library(GEOquery); library(limma); library(DESeq2); library(TCGAbiolinks)
library(SummarizedExperiment); library(clusterProfiler); library(org.Hs.eg.db)
library(enrichplot); library(DOSE)







#------------------------[ 2. GEO]------------------------
gse <- getGEO("GSE62452", GSEMatrix = TRUE)[[1]]
expr_geo <- exprs(gse)
pheno <- pData(gse)


#========
View(expr_geo)
dim(expr_geo)  # Genes × Samples data
head(rownames(expr_geo), 5)
head(colnames(expr_geo), 5)
head(pheno[, 1:5], 3) 
#========


labels_geo <- ifelse(pheno$`tissue:ch1` == "Pancreatic tumor", 1, 0)


#==================
names(pheno)
table(pheno$`tissue:ch1`)
length(labels_geo) == ncol(expr_geo)#verify cheyan vendi
#===================




group_geo <- factor(labels_geo)
design_geo <- model.matrix(~group_geo)
fit_geo <- lmFit(expr_geo, design_geo)
fit_geo <- eBayes(fit_geo)
deg_geo <- topTable(fit_geo, coef = 2, number = Inf, adjust.method = "fdr")
deg_geo_sig <- deg_geo[deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 0.5, ]









#==================Extra views ahn==============================================
cat("Number of significant DEGs:", nrow(deg_geo_sig), "\n")
# Direction: Up and Down
table(sign(deg_geo_sig$logFC))  # -1 = down, +1 = up
summary(deg_geo_sig)
sum(is.na(deg_geo_sig$logFC))
sum(is.na(deg_geo_sig$adj.P.Val))



library(ggplot2)# express change kanan vendi ahn this volvcano plot

deg_geo$threshold <- as.factor(deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 0.5)

ggplot(deg_geo, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of GEO DEGs",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")



library(pheatmap)

#Top 50
top_genes <- rownames(deg_geo_sig)[1:min(50, nrow(deg_geo_sig))]
heat_data <- expr_geo[top_genes, ]

#z-score norm
heat_data_scaled <- t(scale(t(heat_data)))
#annotation 
annotation_col <- data.frame(Group = factor(labels_geo))
rownames(annotation_col) <- colnames(heat_data_scaled)
pheatmap(heat_data_scaled,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         main = "Top 50 DEGs Heatmap - GEO")

#==============================extra end========================================







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


#===============================================================================
#==================================GEO ending===================================
#===============================================================================








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
gene_stats <- deg_tcga$logFC
names(gene_stats) <- rownames(deg_tcga)

head(rownames(deg_tcga), 10)


library(clusterProfiler)

gsea_kegg <- gseKEGG(geneList = sort(gene_stats, decreasing = TRUE),
                     organism = "hsa",
                     pvalueCutoff = 0.2)

dotplot(gsea_kegg, showCategory = 15, title = "KEGG GSEA (TCGA-PAAD)")

#===============================================================================
clinical <- read.delim("D:/MSC/MiniProject/Dataset/firebrowse/clini/All_CDEs.txt", stringsAsFactors = FALSE, sep = "\t")
clinical2 <- read.delim("D:/MSC/MiniProject/Dataset/firebrowse/clini/PAAD.clin.merged.picked.txt", stringsAsFactors = FALSE, sep = "\t")

colnames(clinical)
colnames(clinical2)
View(clinical
head(colnames(deg_tcga), 5)
View(deg_tcga_sig)

head(clinical$patient, 5)


ls()

head(colnames(expr_matrix), 5)
colnames(clinical)
# 1. Extract patient IDs from expression column names
expr_patients <- substr(colnames(expr_matrix), 1, 12)         # Gets "TCGA.2J.AAB1"
expr_patients <- gsub("\\.", "-", expr_patients)              # Converts to "TCGA-2J-AAB1"

# 2. Pull matching clinical data
matched_idx <- match(expr_patients, clinical$bcr_patient_barcode)

# 3. Filter expression + clinical to only matched samples
expr_matrix_matched <- expr_matrix[, !is.na(matched_idx)]
clinical_matched <- clinical[matched_idx[!is.na(matched_idx)], ]

# 4. Confirm match
all.equal(substr(colnames(expr_matrix_matched), 1, 12), gsub("-", ".", clinical_matched$bcr_patient_barcode))











     #===============================================================================
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
