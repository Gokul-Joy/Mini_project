
raw <- read.delim('paad.txt', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Extract real sample IDs from first row (excluding gene_id column)
sample_ids <- as.character(unlist(raw[1, -1]))  # remove first column

# Data starts from 3rd row onward — remove 1st (sample names) and 2nd row (repeated "normalized_counts")
expr_data <- raw[-1, ]



colnames(expr_data)[1]<-'gene_id'

# Remove repeated "gene_id" in expression data rows
expr_data$gene_id <- sub(".*\\|", "", expr_data$gene_id)

# Set gene_id as rownames and drop that column
rownames(expr_data) <- expr_data$gene_id
expr_data <- expr_data[, -1]





# Convert expression values to numeric
expr_matrix <- as.data.frame(lapply(expr_data, as.numeric))
expr_matrix <- as.matrix(expr_matrix)
rownames(expr_matrix) <- rownames(expr_data)

dim(expr_matrix)


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
labels<-labels_tcga
sample_ids <- colnames(expr_matrix)






# Load limma package
library(limma)

# Create design matrix: tumor (1) vs normal (0)
group <- factor(labels, levels = c(0, 1))  # 0 = normal, 1 = tumor
design <- model.matrix(~ group)


#log transform================================================================#|
log_mat <- log2(expr_matrix + 1)                                              #|
#=============================================================================#|




# Apply limma=================
fit <- lmFit(log_mat, design)
fit <- eBayes(fit,trend = TRUE)
#=============================

plotSA(fit, main="Mean-variance with log transform")

# Extract DEGs
deg_tcga <- topTable(fit, coef = 2, number = Inf, adjust = "fdr")





#===============================================================================
#To see the histogram of pval desgs

# Filter significant DEGs
sig_degs <- deg_tcga[deg_tcga$adj.P.Val < 0.05, ]


# Plot histogram of logFC values
hist(
  sig_degs$logFC,
  breaks = 50,                    # Number of bins
  col = "skyblue",                # Color of bars
  border = "white",               # No border color
  main = "Distribution of log2 Fold Change (Significant DEGs)",
  xlab = "log2 Fold Change (logFC)",
  ylab = "Number of Genes"
)

# Add vertical lines for thresholds (optional)
abline(v = c(-1.5, -1, -0.5, 0.5, 1, 1.5), col = "red", lty = 2)

summary(sig_degs$logFC)
range(sig_degs$logFC, na.rm = TRUE)
#===============================================================================







# Filter significant DEGs===========================================================|
deg_tcga_sig_1.5 <- subset(deg_tcga, adj.P.Val < 0.05 & abs(logFC) > 1.5)          #|
deg_tcga_sig_1 <- subset(deg_tcga, adj.P.Val < 0.05 & abs(logFC) > 1)              #|
deg_tcga_sig_0.5 <- subset(deg_tcga, adj.P.Val < 0.05 & abs(logFC) > 0.5)          #|
deg_tcga_sig_0.05 <- subset(deg_tcga, adj.P.Val < 0.05 & abs(logFC) > 0.05)        #|   
deg_tcga_sig_0.05_pval_0.1 <- subset(deg_tcga, adj.P.Val < 0.1 & abs(logFC) > 0.05)#|
#===================================================================================|



# Summary=======================================================================================|
cat("Total DEGs for 1.5:", nrow(deg_tcga_sig_1.5), "\n")                                       #|
cat("Total DEGs for 1:", nrow(deg_tcga_sig_1), "\n")                                           #|
cat("Total DEGs for 0.5:", nrow(deg_tcga_sig_0.5), "\n")                                       #|
cat("Total DEGs for fold >0.05:", nrow(deg_tcga_sig_0.05), "\n")                                #|
cat("Total DEGs for logFC > 0.05 and adj.P.Val < 0.1:", nrow(deg_tcga_sig_0.05_pval_0.1), "\n")#|
#==============================================================================================#|

View(deg_tcga_sig)
View(deg_tcga_sig_1.5)


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







# Load all libraries,almost ig
library(GEOquery); library(limma); library(DESeq2); library(TCGAbiolinks)
library(SummarizedExperiment); library(clusterProfiler); library(org.Hs.eg.db)
library(enrichplot); library(DOSE)







#------------------------[ 2. GEO]------------------------

#====================================================================================#|
gse <- getGEO(filename = 'D:/MSC/MiniProject/Dataset/GSE62452_series_matrix.txt.gz') #|
#====================================================================================#|

expr_geo <- exprs(gse)
pheno <- pData(gse)


#========

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
dim(deg_geo)
dim(deg_tcga)
#===================================================================================|
deg_geo_sig_1.5 <- deg_geo[deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 1.5, ]  #|
deg_geo_sig_0.5 <- deg_geo[deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 0.5, ]  #|
deg_geo_sig_0.05 <- deg_geo[deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 0.05, ]#|
#==================================================================================#|




#==============================================================================|
cat("Number of significant DEGs 1.5:", nrow(deg_geo_sig_1.5), "\n")
cat("Number of significant DEGs 0.5:", nrow(deg_geo_sig_0.5), "\n")
cat("Number of significant DEGs 0.05:", nrow(deg_geo_sig_0.05), "\n")

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

#Top 50
top_genes <- rownames(deg_geo_sig_0.05)[1:min(50, nrow(deg_geo_sig_0.05))]
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
#|
#|
#|
#|
#Trying two set of LOGFC
#--------------------
deg_geo_sig_1.5$SYMBOL <- probe2gene$SYMBOL[match(rownames(deg_geo_sig_1.5), probe2gene$PROBEID)]
deg_geo_sig_1.5 <- deg_geo_sig_1.5[!is.na(deg_geo_sig_1.5$SYMBOL), ]
gene_geo_1.5 <- deg_geo_sig_1.5$SYMBOL
#-------------------
deg_geo_sig_0.05$SYMBOL <- probe2gene$SYMBOL[match(rownames(deg_geo_sig_0.05), probe2gene$PROBEID)]
deg_geo_sig_0.05 <- deg_geo_sig_0.05[!is.na(deg_geo_sig_0.05$SYMBOL), ]
gene_geo_0.05 <- deg_geo_sig_0.05$SYMBOL

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
gene_stats_1.5<- deg_tcga$logFC
names(gene_stats_1.5) <- rownames(deg_tcga)
#|

head(rownames(deg_tcga), 10)


library(clusterProfiler)

gsea_kegg_1.5 <- gseKEGG(geneList = sort(gene_stats_1.5, decreasing = TRUE),
                     organism = "hsa",
                     pvalueCutoff = 0.2)



dotplot(gsea_kegg_1.5, showCategory = 15, title = "KEGG GSEA logfc1.5 (TCGA-PAAD)")


#===============================================================================







#  column names
cat("Expression matrix sample IDs:\n")
print(head(colnames(expr_matrix)))

#  clinical file again
clin_raw <- read.delim("D:/MSC/MiniProject/Dataset/firebrowse/clini/PAAD.clin.merged.picked.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# =============================


# Extract base sample ids TCGA2JAAB1
expr_ids <- gsub("\\.", "", substr(colnames(expr_matrix), 1, 12))  # remove dots
cat("Cleaned expr IDs:\n")
print(head(expr_ids))

# =============================
# =============================


clin_data <- clin_raw[-1, ]
var_names <- clin_data[[1]]

sample_ids_clin <- toupper(gsub("-", "", gsub("TCGA-", "", colnames(clin_data)[-1])))  # already like "TCGA2JAABR"
clin_matrix <- as.matrix(clin_data[, -1])
rownames(clin_matrix) <- var_names
colnames(clin_matrix) <- sample_ids_clin

clin_df <- as.data.frame(t(clin_matrix), stringsAsFactors = FALSE)
clin_df$SampleID <- rownames(clin_df)

cat("Cleaned clinical SampleIDs:\n")
print(head(clin_df$SampleID))

# =============================
# 3. Match samples
# =============================
matched_ids <- intersect(expr_ids, clin_df$SampleID)
cat("✅ Matched samples: ", length(matched_ids), "\n")

# Subset and reorder
expr_matrix <- expr_matrix[, expr_ids %in% matched_ids]
clin_df <- clin_df[clin_df$SampleID %in% matched_ids, ]
clin_df <- clin_df[match(gsub("\\.", "", substr(colnames(expr_matrix), 1, 12)), clin_df$SampleID), ]

# =============================
# 4. Coerce and create survival objects
# =============================
clin_df$vital_status <- as.numeric(clin_df$vital_status)
clin_df$days_to_death <- as.numeric(clin_df$days_to_death)
clin_df$days_to_last_followup <- as.numeric(clin_df$days_to_last_followup)

surv_time <- ifelse(is.na(clin_df$days_to_death), clin_df$days_to_last_followup, clin_df$days_to_death)
surv_status <- clin_df$vital_status

# Final check
summary(surv_time)
table(surv_status)
dim(expr_matrix)







#===============================================================================
#LASSO + Cox Regression Modeling (on TCGA expression)
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

#==============
#fit_1.5 <- cv.glmnet(x_1.5, y, family = "cox", alpha = 1, nfolds = 10)!!!!!!!!!
fit_0.05 <- cv.glmnet(x_0.05, y, family = "cox", alpha = 1, nfolds = 10)
#==============
#|
#==============
#plot(fit_1.5)!!!!!!
plot(fit_0.05)
#==============
# Extract non-zero genes==================================
#lasso_coef_1.5 <- coef(fit_1.5, s = fit_1.5$lambda.min)!!!!!!!!
lasso_coef_0.05 <- coef(fit_0.05, s = fit_0.05$lambda.min)
#|
#lasso_genes_1.5 <- rownames(lasso_coef_1.5)[lasso_coef_1.5[, 1] != 0]!!!!!!!!!!
lasso_genes_0.05 <- rownames(lasso_coef_0.05)[lasso_coef_0.05[, 1] != 0]
#|
#|
cat("✅ Selected   biomarker  for logfc 1.5 and 0.05:\n")
#print(lasso_genes_1.5)!!!
#length(lasso_genes_1.5)!!!
print(lasso_genes_0.05)
length(lasso_genes_0.05)
#=======================

dim(x_1.5)
dim(x_0.05)
#==========
#|
#|
#|
# Genes from lambda.min=========================================================
#genes_min_1.5 <- rownames(coef(fit_1.5, s = "lambda.min"))[coef(fit_1.5, s = "lambda.min")[,1] != 0]
genes_min_0.05 <- rownames(coef(fit_0.05, s = "lambda.min"))[coef(fit_0.05, s = "lambda.min")[,1] != 0]
#|
#|
# Genes from lambda.1se
#genes_1se_1.5 <- rownames(coef(fit_1.5, s = "lambda.1se"))[coef(fit_1.5, s = "lambda.1se")[,1] != 0]
genes_1se_0.05 <- rownames(coef(fit_0.05, s = "lambda.1se"))[coef(fit_0.05, s = "lambda.1se")[,1] != 0]
#|
#|
# Print comparison
#cat("Genes (lambda.min) for logfc 1.5:", length(genes_min_1.5), "\n", genes_min_1.5, "\n\n")
cat("Genes (lambda.min) for logfc 0.05:", length(genes_min_0.05), "\n", genes_min_0.05, "\n\n")
#cat("Genes (lambda.1se) for logcf 1.5:", length(genes_1se_1.5), "\n", genes_1se_1.5)
cat("Genes (lambda.1se) for logcf 0.05:", length(genes_1se_0.05), "\n", genes_1se_0.05)
#===============================================================================









#===============================================================================
#Risk Score Calculation + Grouping
#===============================================================================

# Get coefficients (as named vector)
#coef_vector_1.5 <- as.vector(lasso_coef_1.5[lasso_coef_1.5[, 1] != 0])!!!!
coef_vector_0.05 <- as.vector(lasso_coef_0.05[lasso_coef_0.05[, 1] != 0])
#names(coef_vector_1.5) <- lasso_genes_1.5!!!!!!!
names(coef_vector_0.05) <- lasso_genes_0.05
#|
# Subset only selected genes from x
#x_lasso_1.5 <- x_1.5[, lasso_genes_1.5]!!
x_lasso_0.05<- x_0.05[, lasso_genes_0.05]
#|
#|
#|
#|
#Risk score per patient (sample)
#risk_score_1.5 <- as.numeric(x_lasso_1.5 %*% coef_vector_1.5)!!!!
risk_score_0.05<- as.numeric(x_lasso_0.05 %*% coef_vector_0.05)
#|
# Median split into high- and low-risk
#median_cutoff_1.5 <- median(risk_score_1.5)!!!!!!!!!!!!
median_cutoff_0.05 <- median(risk_score_0.05)
#|
#risk_group_1.5 <- ifelse(risk_score_1.5 > median_cutoff_1.5, "High", "Low")!!!!!!!!!!!!
risk_group_0.05 <- ifelse(risk_score_0.05> median_cutoff_0.05, "High", "Low")
#|
# Convert to factor
#risk_group_1.5<- factor(risk_group_1.5, levels = c("Low", "High"))!!!!
risk_group_0.05 <- factor(risk_group_0.05, levels = c("Low", "High"))
#===============================================================================






#KAplan meirer plot=============================================================
library(survival)
library(survminer)

# Survival object
surv_obj <- Surv(surv_time, surv_status)

# Fit KM
#fit_km_1.5 <- survfit(surv_obj ~ risk_group_1.5)!!!!!!!
fit_km_0.05 <- survfit(surv_obj ~ risk_group_0.05)
#|
#|
# Plot
#ggsurvplot(fit_km_1.5,
 #          data = data.frame(risk_group_1.5),
 #          pval = TRUE,
  #         risk.table = TRUE,
  #         title = "Survival Curve: High vs Low Risk logfc 1.5",
   #        palette = c("blue", "red"))
#
#|
#|
#|
#|

ggsurvplot(fit_km_0.05,
           data = data.frame(risk_group_0.05),
           pval = TRUE,
           risk.table = TRUE,
           title = "Survival Curve: High vs Low Risk logfc 0.05",
           palette = c("blue", "red"))



#ROC Curve (1/3/5-year Survival)================================================


library(timeROC)

roc_obj <- timeROC(T = surv_time,
                   delta = surv_status,
                   marker = risk_score_0.05,
                   cause = 1,
                   times = c(365, 1095, 1825),
                   iid = TRUE)

plot(roc_obj, time = 365, col = "blue", title = TRUE)
plot(roc_obj, time = 1095, add = TRUE, col = "green")
plot(roc_obj, time = 1825, add = TRUE, col = "red")
legend("bottomright", legend = c("1-year", "3-year", "5-year"),
       col = c("blue", "green", "red"), lwd = 2)



#===============================================================================
#|                                                                            #|
#|                                                                            #|
#|                                                                            #|
#|                                                                            #|
#|                                                                            #|
#|                                                                            #|
#|                                                                            #|
#|                                                                            #|
#=================================Change of regression==========================

gene_pvals_1.5 <- apply(x_1.5, 2, function(g) {
  summary(coxph(Surv(surv_time[valid_samples], surv_status[valid_samples]) ~ g))$coefficients[,"Pr(>|z|)"]
})
gene_pvals_0.05 <- apply(x_0.05, 2, function(g) {
  summary(coxph(Surv(surv_time[valid_samples], surv_status[valid_samples]) ~ g))$coefficients[,"Pr(>|z|)"]
})



sort(gene_pvals_1.5)[1:10]
sort(gene_pvals_0.05)[1:10]
length(gene_pvals_0.05)
#===============================================================================













#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Extract genes & check stability across runs====================================
library(glmnet)

# Fit your model
fit <- cv.glmnet(x_0.05, y, family = "cox", alpha = 1, nfolds = 10)

# Extract coefficients at optimal lambda
coef_opt <- coef(fit, s = "lambda.min")   # s = "lambda.1se" is more conservative
selected_genes <- rownames(coef_opt)[which(coef_opt != 0)]

selected_genes



n_runs <- 50
gene_list <- list()

for (i in 1:n_runs) {
  fit_tmp <- cv.glmnet(x_0.05, y, family = "cox", alpha = 1, nfolds = 10)
  coef_tmp <- coef(fit_tmp, s = "lambda.min")
  genes_tmp <- rownames(coef_tmp)[which(coef_tmp != 0)]
  gene_list[[i]] <- genes_tmp
}

# Frequency table: how many times each gene was selected
gene_freq <- sort(table(unlist(gene_list)), decreasing = TRUE)

gene_freq






#===============================================================================
library(glmnet)
library(tibble)
library(dplyr)
library(ggplot2)

# Parameters
n_runs <- 100     # Number of repeated LASSO CV runs
alpha_val <- 1    # LASSO penalty
nfolds_val <- 10  # CV folds
seed_base <- 123  # Base seed for reproducibility

gene_list <- list()

for (i in 1:n_runs) {
  set.seed(seed_base + i)  # Different seed each run for fold variability
  
  fit <- cv.glmnet(
    x_0.05, y,
    family = "cox",
    alpha = alpha_val,
    nfolds = nfolds_val
  )
  
  coef_opt <- coef(fit, s = "lambda.min")
  selected_genes <- rownames(coef_opt)[which(coef_opt != 0)]
  selected_genes <- selected_genes[selected_genes != "(Intercept)"]  # Remove intercept
  
  gene_list[[i]] <- selected_genes
}



# Frequency table
gene_freq <- sort(table(unlist(gene_list)), decreasing = TRUE)
gene_stability <- gene_freq / n_runs  # proportion

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

