if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma"))

#----------------------------GEO------------------------------------------------

library(GEOquery)
library(limma)


gse <- getGEO("GSE62452", GSEMatrix = TRUE)
length(gse)
names(gse)
class(gse[[1]])
View(gse[[1]])



#  expression matrix
expr_data <- exprs(gse[[1]])
dim(expr_data)
View(expr_data)

# sample metadata
pheno <- pData(gse[[1]])
View(head(pheno))

#exploring
colnames(fData(gse[[1]]))
colnames(pData(gse[[1]]))


#  phenotype labels ir is to chk tumor and non tumor
table(pheno$`tissue:ch1`) 

# Create binary labels: 1 = tumor, 0 = normal
labels_geo <- ifelse(pheno$`tissue:ch1` == "Pancreatic tumor", 1, 0)

#to check if log2-normalized
summary(expr_data)
boxplot(expr_data,main="Pancretic cancer dataset",las=2)


#------------------------------TCGA---------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")


library(TCGAbiolinks)
library(SummarizedExperiment)




# Create a query for HTSeq - Counts (STAR aligned)
query_exp <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)


# Download data (uses GDC API) but dont run again if u have it already
GDCdownload(query_exp, method = "api", files.per.chunk = 20)



# Prepare expression matrix and metadata
data_exp <- GDCprepare(query = query_exp)

# Extract count matrix
expr_tcga <- assay(data_exp)

# Extract clinical info
clinical_tcga <- colData(data_exp)

# View structure
dim(expr_tcga)
head(clinical_tcga$shortLetterCode)



# TP = Primary Tumor, NT = Solid Tissue Normal
labels_tcga <- ifelse(clinical_tcga$shortLetterCode == "TP", 1, 0)

# Confirm balance
table(labels_tcga)

# Keep genes with counts >10 in at least 10 samples
keep <- rowSums(expr_tcga > 10) >= 10
expr_tcga_filtered <- expr_tcga[keep, ]

#-------------------------------------------------------------------------------



if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma")

library(limma)

# expr_data_filtered (genes x samples)
# labels_geo (0 = normal, 1 = tumor)

group_geo <- factor(labels_geo)
design_geo <- model.matrix(~group_geo)


# Use correct GEO expression data and labels
expr_geo <- exprs(gse[[1]])
labels_geo <- ifelse(pheno$`tissue:ch1` == "Pancreatic tumor", 1, 0)

# Fit design and model
group_geo <- factor(labels_geo)
design_geo <- model.matrix(~group_geo)

fit_geo <- lmFit(expr_geo, design_geo)
fit_geo <- eBayes(fit_geo)

# Get full DEG list
deg_geo <- topTable(fit_geo, coef = 2, number = Inf, adjust.method = "fdr")

# ✅ Define significant DEGs
deg_geo_sig <- deg_geo[deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 1, ]





BiocManager::install("DESeq2")
library(DESeq2)



# Expression matrix
# expr_tcga_filtered (genes x samples)
# labels_tcga (0 = normal, 1 = tumor)

unique(labels_tcga)


# Make colData
col_data <- data.frame(
  row.names = colnames(expr_tcga_filtered),
  condition = factor(labels_tcga)
)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = expr_tcga_filtered,
  colData = col_data,
  design = ~ condition
)



dds <- DESeq(dds)
res <- results(dds)

# Order by p-value
res_ordered <- res[order(res$padj), ]

# DEG Filtering - TCGA
deg_tcga_sig <- subset(res_ordered, padj < 0.05 & abs(log2FoldChange) > 0.5)

# Clean Ensembl IDs (remove version suffix)
clean_ensembl_ids <- gsub("\\..*", "", rownames(deg_tcga_sig))

# Map Ensembl IDs to gene symbols
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

tcga_symbols <- mapIds(
  org.Hs.eg.db,
  keys = clean_ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
tcga_symbols <- na.omit(tcga_symbols)
genes_tcga <- unname(tcga_symbols)  # Remove names, keep symbols only


# DEG Filtering - GEO
deg_geo_sig <- deg_geo[deg_geo$adj.P.Val < 0.05 & abs(deg_geo$logFC) > 0.5, ]

# Map probe IDs to gene symbols using GPL platform (GPL6244)
gpl <- getGEO("GPL6244", AnnotGPL = TRUE)
gpl_table <- Table(gpl)

# Clean probe-to-gene mapping
probe2gene <- gpl_table[, c("ID", "Gene symbol")]
colnames(probe2gene) <- c("PROBEID", "SYMBOL")
probe2gene$SYMBOL <- sapply(strsplit(probe2gene$SYMBOL, "///", fixed = TRUE), `[`, 1)
probe2gene <- probe2gene[!is.na(probe2gene$SYMBOL) & probe2gene$SYMBOL != "", ]

# Add gene symbols to GEO DEGs
deg_geo_sig$SYMBOL <- probe2gene$SYMBOL[match(rownames(deg_geo_sig), probe2gene$PROBEID)]
deg_geo_sig <- deg_geo_sig[!is.na(deg_geo_sig$SYMBOL) & deg_geo_sig$SYMBOL != "", ]
genes_geo <- deg_geo_sig$SYMBOL




# 🧬 Identify Common DEGs between GEO and TCGA
common_genes <- intersect(genes_geo, genes_tcga)

# View
length(common_genes)
head(common_genes)



#-------------------------------------------------------------------------------
#Functional Enrichment (GO & KEGG) of Common DEGs

# Install if not done yet
BiocManager::install("GO.db")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

# Convert gene symbols to Entrez IDs fr the clusterprrofiler to work
entrez_ids <- bitr(
  common_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# mapped genes
head(entrez_ids)





#enrichment
ego <- enrichGO(
  gene = entrez_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",  # You can also use "BP", "MF", or "CC"
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(ego)




#KEGG 
ekegg <- enrichKEGG(
  gene = entrez_ids$ENTREZID,
  organism = 'hsa',
  pvalueCutoff = 0.05
)

head(ekegg)
head(ekegg@result[, c("ID", "Description", "pvalue", "geneID")], 20)





#visual of top enriched terms
# GO bar plot
barplot(ego, showCategory = 15, title = "GO Enrichment")

# KEGG dot plot
dotplot(ekegg, showCategory = 15, title = "KEGG Pathway Enrichment")



#                      ^
#                      |
#the resut above from the plot is not what we should be getting i think
#need more tweaking above the codes













#-------------------------------------------------------------------------------
#-------------------------------TRoubleshooting---------------------------------



#-----------Pca focused plots
# Reuse your existing entrez_ids
# entrez_ids <- bitr(...)


# Pathways of interest: Pancreatic & general cancer-related ones
cancer_pathways <- c("hsa05212", "hsa05200", "hsa05206", "hsa05205", "hsa05204") 

# Filter KEGG result for only these
ekegg_filtered <- ekegg@result[ekegg@result$ID %in% cancer_pathways, ]

# Optional: Clean if empty
if (nrow(ekegg_filtered) == 0) {
  message("No pancreatic cancer-related pathways found in enrichment results.")
} else {
  # Convert to enrichmentResult object for plotting
  ekegg_filtered_obj <- enrichResult(ekegg_filtered)
  
  # Visualize filtered KEGG pathways
  dotplot(ekegg_filtered_obj, showCategory = 10, title = "Pancreatic Cancer KEGG Enrichment")
}





# Check if hsa05212 is present at all (regardless of p-value)
ekegg_df <- ekegg@result
ekegg_df[ekegg_df$ID == "hsa05212", ]

known_pca_genes <- c("KRAS", "CDKN2A", "TP53", "SMAD4", "BRCA2", "ARID1A")

intersect(common_genes, known_pca_genes)








# Use TCGA only (for KEGG)
entrez_tcga <- bitr(genes_tcga, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ekegg_tcga <- enrichKEGG(gene = entrez_tcga$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
dotplot(ekegg_tcga, showCategory = 15, title = "TCGA-PAAD KEGG Enrichment")












# Start fresh from res_ordered and clean_ensembl_ids
gene_stats <- res_ordered$log2FoldChange
names(gene_stats) <- clean_ensembl_ids

# Map Ensembl to Entrez
ensembl2entrez <- mapIds(
  org.Hs.eg.db,
  keys = clean_ensembl_ids,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Final clean: Remove any NA-named or NA-valued entries
gene_stats <- gene_stats[!is.na(names(gene_stats)) & !is.na(gene_stats)]

# Also ensure unique names (GSEA requirement)
gene_stats <- gene_stats[!duplicated(names(gene_stats))]

# Sort in decreasing order
gene_stats <- sort(gene_stats, decreasing = TRUE)
gsea_kegg <- gseKEGG(geneList = gene_stats, organism = "hsa", pvalueCutoff = 0.2)
dotplot(gsea_kegg, showCategory = 15, title = "KEGG GSEA (TCGA-PAAD)")


