# Load raw data (no header, since first rows are custom)
file_path <- "paad.txt"
raw <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)


View(raw)


# Extract real sample IDs from first row (excluding gene_id column)
sample_ids <- as.character(unlist(raw[1, -1]))  # remove first column
colnames_data <- c("gene_id", sample_ids)

View(raw_data)


# Data starts from 3rd row onward — remove 1st (sample names) and 2nd row (repeated "normalized_counts")
expr_data <- raw[-c(1, 2), ]
colnames(expr_data) <- colnames_data


View(expr_data)




# Remove repeated "gene_id" in expression data rows
expr_data$gene_id <- sub(".*\\|", "", expr_data$gene_id)

# Drop duplicate gene names if any
expr_data <- expr_data[!duplicated(expr_data$gene_id), ]

# Fix: remove row that was originally the 3rd row in file (still "normalized_count" text)
expr_data <- expr_data[!expr_data$gene_id %in% c("normalized_count", "gene_id"), ]

# Set gene_id as rownames and drop that column
rownames(expr_data) <- expr_data$gene_id
expr_data <- expr_data[, -1]

View(expr_data)





# Convert expression values to numeric
expr_matrix <- as.data.frame(lapply(expr_data, as.numeric))
expr_matrix <- as.matrix(expr_matrix)
rownames(expr_matrix) <- rownames(expr_data)

View(expr_matrix)







unique(sample_types)


sample_ids <- colnames(expr_matrix)  # Save column names first

# Extract the 14th-15th characters (sample type codes)
sample_types <- substr(sample_ids, 14, 15)

# Double check lengths match
length(sample_types) == ncol(expr_matrix)  # should return TRUE





valid_types <- sample_types %in% c("01", "11")

# Now safely subset the matrix and other variables
expr_matrix <- expr_matrix[, valid_types]
sample_types <- sample_types[valid_types]
labels <- ifelse(sample_types == "01", 1, 0)
sample_ids <- colnames(expr_matrix)


length(labels)  == ncol(expr_matrix)      # should equal ncol(expr_matrix)
length(sample_types) ==ncol(expr_matrix)        # sam e
length(sample_ids)== ncol(expr_matrix)          # same


print(dim(expr_matrix))
print(length(sample_types))
print(length(valid_types))



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











#====================================================
# Step 6: Map Ensembl IDs to Gene Symbols
library(org.Hs.eg.db)

# Use Ensembl IDs from your DEGs
ensembl_ids <- rownames(deg_tcga_sig)

library(org.Hs.eg.db)

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,           # your Entrez IDs
                       column = "SYMBOL",
                       keytype = "ENTREZID",          # <- important change here
                       multiVals = "first")


# Add mapped symbols as new column
deg_tcga_sig$symbol <- gene_symbols

# Remove rows without gene symbol mapping
deg_tcga_sig <- deg_tcga_sig[!is.na(deg_tcga_sig$symbol), ]

# Optional: Remove duplicated gene symbols (keep highest significance)
deg_tcga_sig <- deg_tcga_sig[!duplicated(deg_tcga_sig$symbol), ]

# Check top results
head(deg_tcga_sig)

gene_tcga <- deg_tcga_sig$symbol


#-==============================================================================






library(GEOquery)

# Download the dataset
gse <- getGEO("GSE62452", GSEMatrix = TRUE)[[1]]

# Expression matrix (already log2 normalized)
expr_geo <- exprs(gse)

# Get phenotype info
pdata <- pData(gse)
group_geo <- pdata$`tissue type:ch1`

# Convert labels: tumor = 1, normal = 0
labels_geo <- ifelse(group_geo == "tumor", 1, 0)





#DEG ANALYSISS



library(limma)

gse <- getGEO("GSE62452", GSEMatrix = TRUE)[[1]]
expr_geo <- exprs(gse)
pheno <- pData(gse)

labels_geo <- ifelse(pheno$`tissue:ch1` == "Pancreatic tumor", 1, 0)


names(pheno)
table(pheno$`tissue:ch1`)
length(labels_geo) == ncol(expr_geo)#verify cheyan vendi



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

table(labels_geo)  # Should give counts of 0 and 1

#===============================================================










head(deg_geo_common$SYMBOL)
head(deg_tcga_common$symbol)
intersect(deg_geo_common$SYMBOL, deg_tcga_common$symbol)



# Clean gene symbols
deg_geo_common$SYMBOL <- toupper(trimws(deg_geo_common$SYMBOL))
deg_tcga_common$symbol <- toupper(trimws(deg_tcga_common$symbol))

# Remove any rows with NA in symbols
deg_geo_common <- deg_geo_common[!is.na(deg_geo_common$SYMBOL), ]
deg_tcga_common <- deg_tcga_common[!is.na(deg_tcga_common$symbol), ]

# Remove duplicates — very important!
deg_geo_common <- deg_geo_common[!duplicated(deg_geo_common$SYMBOL), ]
deg_tcga_common <- deg_tcga_common[!duplicated(deg_tcga_common$symbol), ]

# Set rownames
rownames(deg_geo_common) <- deg_geo_common$SYMBOL
rownames(deg_tcga_common) <- deg_tcga_common$symbol

# Now intersect again
common_genes <- intersect(rownames(deg_geo_common), rownames(deg_tcga_common))

# Subset final DEG tables
deg_geo_final <- deg_geo_common[common_genes, ]
deg_tcga_final <- deg_tcga_common[common_genes, ]

# Check
length(common_genes)










#================================================GO/kegg==================


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)





# Use common_genes or rownames(deg_tcga_final) or deg_geo_final if that is correct
# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(common_genes, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

# Remove NAs (if any)
gene_ids <- gene_entrez$ENTREZID


ego <- enrichGO(gene         = gene_ids,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)



ekegg <- enrichKEGG(gene         = gene_ids,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)


dotplot(ego, showCategory = 20) + ggtitle("GO Biological Process Enrichment")
dotplot(ekegg, showCategory = 20) + ggtitle("KEGG Pathway Enrichment")


barplot(ego, showCategory = 15, title = "GO BP Enrichment", drop = TRUE)
barplot(ekegg, showCategory = 15, title = "KEGG Pathway Enrichment", drop = TRUE)

# View top 10 GO terms
head(ego, 10)

# View top 10 KEGG pathways
head(ekegg, 10)

# Number of significant terms
nrow(ego)      # For GO
nrow(ekegg)    # For KEGG



# For GO
cnetplot(ego, showCategory = 5, circular = FALSE, colorEdge = TRUE)

# For KEGG
cnetplot(ekegg, showCategory = 5, circular = FALSE, colorEdge = TRUE)


heatplot(ego, showCategory = 10)
heatplot(ekegg, showCategory = 10)

go_result_df <- as.data.frame(ego)
kegg_result_df <- as.data.frame(ekegg)

#Print top 10 GO terms  
head(go_result_df[, c("ID", "Description", "p.adjust", "GeneRatio", "Count")], 10)
#Print top 10 KEGG pathways
head(kegg_result_df[, c("ID", "Description", "p.adjust", "GeneRatio", "Count")], 10)

# how many terms were significantly enriched:
sum(go_result_df$p.adjust < 0.05)
sum(kegg_result_df$p.adjust < 0.05)










#GESEA===================================================================
# Step 1.1: Clean TCGA DEG table (already done earlier as deg_tcga)
deg_tcga_ranks <- deg_tcga
deg_tcga_ranks <- deg_tcga_ranks[!is.na(deg_tcga_ranks$log2FoldChange), ]


str(deg_tcga)

# Extract logFC values
gene_list_tcga <- deg_tcga$logFC

# Set gene names from rownames
names(gene_list_tcga) <- toupper(trimws(rownames(deg_tcga)))

# Sort in decreasing order (important for GSEA)
gene_list_tcga <- sort(gene_list_tcga, decreasing = TRUE)

# Preview
head(gene_list_tcga)





#go based GSEA===========================================
# Step 1: Load library and DB
library(clusterProfiler)
library(org.Hs.eg.db)

# Step 2: Rank genes
deg_tcga_ranks <- deg_tcga[order(deg_tcga$logFC, decreasing = TRUE), ]
gene_symbols <- rownames(deg_tcga_ranks)  # these should be SYMBOLs
gene_list_tcga <- deg_tcga_ranks$logFC
names(gene_list_tcga) <- toupper(trimws(gene_symbols))  # force uppercase
gene_list_tcga <- sort(gene_list_tcga, decreasing = TRUE)

# Step 3: Run GSEA
gsea_go <- gseGO(
  geneList     = gene_list_tcga,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  pAdjustMethod = "BH",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = TRUE  # you can keep this TRUE to debug
)





# Convert to data frame
gsea_go_df <- as.data.frame(gsea_go)

# View top 10 results
head(gsea_go_df[, c("ID", "Description", "p.adjust", "NES", "core_enrichment")], 10)

