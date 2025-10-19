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
genes_geo <- deg_geo_sig$SYMBOL

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
keep <- rowSums(expr_tcga > 10) >= 10
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
gene_stats <- sort(gene_stats, decreasing = TRUE)

gsea_kegg <- gseKEGG(geneList = gene_stats, organism = "hsa", pvalueCutoff = 0.2)
dotplot(gsea_kegg, showCategory = 15, title = "KEGG GSEA (TCGA-PAAD)")
