
# (1) Select cells from Sample_1 and Sample_2 whose CD34 gene 
#     expression is greater than 0 while CD38 gene expression is 0, 
#     and merge them into the same object/data matrix.

library(Seurat)
library(fgsea)
library(DESeq2)

Sample_1 <- load("data/Sample_1.RData")
genes_1 <- rownames(coculture_CD34)
cells_1 <- WhichCells(coculture_CD34, expression = ((CD34 > 0) & (CD38 == 0)))
gene_expr_1 <- coculture_CD34.data[["Gene Expression"]][, cells_1]
cells_1 <- DataFrame(cells_1)

Sample_2 <- load("data/Sample_2.RData")
genes_2 <- rownames(monoculture_CD34)
cells_2 <- WhichCells(monoculture_CD34, expression = ((CD34 > 0) & (CD38 == 0)))
gene_expr_2 <- monoculture_CD34.data[["Gene Expression"]][, cells_2]
cells_2 <- DataFrame(cells_2)

# cells = union(cells_1, cells_2)
# cells

# (2) Between the two groups of cells (one group from Sample_1, 
#     one group from Sample_2), run GSEA to check pathway enrichment 
#     in the overexpressed genes in each of the groups against the c2 
#     curated pathways (https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=C2). 
#     For this, we provide GSEA with the gene expression matrix, group labels, 
#     and pathway files containing gene members. For the gene ranking metric in GSEA, 
#     we can select the t-test.

# (3) Correct p-values for multiple tests.

dds_1 <- DESeqDataSetFromMatrix(countData = gene_expr_1,
                              colData = cells_1,
                              design = ~ 1)
dds_1 <- DESeq(dds_1)
res_1 <- results(dds_1)
res_1
res_1 <- subset(res_1, padj < 1.00e-06)
DEGs_1 <- rownames(res_1)
DEGs_1
DEGs_1 <- DEGs_1[order(res_1$padj)]
res_1 <- res_1[complete.cases(res_1$padj), ]
stats_1 <- res_1$padj
names(stats_1) <- rownames(res_1)
stats_1 <- sort(stats_1)
stats_1 <- head(stats_1, 10)

dds_2 <- DESeqDataSetFromMatrix(countData = gene_expr_2,
                                colData = cells_2,
                                design = ~ 1)
dds_2 <- DESeq(dds_2)
res_2 <- results(dds_2)
res_2
res_2 <- subset(res_2, padj < 1.00e-06)
DEGs_2 <- rownames(res_2)
DEGs_2
DEGs_2 <- DEGs_2[order(res_2$padj)]
res_2 <- res_2[complete.cases(res_2$padj), ]
stats_2 <- res_2$padj
names(stats_2) <- rownames(res_2)
stats_2 <- sort(stats_2)
stats_2 <- head(stats_2, 10)

c2_pathways <- gmtPathways("c2.all.v2023.2.Hs.symbols.gmt")

fgsea_res_1 <- fgsea(pathways = c2_pathways, 
                     stats    = stats_1,
                     minSize  = 0,
                     maxSize  = 500)

fgsea_res_2 <- fgsea(pathways = c2_pathways, 
                     stats    = stats_2,
                     minSize  = 0,
                     maxSize  = 500)

fgsea_res_1 <- subset(fgsea_res_1, padj < 1.00e-04)
head(fgsea_res_1, 100)
fgsea_res_1 <- fgsea_res_1[order(fgsea_res_1$padj)]

fgsea_res_2 <- subset(fgsea_res_2, padj < 1.00e-04)
head(fgsea_res_2, 100)
fgsea_res_2 <- fgsea_res_2[order(fgsea_res_2$padj)]


# (4) For the pathways listed in the document (below the green part), 
#     visualize their p-values in a pathway by group matrix, which will 
#     contain two columns, one for each group.

# (5) Repeat the above steps between Sample_4 and Sample_5.