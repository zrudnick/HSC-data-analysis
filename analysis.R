
# (1) Select cells from Sample_1 and Sample_2 whose CD34 gene 
#     expression is greater than 0 while CD38 gene expression is 0, 
#     and merge them into the same object/data matrix.

library(Seurat)
library(fgsea)

Sample_1 <- load("data/Sample_1.RData")
genes_1 <- rownames(coculture_CD34)
cells_1 <- WhichCells(coculture_CD34, expression = ((CD34 > 0) & (CD38 == 0)))
cells_1
gene_expr_1 <- coculture_CD34.data[["Gene Expression"]][, cells_1]


Sample_2 <- load("data/Sample_2.RData")
genes_2 <- rownames(monoculture_CD34)
cells_2 <- WhichCells(monoculture_CD34, expression = ((CD34 > 0) & (CD38 == 0)))
cells_2
gene_expr_2 <- monoculture_CD34.data[["Gene Expression"]][, cells_2]

cells = union(cells_1, cells_2)
cells

# (2) Between the two groups of cells (one group from Sample_0, 
#     one group from Sample_1), run GSEA to check pathway enrichment 
#     in the overexpressed genes in each of the groups against the c2 
#     curated pathways (https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=C2). 
#     For this, we provide GSEA with the gene expression matrix, group labels, 
#     and pathway files containing gene members. For the gene ranking metric in GSEA, 
#     we can select the t-test.

data(examplePathways)
data(exampleRanks)
c2_pathways <- gmtPathways("c2.all.v2023.2.Hs.symbols.gmt")

fgsea_res_1 <- fgsea(pathways = c2_pathways, 
                     stats    = gene_expr_1,
                     minSize  = 0,
                     maxSize  = 500)

fgsea_res_2 <- fgsea(pathways = c2_pathways, 
                     stats    = gene_expr_2,
                     minSize  = 0,
                     maxSize  = 500)

head(fgseaRes)

# (3) Correct p-values for multiple tests.

# (4) For the pathways listed in the document (below the green part), 
#     visualize their p-values in a pathway by group matrix, which will 
#     contain two columns, one for each group.

# (5) Repeat the above steps between Sample_4 and Sample_5.