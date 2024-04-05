
library(Seurat)
library(fgsea)
library(DESeq2)

# Select cells from Sample_1 and Sample_2 whose CD34 gene 
# expression is greater than 0 while CD38 gene expression is 0, 
# and merge them into the same data matrix.

# Load data
Sample_1 <- load("data/Sample_1.RData")
Sample_2 <- load("data/Sample_2.RData")

# Get gene expression data
genes_1 <- rownames(coculture_CD34)
genes_2 <- rownames(monoculture_CD34)
cells_1 <- WhichCells(coculture_CD34, expression = ((CD34 > 0) & (CD38 == 0)))
cells_2 <- WhichCells(monoculture_CD34, expression = ((CD34 > 0) & (CD38 == 0)))
gene_expr_1 <- coculture_CD34.data[["Gene Expression"]][, cells_1]
gene_expr_2 <- monoculture_CD34.data[["Gene Expression"]][, cells_2]

# Between the two groups of cells (one group from Sample_1, 
# one group from Sample_2), run GSEA to check pathway enrichment 
# in the overexpressed genes in each of the groups against the c2 
# curated pathways (https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=C2). 
# For this, we provide GSEA with the gene expression matrix, group labels, 
# and pathway files containing gene members. For the gene ranking metric in GSEA, 
# we can select the t-test.
# Correct p-values for multiple tests.

# Run DESeq
dds_1 <- DESeqDataSetFromMatrix(countData = gene_expr_1,
                              colData = DataFrame(cells_1),
                              design = ~ 1)
dds_2 <- DESeqDataSetFromMatrix(countData = gene_expr_2,
                                colData = DataFrame(cells_2),
                                design = ~ 1)
dds_1 <- DESeq(dds_1)
dds_2 <- DESeq(dds_2)

# Get results of DESeq
res_1 <- subset(results(dds_1), padj < 1.00e-06)
res_2 <- subset(results(dds_2), padj < 1.00e-06)
res_1 <- res_1[complete.cases(res_1$padj), ]
res_2 <- res_2[complete.cases(res_2$padj), ]

# Find DEGs
DEGs_1 <- rownames(res_1)
DEGs_2 <- rownames(res_2)
DEGs_1 <- DEGs_1[order(res_1$padj)]
DEGs_2 <- DEGs_2[order(res_2$padj)]

# Get stats for GSEA
stats_1 <- res_1$padj
stats_2 <- res_2$padj
names(stats_1) <- rownames(res_1)
names(stats_2) <- rownames(res_2)
stats_1 <- head(sort(stats_1), 10)
stats_2 <- head(sort(stats_2), 10)

# Get C2 pathways
c2_pathways <- gmtPathways("c2.all.v2023.2.Hs.symbols.gmt")

# Run GSEA
fgsea_res_1 <- fgsea(pathways = c2_pathways, 
                     stats    = stats_1,
                     minSize  = 0,
                     maxSize  = 500)

fgsea_res_2 <- fgsea(pathways = c2_pathways, 
                     stats    = stats_2,
                     minSize  = 0,
                     maxSize  = 500)

# Analyze GSEA results
fgsea_res_1 <- subset(fgsea_res_1, padj < 1.00e-04)
fgsea_res_2 <- subset(fgsea_res_2, padj < 1.00e-04)
fgsea_res_1 <- fgsea_res_1[order(fgsea_res_1$padj)]
fgsea_res_2 <- fgsea_res_2[order(fgsea_res_2$padj)]

head(fgsea_res_1, 100)
head(fgsea_res_2, 100)


# For the pathways listed in the document (below the green part), 
# visualize their p-values in a pathway by group matrix, which will 
# contain two columns, one for each group.

# Pathways:
# Glycolysis: KEGG_GLYCOLYSIS_GLUCONEOGENESIS
# Citrate Cycle: KEGG_CITRATE_CYCLE_TCA_CYCLE
# Oxidative Phosphorylation: KEGG_OXIDATIVE_PHOSPHORYLATION
# Translation: KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION
# Folate Biosynthesis: KEGG_FOLATE_BIOSYNTHESIS
# Pentose phosphate pathway: KEGG_PENTOSE_PHOSPHATE_PATHWAY
# Pyruvate metabolism: KEGG_PYRUVATE_METABOLISM
# DNA Replication: KEGG_DNA_REPLICATION
# Cell cycle: KEGG_CELL_CYCLE
# Glyoxylate and dicarboxylate metabolism: KEGG_GLYOXYLATE_AND_DICARBOXYLATE_METABOLISM
# Sulfur metabolism: KEGG_SULFUR_METABOLISM
# Arginine and proline metabolism: KEGG_ARGININE_AND_PROLINE_METABOLISM
# Galactose metabolism: KEGG_GALACTOSE_METABOLISM
# Amino-acyl tRNA biosynthesis: KEGG_AMINOACYL_TRNA_BIOSYNTHESIS
# Alanine, aspartate, and glutamate metabolism: KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM
# Histidine metabolism: KEGG_HISTIDINE_METABOLISM
# Tryptophan metabolism: KEGG_TRYPTOPHAN_METABOLISM
# Sphingolipid metabolism: KEGG_SPHINGOLIPID_METABOLISM
# Purine metabolism: KEGG_PURINE_METABOLISM
# Glycine, serine and threonine metabolism: KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM

# Mitochondrial translation: REACTOME_MITOCHONDRIAL_TRANSLATION
# Respiratory electron transport: REACTOME_RESPIRATORY_ELECTRON_TRANSPORT
# Metabolism of RNA: REACTOME_METABOLISM_OF_RNA
# IRE1alpha activates the chaperones: REACTOME_IRE1ALPHA_ACTIVATES_CHAPERONES
# Unfolded protein response: REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR
# Adaptive immune system: REACTOME_ADAPTIVE_IMMUNE_SYSTEM
# Interferon signaling: REACTOME_INTERFERON_SIGNALING
# Cell cycle mitotic: REACTOME_CELL_CYCLE_MITOTIC
# D-glutamine and D-glutamate metabolism: REACTOME_GLUTAMATE_AND_GLUTAMINE_METABOLISM

# XBP1(S) activates the chaperone genes
# Activation of repair proteins at DNA DSB
# G2M DNA damage checkpoint
# G2M checkpoints
# Vitamin B6 metabolism

# Repeat the above steps between Sample_4 and Sample_5.