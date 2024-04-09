
library(Seurat)
library(fgsea)
library(DESeq2)

# Sample 1 and Sample 2

# Load data
Sample_1 <- load("data/Sample_1.RData")
Sample_2 <- load("data/Sample_2.RData")

# Get gene expression data
genes_1 <- rownames(coculture_CD34)
genes_2 <- rownames(monoculture_CD34)
cells_1 <- WhichCells(coculture_CD34, expression = ((CD34 > 0) & (CD38 == 0)))
cells_2 <- WhichCells(monoculture_CD34, expression = ((CD34 > 0) & (CD38 == 0)))
#gene_expr_1 <- coculture_CD34@assays[["RNA"]]@counts[, cells_1]
gene_expr_1 <- coculture_CD34.data[["Gene Expression"]][, cells_1]
cond_1 <- matrix("coculture", nrow = ncol(gene_expr_1), ncol = 1)
#gene_expr_2 <- monoculture_CD34@assays[["RNA"]]@counts[, cells_2]
gene_expr_2 <- monoculture_CD34.data[["Gene Expression"]][, cells_2]
cond_2 <- matrix("monoculture", nrow = ncol(gene_expr_2), ncol = 1)

gene_expr_1_2 <- cbind(gene_expr_1, gene_expr_2)
cells_1_2 <- colnames(gene_expr_1_2)
cond_1_2 <- rbind(cond_1, cond_2)
cond_1_2 <- cbind(cells_1_2, cond_1_2)
colnames(cond_1_2) <- c("cell", "culture")

# Run DESeq
dds_1_2 <- DESeqDataSetFromMatrix(countData = gene_expr_1_2,
                                  colData = cond_1_2,
                                  design = ~ culture)
dds_1_2 <- DESeq(dds_1_2)

# Get results of DESeq
res_1_2 <- subset(results(dds_1_2), padj < 0.005)
res_1_2 <- res_1_2[complete.cases(res_1_2$padj), ]

res_1 <- subset(res_1_2, log2FoldChange < 0) #coculture overexpressed compared to monoculture
res_2 <- subset(res_1_2, log2FoldChange > 0) #monoculture overexpressed compared to coculture

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
stats_1 <- sort(stats_1)
stats_2 <- sort(stats_2)

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
fgsea_res_1 <- fgsea_res_1[order(fgsea_res_1$padj)]
fgsea_res_2 <- fgsea_res_2[order(fgsea_res_2$padj)]

head(fgsea_res_1, 10)
head(fgsea_res_2, 10)

label_1 <- matrix("coculture", nrow = nrow(fgsea_res_1), ncol = 1)
data_1 <- cbind(fgsea_res_1$pathway, label_1, fgsea_res_1$padj, fgsea_res_1$ES)

label_2 <- matrix("monoculture", nrow = nrow(fgsea_res_2), ncol = 1)
data_2 <- cbind(fgsea_res_2$pathway, label_2, fgsea_res_2$padj, fgsea_res_2$ES)

data_1_2 <- rbind(data_1, data_2)
colnames(data_1_2) <- c("pathway", "label", "padj", "enrichment score")
write.csv(data_1_2, file = "data_1_2.csv", row.names = FALSE)

# Sample 4 and Sample 5

# Load data
Sample_4 <- load("data/Sample_4.RData")
Sample_5 <- load("data/Sample_5.RData")

# Get gene expression data
genes_4 <- rownames(LowPctCytoCoculture)
genes_5 <- rownames(LowPctCytoMonoculture)
cells_4 <- WhichCells(LowPctCytoCoculture, expression = ((CD34 > 0) & (CD38 == 0)))
cells_5 <- WhichCells(LowPctCytoMonoculture, expression = ((CD34 > 0) & (CD38 == 0)))
gene_expr_4 <- LowPctCytoCoculture@assays[["RNA"]]@counts[, cells_4]
cond_4 <- matrix("low pct cyto coculture", nrow = ncol(gene_expr_4), ncol = 1)
gene_expr_5 <- LowPctCytoMonoculture@assays[["RNA"]]@counts[, cells_5]
cond_5 <- matrix("low pct cyto monoculture", nrow = ncol(gene_expr_5), ncol = 1)

gene_expr_4_5 <- cbind(gene_expr_4, gene_expr_5)
cells_4_5 <- colnames(gene_expr_4_5)
cond_4_5 <- rbind(cond_4, cond_5)
cond_4_5 <- cbind(cells_4_5, cond_4_5)
colnames(cond_4_5) <- c("cell", "culture")

# Run DESeq
dds_4_5 <- DESeqDataSetFromMatrix(countData = gene_expr_4_5,
                                  colData = cond_4_5,
                                  design = ~ culture)
dds_4_5 <- DESeq(dds_4_5)

# Get results of DESeq
res_4_5 <- subset(results(dds_4_5), padj < 0.005)
res_4_5 <- res_4_5[complete.cases(res$padj), ]

res_4 <- subset(res_4_5, log2FoldChange < 0) #coculture overexpressed compared to monoculture
res_5 <- subset(res_4_5, log2FoldChange > 0) #monoculture overexpressed compared to coculture

# Find DEGs
DEGs_4 <- rownames(res_4)
DEGs_5 <- rownames(res_5)
DEGs_4 <- DEGs_4[order(res_4$padj)]
DEGs_5 <- DEGs_5[order(res_5$padj)]

# Get stats for GSEA
stats_4 <- res_4$padj
stats_5 <- res_5$padj
names(stats_4) <- rownames(res_4)
names(stats_5) <- rownames(res_5)
stats_4 <- sort(stats_4)
stats_5 <- sort(stats_5)

# Get C2 pathways
c2_pathways <- gmtPathways("c2.all.v2023.2.Hs.symbols.gmt")

# Run GSEA
fgsea_res_4 <- fgsea(pathways = c2_pathways, 
                     stats    = stats_4,
                     minSize  = 0,
                     maxSize  = 500)

fgsea_res_5 <- fgsea(pathways = c2_pathways, 
                     stats    = stats_5,
                     minSize  = 0,
                     maxSize  = 500)

# Analyze GSEA results
fgsea_res_4 <- fgsea_res_4[order(fgsea_res_4$padj)]
fgsea_res_5 <- fgsea_res_5[order(fgsea_res_5$padj)]

head(fgsea_res_4, 10)
head(fgsea_res_5, 10)

label_4 <- matrix("coculture", nrow = nrow(fgsea_res_4), ncol = 1)
data_4 <- cbind(fgsea_res_4$pathway, label_4, fgsea_res_4$padj, fgsea_res_4$ES)

label_5 <- matrix("monoculture", nrow = nrow(fgsea_res_5), ncol = 1)
data_5 <- cbind(fgsea_res_5$pathway, label_5, fgsea_res_5$padj, fgsea_res_5$ES)

data_4_5 <- rbind(data_4, data_5)
colnames(data_4_5) <- c("pathway", "label", "padj", "enrichment score")
write.csv(data_4_5, file = "data_4_5.csv", row.names = FALSE)

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