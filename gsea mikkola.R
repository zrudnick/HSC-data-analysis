library(Seurat)
library(fgsea)
library(DESeq2)
library(tibble)

# Mikkola

##########################Load data
load("data/Mikkola_Fig2_23.RData")

# Get gene expression data
sample = UpdateSeuratObject(object = sample) # fix error with image
genes <- rownames(sample)
cells <- sample[genes]

X <- as.matrix(cells[["RNA"]]$counts)

orig.ident <- cells@meta.data[["orig.ident"]]

##########################Run DESeq

dds <- DESeqDataSetFromMatrix(countData = X,
                                  colData = DataFrame(orig.ident),
                                  design = ~ orig.ident)
dds <- DESeq(dds)

# Get results of DESeq:
#agm-4wk-658 vs others
agm_4wk_658 <- results(dds, contrast=c(1, -1/7, -1/7, -1/7, -1/7, -1/7, -1/7, -1/7), tidy=TRUE)
agm_4wk_658_output <- dplyr::select(agm_4wk_658, row, log2FoldChange, padj)
agm_4wk_658_output <- na.omit(agm_4wk_658_output)
agm_4wk_658 <- dplyr::select(agm_4wk_658, row, stat)
agm_4wk_658 <- na.omit(agm_4wk_658)

#agm-5wk-555 vs others
agm_5wk_555 <- results(dds, contrast=c(-1/7, 1, -1/7, -1/7, -1/7, -1/7, -1/7, -1/7), tidy=TRUE)
agm_5wk_555_output <- dplyr::select(agm_5wk_555, row, log2FoldChange, padj)
agm_5wk_555_output <- na.omit(agm_5wk_555_output)
agm_5wk_555 <- dplyr::select(agm_5wk_555, row, stat)
agm_5wk_555 <- na.omit(agm_5wk_555)

#agm-5wk-575 vs others
agm_5wk_575 <- results(dds, contrast=c(-1/7, -1/7, 1, -1/7, -1/7, -1/7, -1/7, -1/7), tidy=TRUE)
agm_5wk_575_output <- dplyr::select(agm_5wk_575, row, log2FoldChange, padj)
agm_5wk_575_output <- na.omit(agm_5wk_575_output)
agm_5wk_575 <- dplyr::select(agm_5wk_575, row, stat)
agm_5wk_575 <- na.omit(agm_5wk_575)

#liver-6wk-563 vs others
liver_6wk_563 <- results(dds, contrast=c(-1/7, -1/7, -1/7, 1, -1/7, -1/7, -1/7, -1/7), tidy=TRUE)
liver_6wk_563_output <- dplyr::select(liver_6wk_563, row, log2FoldChange, padj)
liver_6wk_563_output <- na.omit(liver_6wk_563_output)
liver_6wk_563 <- dplyr::select(liver_6wk_563, row, stat)
liver_6wk_563 <- na.omit(liver_6wk_563)

#liver-8wk-553 vs others
liver_8wk_553 <- results(dds, contrast=c(-1/7, -1/7, -1/7, -1/7, 1, -1/7, -1/7, -1/7), tidy=TRUE)
liver_8wk_553_output <- dplyr::select(liver_8wk_553, row, log2FoldChange, padj)
liver_8wk_553_output <- na.omit(liver_8wk_553_output)
liver_8wk_553 <- dplyr::select(liver_8wk_553, row, stat)
liver_8wk_553 <- na.omit(liver_8wk_553)

#liver-11wk-569 vs others
liver_11wk_569 <- results(dds, contrast=c(-1/7, -1/7, -1/7, -1/7, -1/7, 1, -1/7, -1/7), tidy=TRUE)
liver_11wk_569_output <- dplyr::select(liver_11wk_569, row, log2FoldChange, padj)
liver_11wk_569_output <- na.omit(liver_11wk_569_output)
liver_11wk_569 <- dplyr::select(liver_11wk_569, row, stat)
liver_11wk_569 <- na.omit(liver_11wk_569)

#liver-15wk-101 vs others
liver_15wk_101 <- results(dds, contrast=c(-1/7, -1/7, -1/7, -1/7, -1/7, -1/7, 1, -1/7), tidy=TRUE)
liver_15wk_101_output <- dplyr::select(liver_15wk_101, row, log2FoldChange, padj)
liver_15wk_101_output <- na.omit(liver_15wk_101_output)
liver_15wk_101 <- dplyr::select(liver_15wk_101, row, stat)
liver_15wk_101 <- na.omit(liver_15wk_101)

#cb-40wk-201 vs others
cb_40wk_201 <- results(dds, contrast=c(-1/7, -1/7, -1/7, -1/7, -1/7, -1/7, -1/7, 1), tidy=TRUE)
cb_40wk_201_output <- dplyr::select(cb_40wk_201, row, log2FoldChange, padj)
cb_40wk_201_output <- na.omit(cb_40wk_201_output)
cb_40wk_201 <- dplyr::select(cb_40wk_201, row, stat)
cb_40wk_201 <- na.omit(cb_40wk_201)

##########################Run GSEA
c2_pathways <- gmtPathways("c2.all.v2023.2.Hs.symbols.gmt")

#run gsea for gene rankings: agm_4wk_658 vs others
ranks <- deframe(agm_4wk_658)
head(ranks, 20)

fgsea_agm_4wk_658 <- fgsea(pathways = c2_pathways,
                           stats = ranks,
                           minSize = 5,
                           maxSize = 500,
                           nperm = 10000)

fgsea_agm_4wk_658_tidy <- as_tibble(fgsea_agm_4wk_658)
fgsea_agm_4wk_658_tidy <- dplyr::arrange(fgsea_agm_4wk_658_tidy, desc(NES))
fgsea_agm_4wk_658_tidy <- dplyr::select(fgsea_agm_4wk_658_tidy, pathway, pval, padj, NES)

#Collect deseq results for each pathway
fgsea_agm_4wk_658_tidy$genes <- NA
fgsea_agm_4wk_658_tidy$log2FoldChange <- NA
fgsea_agm_4wk_658_tidy$deseq_padj <- NA
count <- 0
for (pathway_value in fgsea_agm_4wk_658_tidy$pathway) {
  count <- count + 1
  selected_rows <- agm_4wk_658_output[agm_4wk_658_output$row %in% c2_pathways[[pathway_value]] & agm_4wk_658_output$padj<0.05 & agm_4wk_658_output$log2FoldChange>0,]
  selected_rows <- dplyr::arrange(selected_rows, desc(log2FoldChange))
  genes_cvm <- paste(selected_rows$row, collapse = ';')
  fold_cvm <- paste(selected_rows$log2FoldChange, collapse = ';')
  pv_cvm <- paste(selected_rows$padj, collapse = ';')
  fgsea_agm_4wk_658_tidy$genes[count] <- genes_cvm
  fgsea_agm_4wk_658_tidy$log2FoldChange[count] <- fold_cvm
  fgsea_agm_4wk_658_tidy$deseq_padj[count] <- pv_cvm
}

#run gsea for gene rankings: agm_5wk_555 vs others
ranks <- deframe(agm_5wk_555)
head(ranks, 20)

fgsea_agm_5wk_555 <- fgsea(pathways = c2_pathways,
                           stats = ranks,
                           minSize = 5,
                           maxSize = 500,
                           nperm = 10000)

fgsea_agm_5wk_555_tidy <- as_tibble(fgsea_agm_5wk_555)
fgsea_agm_5wk_555_tidy <- dplyr::arrange(fgsea_agm_5wk_555_tidy, desc(NES))
fgsea_agm_5wk_555_tidy <- dplyr::select(fgsea_agm_5wk_555_tidy, pathway, pval, padj, NES)

#Collect deseq results for each pathway
fgsea_agm_5wk_555_tidy$genes <- NA
fgsea_agm_5wk_555_tidy$log2FoldChange <- NA
fgsea_agm_5wk_555_tidy$deseq_padj <- NA
count <- 0
for (pathway_value in fgsea_agm_5wk_555_tidy$pathway) {
  count <- count + 1
  selected_rows <- agm_5wk_555_output[agm_5wk_555_output$row %in% c2_pathways[[pathway_value]] & agm_5wk_555_output$padj<0.05 & agm_5wk_555_output$log2FoldChange>0,]
  selected_rows <- dplyr::arrange(selected_rows, desc(log2FoldChange))
  genes_cvm <- paste(selected_rows$row, collapse = ';')
  fold_cvm <- paste(selected_rows$log2FoldChange, collapse = ';')
  pv_cvm <- paste(selected_rows$padj, collapse = ';')
  fgsea_agm_5wk_555_tidy$genes[count] <- genes_cvm
  fgsea_agm_5wk_555_tidy$log2FoldChange[count] <- fold_cvm
  fgsea_agm_5wk_555_tidy$deseq_padj[count] <- pv_cvm
}

#run gsea for gene rankings: agm_5wk_575 vs others
ranks <- deframe(agm_5wk_575)
head(ranks, 20)

fgsea_agm_5wk_575 <- fgsea(pathways = c2_pathways,
                           stats = ranks,
                           minSize = 5,
                           maxSize = 500,
                           nperm = 10000)

fgsea_agm_5wk_575_tidy <- as_tibble(fgsea_agm_5wk_575)
fgsea_agm_5wk_575_tidy <- dplyr::arrange(fgsea_agm_5wk_575_tidy, desc(NES))
fgsea_agm_5wk_575_tidy <- dplyr::select(fgsea_agm_5wk_575_tidy, pathway, pval, padj, NES)

#Collect deseq results for each pathway
fgsea_agm_5wk_575_tidy$genes <- NA
fgsea_agm_5wk_575_tidy$log2FoldChange <- NA
fgsea_agm_5wk_575_tidy$deseq_padj <- NA
count <- 0
for (pathway_value in fgsea_agm_5wk_575_tidy$pathway) {
  count <- count + 1
  selected_rows <- agm_5wk_575_output[agm_5wk_575_output$row %in% c2_pathways[[pathway_value]] & agm_5wk_575_output$padj<0.05 & agm_5wk_575_output$log2FoldChange>0,]
  selected_rows <- dplyr::arrange(selected_rows, desc(log2FoldChange))
  genes_cvm <- paste(selected_rows$row, collapse = ';')
  fold_cvm <- paste(selected_rows$log2FoldChange, collapse = ';')
  pv_cvm <- paste(selected_rows$padj, collapse = ';')
  fgsea_agm_5wk_575_tidy$genes[count] <- genes_cvm
  fgsea_agm_5wk_575_tidy$log2FoldChange[count] <- fold_cvm
  fgsea_agm_5wk_575_tidy$deseq_padj[count] <- pv_cvm
}

#run gsea for gene rankings: liver_6wk_563 vs others
ranks <- deframe(liver_6wk_563)
head(ranks, 20)

fgsea_liver_6wk_563 <- fgsea(pathways = c2_pathways,
                             stats = ranks,
                             minSize = 5,
                             maxSize = 500,
                             nperm = 10000)

fgsea_liver_6wk_563_tidy <- as_tibble(fgsea_liver_6wk_563)
fgsea_liver_6wk_563_tidy <- dplyr::arrange(fgsea_liver_6wk_563_tidy, desc(NES))
fgsea_liver_6wk_563_tidy <- dplyr::select(fgsea_liver_6wk_563_tidy, pathway, pval, padj, NES)

#Collect deseq results for each pathway
fgsea_liver_6wk_563_tidy$genes <- NA
fgsea_liver_6wk_563_tidy$log2FoldChange <- NA
fgsea_liver_6wk_563_tidy$deseq_padj <- NA
count <- 0
for (pathway_value in fgsea_liver_6wk_563_tidy$pathway) {
  count <- count + 1
  selected_rows <- liver_6wk_563_output[liver_6wk_563_output$row %in% c2_pathways[[pathway_value]] & liver_6wk_563_output$padj<0.05 & liver_6wk_563_output$log2FoldChange>0,]
  selected_rows <- dplyr::arrange(selected_rows, desc(log2FoldChange))
  genes_cvm <- paste(selected_rows$row, collapse = ';')
  fold_cvm <- paste(selected_rows$log2FoldChange, collapse = ';')
  pv_cvm <- paste(selected_rows$padj, collapse = ';')
  fgsea_liver_6wk_563_tidy$genes[count] <- genes_cvm
  fgsea_liver_6wk_563_tidy$log2FoldChange[count] <- fold_cvm
  fgsea_liver_6wk_563_tidy$deseq_padj[count] <- pv_cvm
}

#run gsea for gene rankings: liver_8wk_553 vs others
ranks <- deframe(liver_8wk_553)
head(ranks, 20)

fgsea_liver_8wk_553 <- fgsea(pathways = c2_pathways,
                             stats = ranks,
                             minSize = 5,
                             maxSize = 500,
                             nperm = 10000)

fgsea_liver_8wk_553_tidy <- as_tibble(fgsea_liver_8wk_553)
fgsea_liver_8wk_553_tidy <- dplyr::arrange(fgsea_liver_8wk_553_tidy, desc(NES))
fgsea_liver_8wk_553_tidy <- dplyr::select(fgsea_liver_8wk_553_tidy, pathway, pval, padj, NES)

#Collect deseq results for each pathway
fgsea_liver_8wk_553_tidy$genes <- NA
fgsea_liver_8wk_553_tidy$log2FoldChange <- NA
fgsea_liver_8wk_553_tidy$deseq_padj <- NA
count <- 0
for (pathway_value in fgsea_liver_8wk_553_tidy$pathway) {
  count <- count + 1
  selected_rows <- liver_8wk_553_output[liver_8wk_553_output$row %in% c2_pathways[[pathway_value]] & liver_8wk_553_output$padj<0.05 & liver_8wk_553_output$log2FoldChange>0,]
  selected_rows <- dplyr::arrange(selected_rows, desc(log2FoldChange))
  genes_cvm <- paste(selected_rows$row, collapse = ';')
  fold_cvm <- paste(selected_rows$log2FoldChange, collapse = ';')
  pv_cvm <- paste(selected_rows$padj, collapse = ';')
  fgsea_liver_8wk_553_tidy$genes[count] <- genes_cvm
  fgsea_liver_8wk_553_tidy$log2FoldChange[count] <- fold_cvm
  fgsea_liver_8wk_553_tidy$deseq_padj[count] <- pv_cvm
}

#run gsea for gene rankings: liver_11wk_569 vs others
ranks <- deframe(liver_11wk_569)
head(ranks, 20)

fgsea_liver_11wk_569 <- fgsea(pathways = c2_pathways,
                              stats = ranks,
                              minSize = 5,
                              maxSize = 500,
                              nperm = 10000)

fgsea_liver_11wk_569_tidy <- as_tibble(fgsea_liver_11wk_569)
fgsea_liver_11wk_569_tidy <- dplyr::arrange(fgsea_liver_11wk_569_tidy, desc(NES))
fgsea_liver_11wk_569_tidy <- dplyr::select(fgsea_liver_11wk_569_tidy, pathway, pval, padj, NES)

#Collect deseq results for each pathway
fgsea_liver_11wk_569_tidy$genes <- NA
fgsea_liver_11wk_569_tidy$log2FoldChange <- NA
fgsea_liver_11wk_569_tidy$deseq_padj <- NA
count <- 0
for (pathway_value in fgsea_liver_11wk_569_tidy$pathway) {
  count <- count + 1
  selected_rows <- liver_11wk_569_output[liver_11wk_569_output$row %in% c2_pathways[[pathway_value]] & liver_11wk_569_output$padj<0.05 & liver_11wk_569_output$log2FoldChange>0,]
  selected_rows <- dplyr::arrange(selected_rows, desc(log2FoldChange))
  genes_cvm <- paste(selected_rows$row, collapse = ';')
  fold_cvm <- paste(selected_rows$log2FoldChange, collapse = ';')
  pv_cvm <- paste(selected_rows$padj, collapse = ';')
  fgsea_liver_11wk_569_tidy$genes[count] <- genes_cvm
  fgsea_liver_11wk_569_tidy$log2FoldChange[count] <- fold_cvm
  fgsea_liver_11wk_569_tidy$deseq_padj[count] <- pv_cvm
}

#run gsea for gene rankings: liver_15wk_101 vs others
ranks <- deframe(liver_15wk_101)
head(ranks, 20)

fgsea_liver_15wk_101 <- fgsea(pathways = c2_pathways,
                              stats = ranks,
                              minSize = 5,
                              maxSize = 500,
                              nperm = 10000)

fgsea_liver_15wk_101_tidy <- as_tibble(fgsea_liver_15wk_101)
fgsea_liver_15wk_101_tidy <- dplyr::arrange(fgsea_liver_15wk_101_tidy, desc(NES))
fgsea_liver_15wk_101_tidy <- dplyr::select(fgsea_liver_15wk_101_tidy, pathway, pval, padj, NES)

#Collect deseq results for each pathway
fgsea_liver_15wk_101_tidy$genes <- NA
fgsea_liver_15wk_101_tidy$log2FoldChange <- NA
fgsea_liver_15wk_101_tidy$deseq_padj <- NA
count <- 0
for (pathway_value in fgsea_liver_15wk_101_tidy$pathway) {
  count <- count + 1
  selected_rows <- liver_15wk_101_output[liver_15wk_101_output$row %in% c2_pathways[[pathway_value]] & liver_15wk_101_output$padj<0.05 & liver_15wk_101_output$log2FoldChange>0,]
  selected_rows <- dplyr::arrange(selected_rows, desc(log2FoldChange))
  genes_cvm <- paste(selected_rows$row, collapse = ';')
  fold_cvm <- paste(selected_rows$log2FoldChange, collapse = ';')
  pv_cvm <- paste(selected_rows$padj, collapse = ';')
  fgsea_liver_15wk_101_tidy$genes[count] <- genes_cvm
  fgsea_liver_15wk_101_tidy$log2FoldChange[count] <- fold_cvm
  fgsea_liver_15wk_101_tidy$deseq_padj[count] <- pv_cvm
}

#run gsea for gene rankings: cb_40wk_201 vs others
ranks <- deframe(cb_40wk_201)
head(ranks, 20)

fgsea_cb_40wk_201 <- fgsea(pathways = c2_pathways,
                           stats = ranks,
                           minSize = 5,
                           maxSize = 500,
                           nperm = 10000)

fgsea_cb_40wk_201_tidy <- as_tibble(fgsea_cb_40wk_201)
fgsea_cb_40wk_201_tidy <- dplyr::arrange(fgsea_cb_40wk_201_tidy, desc(NES))
fgsea_cb_40wk_201_tidy <- dplyr::select(fgsea_cb_40wk_201_tidy, pathway, pval, padj, NES)

#Collect deseq results for each pathway
fgsea_cb_40wk_201_tidy$genes <- NA
fgsea_cb_40wk_201_tidy$log2FoldChange <- NA
fgsea_cb_40wk_201_tidy$deseq_padj <- NA
count <- 0
for (pathway_value in fgsea_cb_40wk_201_tidy$pathway) {
  count <- count + 1
  selected_rows <- cb_40wk_201_output[cb_40wk_201_output$row %in% c2_pathways[[pathway_value]] & cb_40wk_201_output$padj<0.05 & cb_40wk_201_output$log2FoldChange>0,]
  selected_rows <- dplyr::arrange(selected_rows, desc(log2FoldChange))
  genes_cvm <- paste(selected_rows$row, collapse = ';')
  fold_cvm <- paste(selected_rows$log2FoldChange, collapse = ';')
  pv_cvm <- paste(selected_rows$padj, collapse = ';')
  fgsea_cb_40wk_201_tidy$genes[count] <- genes_cvm
  fgsea_cb_40wk_201_tidy$log2FoldChange[count] <- fold_cvm
  fgsea_cb_40wk_201_tidy$deseq_padj[count] <- pv_cvm
}

#join eight tables
names(fgsea_agm_4wk_658_tidy)[-1] <- paste0(names(fgsea_agm_4wk_658_tidy)[-1], "_agm_4wk_658_vs_others")
names(fgsea_agm_5wk_555_tidy)[-1] <- paste0(names(fgsea_agm_5wk_555_tidy)[-1], "_agm_5wk_555_vs_others")
names(fgsea_agm_5wk_575_tidy)[-1] <- paste0(names(fgsea_agm_5wk_575_tidy)[-1], "_agm_5wk_575_vs_others")
names(fgsea_liver_6wk_563_tidy)[-1] <- paste0(names(fgsea_liver_6wk_563_tidy)[-1], "_liver_6wk_563_vs_others")
names(fgsea_liver_8wk_553_tidy)[-1] <- paste0(names(fgsea_liver_8wk_553_tidy)[-1], "_liver_8wk_553_vs_others")
names(fgsea_liver_11wk_569_tidy)[-1] <- paste0(names(fgsea_liver_11wk_569_tidy)[-1], "_liver_11wk_569_vs_others")
names(fgsea_liver_15wk_101_tidy)[-1] <- paste0(names(fgsea_liver_15wk_101_tidy)[-1], "_liver_15wk_101_vs_others")
names(fgsea_cb_40wk_201_tidy)[-1] <- paste0(names(fgsea_cb_40wk_201_tidy)[-1], "_cb_40wk_201_vs_others")

joined <- dplyr::left_join(fgsea_agm_4wk_658_tidy, fgsea_agm_5wk_555_tidy, by = "pathway")
joined <- dplyr::left_join(joined, fgsea_agm_5wk_575_tidy, by = "pathway")
joined <- dplyr::left_join(joined, fgsea_liver_6wk_563_tidy, by = "pathway")
joined <- dplyr::left_join(joined, fgsea_liver_8wk_553_tidy, by = "pathway")
joined <- dplyr::left_join(joined, fgsea_liver_11wk_569_tidy, by = "pathway")
joined <- dplyr::left_join(joined, fgsea_liver_15wk_101_tidy, by = "pathway")
joined <- dplyr::left_join(joined, fgsea_cb_40wk_201_tidy, by = "pathway")

joined <- na.omit(joined)

write.csv(joined, "gsea_mikkola.csv", row.names = FALSE)
