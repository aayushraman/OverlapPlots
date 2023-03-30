#############################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 30th March 2023
#
# Program is used for:
# 1. Overlap analyses for all the Boxer et al Dataset
#############################################################################

## packages 
rm(list = ls())
setwd("~/Desktop/LongGene_Part-II/OverlapPlots/src/")
source("libraries.R")
 
## annotation dataset
load(file = "../dat-info/mm10_ncbi-refSeqGene_Dec2019.RData")
genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))

## KO/WT whole cell dataset --> 2922 DEGs
dat <- read.table("../dat/counts/GSE128178_10WT_10MeCP2_KO_whole_cell_RNAseq_exon_counts.txt",
                  sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
wholeCell.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
mat <- wholeCell.KO$results[,c(8:27,1,3)]
colnames(mat)[21] <- "gene.name"
grp.idx <- WTgrp_kmeans(control_mat = mat[,1:10])
res1 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
                        WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                        WT2.idx = grp.idx$WT.idx2, bin.size = 200, 
                        shift.size = 40)
res1$plot

### Using the DEGs list
degs.dat <- read.table("../dat/DEGs/KO-WT_whole-cell_RNA-seq.txt", sep = "\t", 
                       stringsAsFactors=F, header = TRUE, row.names = 1)
cat("Number of DEGs =",dim(degs.dat)[1],"\n\n")
degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"),
                       refseq, by = "gene.name")
Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
             log2FC = log2(1), comp.between = "")
Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
             log2FC = log2(1.2), comp.between = "")
degs.dat <- inner_join(x = wholeCell.KO$results[wholeCell.KO$results$gene %in% 
                                                degs.dat$gene.name,c(8:27,1,3)],
                       y = degs.dat[,c("gene.name", "logFC")], 
                       by = c("gene" = "gene.name"))
mat <- degs.dat[,c(1:21,23)]
colnames(mat)[21] <- "gene.name"

### logFC from Degs from edgeR (Boxer et al.) and 
### much closer to the Boxer et al. paper
colnames(mat)[22] <- "log2FoldChange" 
res2 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
                        WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
                        WT2.idx = grp.idx$WT.idx2, bin.size = 60, 
                        shift.size = 6, shrink_lfc = T)
res2$plot

### logFC from Degs from DESeq2
mat <- degs.dat[,c(1:22)]
colnames(mat)[21] <- "gene.name" 
res3 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
                        WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
                        WT2.idx = grp.idx$WT.idx2, bin.size = 60, 
                        shift.size = 6, shrink_lfc = T)
res3$plot

### overlap_degs_wrapper
# res2 <- overlap_degs_wrapper(degs.dat, count.dat = wholeCell.KO$counts, 
#                              refseq, WT1.idx = grp.idx$WT.idx1, 
#                              WT2.idx = grp.idx$WT.idx2, bin.size = 60, 
#                              shift.size = 6)
# res2$plot
#