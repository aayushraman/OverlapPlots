#############################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Updated Date: 12th July 2023
#
# Program is used for:
# 1. Overlap analyses for all the Boxer et al Dataset
#############################################################################

## packages 
rm(list = ls())
setwd("~/Desktop/LongGene_Part-II/OverlapPlots/src/")
source("libraries.R")

## function
finalRun <- function(count.file, genotypes, degs.file, bin.size, shift.size,
                     title = "MeCP2 KO"){
    ## counts file
    dat <- read.table(file = count.file, sep = "\t", stringsAsFactors = F, 
                      header = T, row.names = 1)
    
    ## running DESeq and all the plots are in dds object
    dds <- DESeqCalculation(dat = dat, genotypes = genotypes, fc=1.15)
    mat <- dds$results[,c(8:27,1,3)]
    colnames(mat)[21] <- "gene.name"
    grp.idx <- WTgrp_kmeans(control_mat = mat[,1:10])
    if(length(grp.idx$WT.idx1) != length(grp.idx$WT.idx2)){
        message("Cluster size is not equal, therefore run same size k-means 
                variation!")
        grp.idx <- WTgrp_kmeans_eqSize(control_mat = mat[,1:10])
    }
    
    ## overlap plot using all the genes
    op1 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
                           WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                           WT2.idx = grp.idx$WT.idx2, bin.size = 200, 
                           shift.size = 40)
    
    ## Using the DEGs list
    degs.dat <- read.table(file = degs.file, sep = "\t", stringsAsFactors = F, 
                           header = TRUE, row.names = 1)
    degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
    cat("Number of DEGs =",dim(degs.dat)[1],"\n\n")
    degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"), 
                           refseq, by = "gene.name")
    
    ## scatter plot (distribution between 2 different logFC)
    scater1 <- Plot.Scatter(log2FC = log2(1), comp.between = "(> log2FC(0))", 
                            dat = degs.dat[,c("gene.name","logFC", "FDR", 
                                              "gene.length")])
    scater2 <- Plot.Scatter(log2FC = log2(1.2), comp.between = "(> log2(1.2))",
                            dat = degs.dat[,c("gene.name","logFC", "FDR", 
                                              "gene.length")],)
    
    ## prep for overlap plots with DEGs genes
    degs.dat <- inner_join(y = degs.dat[,c("gene.name", "logFC")], 
                           x = dds$results[dds$results$gene %in% 
                                               degs.dat$gene.name,c(8:27,1,3)],
                           by = c("gene" = "gene.name"))
    mat <- degs.dat[,c(1:21,23)]
    colnames(mat)[21] <- "gene.name"
    colnames(mat)[22] <- "log2FoldChange"
    
    ## logFC from degs from edgeR (Boxer et al.)
    op2 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
                           WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
                           WT2.idx = grp.idx$WT.idx2, bin.size = bin.size, 
                           shift.size = shift.size, shrink_lfc = T)
    
    ## logFC from degs from DESeq2
    mat <- degs.dat[,c(1:22)]
    colnames(mat)[21] <- "gene.name" 
    op3 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
                           WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
                           WT2.idx = grp.idx$WT.idx2, bin.size = bin.size, 
                           shift.size = shift.size, shrink_lfc = T)
    ## plots
    p1 <- op1$plot + ggtitle("All genes") + 
        theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")) 
    p2  <- op2$plot + ggtitle("edgeR DEGs") + 
        theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
    p3 <- op3$plot + ggtitle("DESeq2 DEGs") + 
        theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))
    q1 <- (p1 | p2 | p3) + plot_annotation(title = title, theme = theme(
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold")))
    
    ## return
    return(list(res = dds, overlapPlots = list(combined = q1, all_genes = op1, 
            degs_edgeR = op2, degs_deseq2 = op3), scatterPlot1 = scater1, 
            scatterPlot2 = scater2))
}

## annotation dataset
load(file = "../dat-info/mm10_ncbi-refSeqGene_Dec2019.RData")
genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))

## KO/WT whole cell dataset (ref: Fig 1D Boxer et al. Mol Cell 2020)
count.file <- "../dat/counts/GSE128178_10WT_10MeCP2_KO_whole_cell_RNAseq_exon_counts.txt"
degs.file <- "../dat/DEGs/KO-WT_whole-cell_RNA-seq.txt"
bin.size <- 60
shift.size <- 6
wholeCell.KO <- finalRun(count.file, genotypes, degs.file, bin.size, shift.size,
                         title = "KO/WT whole cell dataset")

## KO/WT nuclear dataset  (ref: Fig 2E Boxer et al. Mol Cell 2020)
count.file <- "../dat/counts/GSE128178_10WT_10MeCP2_KO_nuclear_RNAseq_gene_body_counts.txt"
degs.file <- "../dat/DEGs/KO-WT_nuclear_RNA-seq.txt"
bin.size <- 40
shift.size <- 4
nuclear.KO <- finalRun(count.file, genotypes, degs.file, bin.size, shift.size,
                       title = "KO/WT nuclear dataset")

## KO/WT chromatin dataset (ref: Fig 2F Boxer et al. Mol Cell 2020)
count.file <- "../dat/counts/GSE128178_10WT_10MeCP2_KO_chromatin_associated_RNAseq_gene_body_counts.txt"
degs.file <- "../dat/DEGs/KO-WT_chromatin_RNA-seq.txt"
bin.size <- 40
shift.size <- 4
chromatin.KO <- finalRun(count.file, genotypes, degs.file, bin.size, shift.size,
                         title = "KO/WT chromatin dataset")

## R306C/WT whole cell dataset  (ref: Fig 3F Boxer et al. Mol Cell 2020)
count.file <- "../dat/counts/GSE128178_10WT_10R306C_whole_cell_RNAseq_exon_counts.txt"
degs.file <- "../dat/DEGs/R306C-WT_whole-cell_RNA-seq.txt"
bin.size <- 40
shift.size <- 4
wholeCell.R306C <- finalRun(count.file, genotypes, degs.file, bin.size, 
                            shift.size, title = "R306C/WT whole cell dataset")

## R306C/WT nuclear dataset 
count.file <- "../dat/counts/GSE128178_10WT_10R306C_nuclear_RNAseq_gene_body_counts.txt"
degs.file <- "../dat/DEGs/R306C-WT_nuclear_RNA-seq.txt"
bin.size <- 40
shift.size <- 4
nuclear.R306C <- finalRun(count.file, genotypes, degs.file, bin.size, 
                          shift.size, title = "R306C/WT nuclear dataset")

## R306C/WT chromatin dataset
count.file <- "../dat/counts/GSE128178_10WT_10R306C_chromatin_associated_RNAseq_gene_body_counts.txt"
degs.file <- "../dat/DEGs/R306C-WT_chromatin_RNA-seq.txt"
bin.size <- 40
shift.size <- 4
chromatin.R306C <- finalRun(count.file, genotypes, degs.file, bin.size, 
                            shift.size, title = "R306C/WT chromatin dataset")



## for checking if the results are identical --
# dat <- read.table("../dat/counts/GSE128178_10WT_10MeCP2_KO_whole_cell_RNAseq_exon_counts.txt",
#                   sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
# wholeCell.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
# mat <- wholeCell.KO$results[,c(8:27,1,3)]
# colnames(mat)[21] <- "gene.name"
# 
# ## same cluster sizes
# grp.idx <- WTgrp_kmeans(control_mat = mat[,1:10])
# res1 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
#                         WT2.idx = grp.idx$WT.idx2, bin.size = 200, 
#                         shift.size = 40)
# 
# ### Using the DEGs list
# degs.dat <- read.table("../dat/DEGs/KO-WT_whole-cell_RNA-seq.txt", sep = "\t", 
#                        stringsAsFactors=F, header = TRUE, row.names = 1)
# degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
# cat("Number of DEGs =",dim(degs.dat)[1],"\n\n")
# degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"),
#                        refseq, by = "gene.name")
# Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
#              log2FC = log2(1), comp.between = "")
# Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
#              log2FC = log2(1.2), comp.between = "")
# degs.dat <- inner_join(x = wholeCell.KO$results[wholeCell.KO$results$gene %in% 
#                                                 degs.dat$gene.name,c(8:27,1,3)],
#                        y = degs.dat[,c("gene.name", "logFC")], 
#                        by = c("gene" = "gene.name"))
# mat <- degs.dat[,c(1:21,23)]
# colnames(mat)[21] <- "gene.name"
# 
# ### logFC from Degs from edgeR (Boxer et al.) and 
# ### much closer to the Boxer et al. paper
# colnames(mat)[22] <- "log2FoldChange" 
# res2 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
#                         WT2.idx = grp.idx$WT.idx2, bin.size = 60, 
#                         shift.size = 6, shrink_lfc = T)
# 
# ### logFC from Degs from DESeq2
# mat <- degs.dat[,c(1:22)]
# colnames(mat)[21] <- "gene.name" 
# res3 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
#                         WT2.idx = grp.idx$WT.idx2, bin.size = 60, 
#                         shift.size = 6, shrink_lfc = T)
# 
# ## plots
# p1 <- res1$plot + ggtitle("All genes") + 
#     theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")) 
# p2  <- res2$plot + ggtitle("edgeR DEGs") + 
#     theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
# p3 <- res3$plot + ggtitle("DESeq2 DEGs") + 
#         theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))
# q1 <- (p1 | p2 | p3) + plot_annotation(title = "KO/WT whole cell dataset", 
#                                        theme = theme(plot.title = element_text(
#                                         size = 18, hjust = 0.5, face = "bold")))
# rm(dat, wholeCell.KO, mat, grp.idx, degs.dat)
# 
# #######################################################################
# #
# ## KO/WT nuclear dataset  (ref: Fig 2E Boxer et al. Mol Cell 2020)
# #
# #######################################################################
# 
# dat <- read.table("../dat/counts/GSE128178_10WT_10MeCP2_KO_nuclear_RNAseq_gene_body_counts.txt",
#                   sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
# nuclear.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
# mat <- nuclear.KO$results[,c(8:27,1,3)]
# colnames(mat)[21] <- "gene.name"
# grp.idx <- WTgrp_kmeans(control_mat = mat[,1:10])
# if(length(grp.idx$WT.idx1) != length(grp.idx$WT.idx2)){
#     message("Cluster size is not equal, therefore run same size k-means variation!")
#     grp.idx <- WTgrp_kmeans_eqSize(control_mat = mat[,1:10])
# }
# res4 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
#                         WT2.idx = grp.idx$WT.idx2, bin.size = 200, 
#                         shift.size = 40)
# 
# ### Using the DEGs list
# degs.dat <- read.table("../dat/DEGs/KO-WT_nuclear_RNA-seq.txt", sep = "\t", 
#                        stringsAsFactors=F, header = TRUE, row.names = 1)
# degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
# cat("Number of DEGs =",dim(degs.dat)[1],"\n\n")
# degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"),
#                        refseq, by = "gene.name")
# Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
#              log2FC = log2(1), comp.between = "")
# Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
#              log2FC = log2(1.2), comp.between = "")
# degs.dat <- inner_join(x = nuclear.KO$results[nuclear.KO$results$gene %in% 
#                                                     degs.dat$gene.name,c(8:27,1,3)],
#                        y = degs.dat[,c("gene.name", "logFC")], 
#                        by = c("gene" = "gene.name"))
# mat <- degs.dat[,c(1:21,23)]
# colnames(mat)[21] <- "gene.name"
# 
# ### logFC from Degs from edgeR (Boxer et al.) and 
# ### much closer to the Boxer et al. paper
# colnames(mat)[22] <- "log2FoldChange" 
# res5 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
#                         WT2.idx = grp.idx$WT.idx2, bin.size = 40, 
#                         shift.size = 4, shrink_lfc = T)
# 
# ### logFC from Degs from DESeq2
# mat <- degs.dat[,c(1:22)]
# colnames(mat)[21] <- "gene.name" 
# res6 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
#                         WT2.idx = grp.idx$WT.idx2, bin.size = 40, 
#                         shift.size = 4, shrink_lfc = T)
# 
# ## plots
# p4 <- res4$plot + ggtitle("All genes") + 
#     theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")) 
# p5  <- res5$plot + ggtitle("edgeR DEGs") + 
#     theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
# p6 <- res6$plot + ggtitle("DESeq2 DEGs") + 
#     theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))
# q2 <- (p4 | p5 | p6) + plot_annotation(title = "KO/WT nuclear dataset", 
#                                        theme = theme(plot.title = element_text(
#                                         size = 18, hjust = 0.5, face = "bold")))
# rm(dat, nuclear.KO, mat, grp.idx, degs.dat)
# 
# #######################################################################
# #
# ## KO/WT chromatin dataset (ref: Fig 2F Boxer et al. Mol Cell 2020)
# #
# #######################################################################
# 
# dat <- read.table("../dat/counts/GSE128178_10WT_10MeCP2_KO_chromatin_associated_RNAseq_gene_body_counts.txt",
#                   sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
# chromatin.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
# mat <- chromatin.KO$results[,c(8:27,1,3)]
# colnames(mat)[21] <- "gene.name"
# grp.idx <- WTgrp_kmeans(control_mat = mat[,1:10])
# if(length(grp.idx$WT.idx1) != length(grp.idx$WT.idx2)){
#     message("Cluster size is not equal, therefore run same size k-means variation!")
#     grp.idx <- WTgrp_kmeans_eqSize(control_mat = mat[,1:10])
# }
# res7 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
#                         WT2.idx = grp.idx$WT.idx2, bin.size = 200, 
#                         shift.size = 40)
# 
# ### Using the DEGs list
# degs.dat <- read.table("../dat/DEGs/KO-WT_chromatin_RNA-seq.txt", sep = "\t", 
#                        stringsAsFactors=F, header = TRUE, row.names = 1)
# degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
# cat("Number of DEGs =",dim(degs.dat)[1],"\n\n")
# degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"),
#                        refseq, by = "gene.name")
# Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
#              log2FC = log2(1), comp.between = "")
# Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
#              log2FC = log2(1.2), comp.between = "")
# degs.dat <- inner_join(x = chromatin.KO$results[chromatin.KO$results$gene %in% 
#                                                   degs.dat$gene.name,c(8:27,1,3)],
#                        y = degs.dat[,c("gene.name", "logFC")], 
#                        by = c("gene" = "gene.name"))
# mat <- degs.dat[,c(1:21,23)]
# colnames(mat)[21] <- "gene.name"
# 
# ### logFC from Degs from edgeR (Boxer et al.) and 
# ### much closer to the Boxer et al. paper
# colnames(mat)[22] <- "log2FoldChange" 
# res8 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
#                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
#                         WT2.idx = grp.idx$WT.idx2, bin.size = 40, 
#                         shift.size = 4, shrink_lfc = T)
# 
# ### logFC from Degs from DESeq2
# mat <- degs.dat[,c(1:22)]
# colnames(mat)[21] <- "gene.name" 
# res9 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
#                         WT2.idx = grp.idx$WT.idx2, bin.size = 40, 
#                         shift.size = 4, shrink_lfc = T)
# 
# ## plots
# p7 <- res7$plot + ggtitle("All genes") + 
#     theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")) 
# p8  <- res8$plot + ggtitle("edgeR DEGs") + 
#     theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
# p9 <- res9$plot + ggtitle("DESeq2 DEGs") + 
#     theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))
# q3 <- (p7 | p8 | p9) + plot_annotation(title = "KO/WT chromatin dataset", 
#                                        theme = theme(plot.title = element_text(
#                                         size = 18, hjust = 0.5, face = "bold")))
# rm(dat, chromatin.KO, mat, grp.idx, degs.dat)
# 
# ##########################################################################
# #
# ## R306C/WT whole cell dataset  (ref: Fig 3F Boxer et al. Mol Cell 2020)
# #
# ##########################################################################
# 
# dat <- read.table("../dat/counts/GSE128178_10WT_10R306C_whole_cell_RNAseq_exon_counts.txt",
#                   sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
# R306C.wholeCell.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
# mat <- R306C.wholeCell.KO$results[,c(8:27,1,3)]
# colnames(mat)[21] <- "gene.name"
# 
# ## same cluster sizes
# grp.idx <- WTgrp_kmeans(control_mat = mat[,1:10])
# if(length(grp.idx$WT.idx1) != length(grp.idx$WT.idx2)){
#     message("Cluster size is not equal, therefore run same size k-means variation!")
#     grp.idx <- WTgrp_kmeans_eqSize(control_mat = mat[,1:10])
# }
# res10 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                          WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
#                          WT2.idx = grp.idx$WT.idx2, bin.size = 200, 
#                          shift.size = 40)
# 
# ### Using the DEGs list
# degs.dat <- read.table("../dat/DEGs/R306C-WT_whole-cell_RNA-seq.txt", sep = "\t", 
#                        stringsAsFactors=F, header = TRUE, row.names = 1)
# degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
# cat("Number of DEGs =",dim(degs.dat)[1],"\n\n")
# degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"),
#                        refseq, by = "gene.name")
# Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
#              log2FC = log2(1), comp.between = "")
# Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
#              log2FC = log2(1.2), comp.between = "")
# degs.dat <- inner_join(x = R306C.wholeCell.KO$results[R306C.wholeCell.KO$results$gene %in% 
#                                                     degs.dat$gene.name,c(8:27,1,3)],
#                        y = degs.dat[,c("gene.name", "logFC")], 
#                        by = c("gene" = "gene.name"))
# mat <- degs.dat[,c(1:21,23)]
# colnames(mat)[21] <- "gene.name"
# 
# ### logFC from Degs from edgeR (Boxer et al.) and 
# ### much closer to the Boxer et al. paper
# colnames(mat)[22] <- "log2FoldChange" 
# res11 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
#                         WT2.idx = grp.idx$WT.idx2, bin.size = 40, 
#                         shift.size = 4, shrink_lfc = T)
# 
# ### logFC from Degs from DESeq2
# mat <- degs.dat[,c(1:22)]
# colnames(mat)[21] <- "gene.name" 
# res12 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
#                         WT2.idx = grp.idx$WT.idx2, bin.size = 40, 
#                         shift.size = 4, shrink_lfc = T)
# 
# ## plots
# p10 <- res10$plot + ggtitle("All genes") + 
#     theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")) 
# p11  <- res11$plot + ggtitle("edgeR DEGs") + 
#     theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
# p12 <- res12$plot + ggtitle("DESeq2 DEGs") + 
#     theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))
# q4 <- (p10 | p11 | p12) + plot_annotation(title = "R306C/WT whole cell dataset", 
#                                         theme = theme(plot.title = element_text(
#                                         size = 18, hjust = 0.5, face = "bold")))
# rm(dat, R306C.wholeCell.KO, mat, grp.idx, degs.dat)
# 
# #######################################
# #
# ## R306C/WT nuclear dataset
# #
# #######################################
# 
# dat <- read.table("../dat/counts/GSE128178_10WT_10R306C_nuclear_RNAseq_gene_body_counts.txt",
#                   sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
# R306C.wholeCell.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
# mat <- R306C.wholeCell.KO$results[,c(8:27,1,3)]
# colnames(mat)[21] <- "gene.name"
# 
# ## same cluster sizes
# grp.idx <- WTgrp_kmeans(control_mat = mat[,1:10])
# if(length(grp.idx$WT.idx1) != length(grp.idx$WT.idx2)){
#     message("Cluster size is not equal, therefore run same size k-means variation!")
#     grp.idx <- WTgrp_kmeans_eqSize(control_mat = mat[,1:10])
# }
# res13 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                          WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
#                          WT2.idx = grp.idx$WT.idx2, bin.size = 200, 
#                          shift.size = 40)
# 
# ### Using the DEGs list
# degs.dat <- read.table("../dat/DEGs/R306C-WT_nuclear_RNA-seq.txt", sep = "\t", 
#                        stringsAsFactors=F, header = TRUE, row.names = 1)
# degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
# cat("Number of DEGs =",dim(degs.dat)[1],"\n\n")
# degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"),
#                        refseq, by = "gene.name")
# Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
#              log2FC = log2(1), comp.between = "")
# Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
#              log2FC = log2(1.2), comp.between = "")
# degs.dat <- inner_join(x = R306C.wholeCell.KO$results[R306C.wholeCell.KO$results$gene %in% 
#                                                           degs.dat$gene.name,c(8:27,1,3)],
#                        y = degs.dat[,c("gene.name", "logFC")], 
#                        by = c("gene" = "gene.name"))
# mat <- degs.dat[,c(1:21,23)]
# colnames(mat)[21] <- "gene.name"
# 
# ### logFC from Degs from edgeR (Boxer et al.) and 
# ### much closer to the Boxer et al. paper
# colnames(mat)[22] <- "log2FoldChange" 
# res14 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                          WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
#                          WT2.idx = grp.idx$WT.idx2, bin.size = 40, 
#                          shift.size = 4, shrink_lfc = T)
# 
# ### logFC from Degs from DESeq2
# mat <- degs.dat[,c(1:22)]
# colnames(mat)[21] <- "gene.name" 
# res15 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                          WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
#                          WT2.idx = grp.idx$WT.idx2, bin.size = 40, 
#                          shift.size = 4, shrink_lfc = T)
# 
# ## plots
# p13 <- res13$plot + ggtitle("All genes") + 
#     theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")) 
# p14  <- res14$plot + ggtitle("edgeR DEGs") + 
#     theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
# p15 <- res15$plot + ggtitle("DESeq2 DEGs") + 
#     theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))
# q5 <- (p13 | p14 | p15) + plot_annotation(title = "R306C/WT nuclear dataset", 
#                                     theme = theme(plot.title = element_text(
#                                         size = 18, hjust = 0.5, face = "bold")))
# rm(dat, R306C.wholeCell.KO, mat, grp.idx, degs.dat)
# 
# #######################################
# #
# ## R306C/WT chromatin dataset
# #
# #######################################
# 
# dat <- read.table("../dat/counts/GSE128178_10WT_10R306C_chromatin_associated_RNAseq_gene_body_counts.txt",
#                   sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
# R306C.wholeCell.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
# mat <- R306C.wholeCell.KO$results[,c(8:27,1,3)]
# colnames(mat)[21] <- "gene.name"
# 
# ## same cluster sizes
# grp.idx <- WTgrp_kmeans(control_mat = mat[,1:10])
# if(length(grp.idx$WT.idx1) != length(grp.idx$WT.idx2)){
#     message("Cluster size is not equal, therefore run same size k-means variation!")
#     grp.idx <- WTgrp_kmeans_eqSize(control_mat = mat[,1:10])
# }
# res16 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                          WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
#                          WT2.idx = grp.idx$WT.idx2, bin.size = 200, 
#                          shift.size = 40)
# 
# ### Using the DEGs list
# degs.dat <- read.table("../dat/DEGs/R306C-WT_chromatin_RNA-seq.txt", sep = "\t", 
#                        stringsAsFactors=F, header = TRUE, row.names = 1)
# degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
# cat("Number of DEGs =",dim(degs.dat)[1],"\n\n")
# degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"),
#                        refseq, by = "gene.name")
# Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
#              log2FC = log2(1), comp.between = "")
# Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
#              log2FC = log2(1.2), comp.between = "")
# degs.dat <- inner_join(x = R306C.wholeCell.KO$results[R306C.wholeCell.KO$results$gene %in% 
#                                                           degs.dat$gene.name,c(8:27,1,3)],
#                        y = degs.dat[,c("gene.name", "logFC")], 
#                        by = c("gene" = "gene.name"))
# mat <- degs.dat[,c(1:21,23)]
# colnames(mat)[21] <- "gene.name"
# 
# ### logFC from Degs from edgeR (Boxer et al.) and 
# ### much closer to the Boxer et al. paper
# colnames(mat)[22] <- "log2FoldChange" 
# res17 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                          WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
#                          WT2.idx = grp.idx$WT.idx2, bin.size = 40, 
#                          shift.size = 4, shrink_lfc = T)
# 
# ### logFC from Degs from DESeq2
# mat <- degs.dat[,c(1:22)]
# colnames(mat)[21] <- "gene.name" 
# res18 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20), 
#                          WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1, 
#                          WT2.idx = grp.idx$WT.idx2, bin.size = 40, 
#                          shift.size = 4, shrink_lfc = T)
# 
# ## plots
# p16 <- res13$plot + ggtitle("All genes") + 
#     theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")) 
# p17  <- res14$plot + ggtitle("edgeR DEGs") + 
#     theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
# p18 <- res15$plot + ggtitle("DESeq2 DEGs") + 
#     theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))
# q6 <- (p16 | p17 | p18) + plot_annotation(title = "R306C/WT chromatin dataset", 
#                                           theme = theme(plot.title = element_text(
#                                               size = 18, hjust = 0.5, face = "bold")))
# rm(dat, R306C.wholeCell.KO, mat, grp.idx, degs.dat)

### overlap_degs_wrapper
# res2 <- overlap_degs_wrapper(degs.dat, count.dat = wholeCell.KO$counts, 
#                              refseq, WT1.idx = grp.idx$WT.idx1, 
#                              WT2.idx = grp.idx$WT.idx2, bin.size = 60, 
#                              shift.size = 6)
# res2$plot
#