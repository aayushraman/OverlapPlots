#############################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Updated Date: 18th July 2023
#
# Program is used for:
# 1. Overlap analyses for all the Boxer et al datasets
#############################################################################

## packages 
rm(list = ls())
setwd("~/Desktop/LongGene_Part-II/OverlapPlots/src/")
source("libraries.R")

## functions
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

finalRun_mCA <- function(count.file, degs.file, mCA, bin.size, shift.size,
                         title = "MeCP2 KO"){
    degs.dat <- read.table(file = degs.file, sep = "\t", stringsAsFactors = F, 
                           header = TRUE, row.names = 1)
    count.file <- count.file[,c(8:27,1,3)]
    colnames(count.file)[21] <- "gene.name"
    grp.idx <- WTgrp_kmeans(control_mat = count.file[,1:10])
    if(length(grp.idx$WT.idx1) != length(grp.idx$WT.idx2)){
        message("Cluster size is not equal, therefore run same size k-means 
                variation!")
        grp.idx <- WTgrp_kmeans_eqSize(control_mat = count.file[,1:10])
    }
    res <- overlap_degs_mCA_wrapper(degs.dat = degs.dat, count.dat = count.file,
                                    refseq = mCA, WT1.idx = grp.idx$WT.idx1, 
                                    WT2.idx = grp.idx$WT.idx2, 
                                    bin.size = bin.size, 
                                    shift.size = shift.size, 
                                    methyl.type = "mCA/CA")
    comb.plot <- ((res$overlapPlot + coord_cartesian(ylim = c(-0.25,0.25))) / 
                      res$mCAvGl_plot) | (res$plotBar / res$diffPlot)
    comb.plot <- comb.plot + plot_annotation(title = title, theme = theme(
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold")))
    return(list(res = res, combined = comb.plot))
}

## annotation dataset
load(file = "../dat-info/mm10_ncbi-refSeqGene_Dec2019.RData")
load(file = "../dat-info/mm10_ncbi-refSeqGene_Dec2019.RData")
genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))
mCA <- data.frame(readRDS("../dat/mCA/CAperGene_10wk_boxer.RDS"), stringsAsFactors = F)
mCA <- mCA[mCA$CA >= 5,]
mCA <- mCA[mCA$gene.length >= 4500,]
mCA$mCA.CA <- mCA$mCA/mCA$CA
mCA.sub <- mCA[,c(1,6,7)]

#######################################################
#
## overlap plots
#
#######################################################
 
## KO/WT whole cell dataset (ref: Fig 1D Boxer et al. Mol Cell 2020)
count.file <- "../dat/counts/GSE128178_10WT_10MeCP2_KO_whole_cell_RNAseq_exon_counts.txt"
degs.file <- "../dat/DEGs/KO-WT_whole-cell_RNA-seq.txt"
bin.size <- 60
shift.size <- 6
wholeCell.KO <- finalRun(count.file, genotypes, degs.file, bin.size, shift.size,
                         title = "KO/WT whole cell dataset")
mCA.wholeCell.KO <- finalRun_mCA(wholeCell.KO$res$results, degs.file, mCA.sub, 
                                 bin.size, shift.size, 
                                 title = "KO/WT whole cell dataset")
wholeCell.KO$overlapPlots$combined
mCA.wholeCell.KO$combined
rm(count.file, degs.file)

## KO/WT nuclear dataset  (ref: Fig 2E Boxer et al. Mol Cell 2020)
count.file <- "../dat/counts/GSE128178_10WT_10MeCP2_KO_nuclear_RNAseq_gene_body_counts.txt"
degs.file <- "../dat/DEGs/KO-WT_nuclear_RNA-seq.txt"
bin.size <- 40
shift.size <- 4
nuclear.KO <- finalRun(count.file, genotypes, degs.file, bin.size, shift.size,
                       title = "KO/WT nuclear dataset")
mCA.nuclear.KO <- finalRun_mCA(nuclear.KO$res$results, degs.file, mCA.sub, 
                               bin.size, shift.size, 
                               title = "KO/WT nuclear dataset")
nuclear.KO$overlapPlots$combined
mCA.nuclear.KO$combined
rm(count.file, degs.file)

## KO/WT chromatin dataset (ref: Fig 2F Boxer et al. Mol Cell 2020)
count.file <- "../dat/counts/GSE128178_10WT_10MeCP2_KO_chromatin_associated_RNAseq_gene_body_counts.txt"
degs.file <- "../dat/DEGs/KO-WT_chromatin_RNA-seq.txt"
bin.size <- 40
shift.size <- 4
chromatin.KO <- finalRun(count.file, genotypes, degs.file, bin.size, shift.size,
                         title = "KO/WT chromatin dataset")
mCA.chromatin.KO <- finalRun_mCA(chromatin.KO$res$results, degs.file, mCA.sub, 
                                 bin.size, shift.size, 
                                 title = "KO/WT chromatin dataset")
chromatin.KO$overlapPlots$combined
mCA.chromatin.KO$combined
rm(count.file, degs.file)

## R306C/WT whole cell dataset  (ref: Fig 3F Boxer et al. Mol Cell 2020)
count.file <- "../dat/counts/GSE128178_10WT_10R306C_whole_cell_RNAseq_exon_counts.txt"
degs.file <- "../dat/DEGs/R306C-WT_whole-cell_RNA-seq.txt"
bin.size <- 40
shift.size <- 4
wholeCell.R306C <- finalRun(count.file, genotypes, degs.file, bin.size, 
                            shift.size, title = "R306C/WT whole cell dataset")
mCA.wholeCell.R306C <- finalRun_mCA(wholeCell.R306C$res$results, degs.file, 
                                    mCA.sub, bin.size, shift.size,
                                    title = "R306C/WT whole cell dataset")
wholeCell.R306C$overlapPlots$combined
mCA.wholeCell.R306C$combined
rm(count.file, degs.file)

## R306C/WT nuclear dataset 
count.file <- "../dat/counts/GSE128178_10WT_10R306C_nuclear_RNAseq_gene_body_counts.txt"
degs.file <- "../dat/DEGs/R306C-WT_nuclear_RNA-seq.txt"
bin.size <- 40
shift.size <- 4
nuclear.R306C <- finalRun(count.file, genotypes, degs.file, bin.size,
                          shift.size, title = "R306C/WT nuclear dataset")
mCA.nuclear.R306C <- finalRun_mCA(nuclear.R306C$res$results, degs.file, mCA.sub, 
                                  bin.size, shift.size, 
                                  title = "R306C/WT nuclear dataset")
nuclear.R306C$overlapPlots$combined
mCA.nuclear.R306C$combined
rm(count.file, degs.file)

## R306C/WT chromatin dataset
count.file <- "../dat/counts/GSE128178_10WT_10R306C_chromatin_associated_RNAseq_gene_body_counts.txt"
degs.file <- "../dat/DEGs/R306C-WT_chromatin_RNA-seq.txt"
bin.size <- 40
shift.size <- 4
chromatin.R306C <- finalRun(count.file, genotypes, degs.file, bin.size, 
                            shift.size, title = "R306C/WT chromatin dataset")
mCA.chromatin.R306C <- finalRun_mCA(chromatin.R306C$res$results, degs.file, 
                                    mCA.sub, bin.size, shift.size, 
                                    title = "R306C/WT chromatin dataset")
chromatin.R306C$overlapPlots$combined
mCA.chromatin.R306C$combined
rm(count.file, degs.file)
