## DESeq Run
DESeqCalculation <- function(dat, genotypes, fc = 1.15){
  colData <- data.frame(samples = colnames(dat), genotypes = factor(genotypes), 
                        row.names = colnames(dat))
  dds <- DESeqDataSetFromMatrix(countData = dat, colData = colData, 
                                design = ~ genotypes)
  dds <- estimateSizeFactors(dds)
  dat <- counts(dds, normalized=TRUE)
  print(dim(dat))

  ## size factor vs read counts
  sizeFactorRatio <- sizeFactors(dds)
  readCounts <- colSums(counts(dds))
  dds.mat <- data.frame(x = sizeFactors(dds), y = colSums(counts(dds))/1e6)
  p1 <- ggplot(dds.mat, aes(x = x, y = y)) + geom_point(aes(colour = dds$genotypes), size = 5) +
            geom_smooth(method = "lm", se = TRUE, colour = "grey30") + xlab("Size Factor") + 
            ylab("Number of Aligned Reads (in million)") + theme_bw() +
            theme(axis.title = element_text(size = 16, face = "bold"),
                  axis.text.x = element_text(size = 16, face = "bold", color = "black"),
                  axis.text.y = element_text(size = 16, face = "bold", color = "black"),
                  plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"), legend.position = c(.9, .1),
                  legend.title = element_blank(), 
                  legend.text = element_text(size = 14, face = "bold", color = "black"))
  
  ## plots
  boxPlot(data = log2(dat+1), samples = factor(genotypes))
  rld.dds <- assay(rlog(dds, blind = FALSE))
  p2 <- PCAplot(data = rld.dds, genotypes = colData$genotypes, conditions = colData$genotypes)
  p2 <- p2 + ggrepel::geom_text_repel(aes(label = colData$samples), size = 6, fontface = "bold", 
                                  color = "black", box.padding = unit(0.35, "lines"), 
                             point.padding = unit(0.2, "lines")) + theme_bw() +
                  theme(axis.title = element_text(size = 22, face = "bold"),
                        axis.text.x = element_text(size = 22, face = "bold", color = "black"),
                        axis.text.y = element_text(size = 22, face = "bold", color = "black"),
                        legend.title = element_text(size = 16, face = "bold", color = "black"),
                        legend.text = element_text(size = 16, face = "bold", color = "black"),
                        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  p3 <- PCAplot(data = log2(dat+1), genotypes = colData$genotypes, conditions = colData$genotypes)
  p3 <- p3 + ggrepel::geom_text_repel(aes(label = colData$samples), size = 6, fontface = "bold", 
                                      color = "black", box.padding = unit(0.35, "lines"), 
                                      point.padding = unit(0.2, "lines")) + theme_bw() +
                  theme(axis.title = element_text(size = 22, face = "bold"),
                        axis.text.x = element_text(size = 22, face = "bold", color = "black"),
                        axis.text.y = element_text(size = 22, face = "bold", color = "black"),
                        legend.title = element_text(size = 16, face = "bold", color = "black"),
                        legend.text = element_text(size = 16, face = "bold", color = "black"),
                        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  p4 <- MDSplot(data = log2(dat+1), genotypes = colData$genotypes, conditions = colData$genotypes)
  p4 <- p4 + ggrepel::geom_text_repel(aes(label = colData$samples), size = 6, fontface = "bold", 
                                      color = "black", box.padding = unit(0.35, "lines"), 
                                      point.padding = unit(0.2, "lines")) + theme_bw() +
                  theme(axis.title = element_text(size = 22, face = "bold"),
                        axis.text.x = element_text(size = 22, face = "bold", color = "black"),
                        axis.text.y = element_text(size = 22, face = "bold", color = "black"),
                        legend.title = element_text(size = 16, face = "bold", color = "black"),
                        legend.text = element_text(size = 16, face = "bold", color = "black"),
                        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  
  ## Expression Test
  dds <- DESeq(dds, betaPrior = TRUE)
  if(sum(unique(genotypes) %in% "KO") > 0){
    res.dds <- results(dds,contrast = c("genotypes","KO","WT"))
  }
  else{
    res.dds <- results(dds,contrast = c("genotypes","MUT","WT"))
  }
  res.dds$norm.counts <- counts(dds, normalized=TRUE)
  message(sum(res.dds$padj < 0.05, na.rm = TRUE))
  results <- as.data.frame(res.dds)
  
  ## Sorted as per adj. P-value
  resSort <- res.dds[order(res.dds$padj),]
  topGenes <- resSort[1:20,]
  topGenes$genes <- rownames(topGenes)
  
  ## Histogram and MA Plot with top 10 genes
  par(mfrow=c(1,1))
  hist(res.dds$pvalue[res.dds$baseMean > 1], breaks=0:20/20, col="grey50", 
       border="white", main="Histogram of p-values with baseMean > 1")
  DESeq2::plotMA(resSort, main = "MA Plot")
  for(i in 1:nrow(topGenes)){
    with(topGenes[i, ],{
      points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
      text(baseMean, log2FoldChange, genes, pos=2, col="dodgerblue")})
  }
  
  ## Volcano Plot
  results <- results %>% tibble::rownames_to_column()
  colnames(results)[1] <- "gene"
  results <- results[which(!is.na(results$padj)),]
  results$Significant <- ifelse(results$log2FoldChange > log2(fc) & results$padj < 0.05, "Up",
                                ifelse(results$log2FoldChange < log2(1/fc) & results$padj < 0.05, "Down","Not Signif"))
  pval.sig <- max(results[which(results$padj < 0.05),"pvalue"])
  p5 <- ggplot(results, aes(x = log2FoldChange, y = -log10(pvalue))) + 
            geom_point(aes(color = Significant)) +  
            scale_color_manual(values = c("green", "grey", "red")) +  
            theme_bw(base_size = 16) + xlab("Log2 Fold Change") + 
            ylab("-Log10 P-value") +
            geom_hline(aes(yintercept = -log10(pval.sig)), color="dodgerblue", 
                       linetype="dashed") + 
            geom_vline(aes(xintercept = log2(fc)), color="dodgerblue", 
                       linetype="dashed") +
            geom_vline(aes(xintercept = log2(1/fc)), color="dodgerblue", 
                       linetype="dashed") +
            theme(axis.title = element_text(size = 22, face = "bold"),
                  axis.text.x = element_text(size = 22, face = "bold", 
                                             color = "black"),
                  axis.text.y = element_text(size = 22, face = "bold", 
                                             color = "black"),
                  legend.title = element_text(size = 16, face = "bold", 
                                              color = "black"),
                  legend.text = element_text(size = 16, face = "bold", 
                                             color = "black"),
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  
  ## normalized values for the heatmaps
  topVarGenes <- order(rowVars(rld.dds),decreasing=TRUE)
  mat <- rld.dds[topVarGenes[1:50], ]
  topVarGenes <- order(rowVars(rld.dds),decreasing=TRUE)
  annot <- data.frame(genotypes = factor(genotypes), row.names = colnames(dat))
  if(sum(levels(annot$genotypes) %in% c("WT", "MUT", "KO")) >= 2){
    annot$genotypes <- relevel(annot$genotypes, ref = "WT")
  }
  p1.top <- pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE, 
                     show_rownames = TRUE, show_colnames = TRUE, 
                     fontsize_row = 9, legend = TRUE, filename = NA, 
                     fontsize_col = 10, scale = "row", #fontface="bold", 
                     clustering_distance_rows = "correlation", 
                     clustering_distance_cols = "euclidean", 
                     annotation_col = annot)
  
  ## upregulated and downregulated genes
  ind.up <- which(results$log2FoldChange > log2(fc) & results$padj < 0.05)
  ind.down <- which(results$log2FoldChange < log2(1/fc) & results$padj < 0.05)
  up.reg <- results[ind.up, ]
  down.reg <- results[ind.down, ]
  
  ## upregulated heatmaps
  up.rld.dds <- rld.dds[rownames(rld.dds) %in% up.reg$gene,]
  up.topVarGenes <- order(rowVars(up.rld.dds),decreasing=TRUE)
  up.mat <- up.rld.dds[up.topVarGenes[1:50], ]
  p2.up <- pheatmap(up.mat, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE,
                    show_colnames = TRUE, fontsize_row = 9, legend = TRUE, filename = NA, 
                    fontsize_col = 10, scale = "row", #fontface="bold", 
                    clustering_distance_rows = "correlation", 
                    clustering_distance_cols = "euclidean", annotation_col = annot)

  ## downregulated heatmaps
  down.rld.dds <- rld.dds[rownames(rld.dds) %in% down.reg$gene,]
  down.topVarGenes <- order(rowVars(down.rld.dds),decreasing=TRUE)
  down.mat <- down.rld.dds[down.topVarGenes[1:50], ]
  p3.down <- pheatmap(down.mat, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE,
                      show_colnames = TRUE, fontsize_row = 9, legend = TRUE, filename = NA,
                      fontsize_col = 10, scale = "row", #fontface="bold", 
                      clustering_distance_rows = "correlation", 
                      clustering_distance_cols = "euclidean", annotation_col = annot)
  
  ## PCA plot on high variable genes
  p6 <- PCAplot(data = log2(dat[rowMeans(dat) > 30, ]+1), genotypes = colData$genotypes)
  # p6 <- p6 + ggrepel::geom_text_repel(aes(label = colData$samples), size = 6, fontface = "bold", 
  #                                     color = "black", box.padding = unit(0.35, "lines"), 
  #                                     point.padding = unit(0.2, "lines")) + theme_bw() +
  #               theme(axis.title = element_text(size = 22, face = "bold"),
  #                     axis.text.x = element_text(size = 22, face = "bold", color = "black"),
  #                     axis.text.y = element_text(size = 22, face = "bold", color = "black"),
  #                     legend.title = element_text(size = 16, face = "bold", color = "black"),
  #                     legend.text = element_text(size = 16, face = "bold", color = "black"),
  #                     plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  print(p6)
  
  ## DEGs plot
  results.list = list(results = results, up.reg = up.reg, down.reg = down.reg, counts = dat,
                      plots = list(sizefactorPlot = p1, PCAplot_rld = p2, PCAplot = p3, MDSplot = p4, 
                                   Volplot = p5, PCAplot_highvar = p6, heatmap1 = p1.top, heatmap2 = p2.up,
                                   heatmap3 = p3.down))
  return(results.list)
}