## wrapper around overlap plots
overlap_wrapper <- function(dat, refseq, KO.idx, WT.idx, WT1.idx, WT2.idx, 
                            bin.size, shift.size, conf_int = 0.50, 
                            shrink_lfc = F){
    dat.annot <- inner_join(x = dat, y = refseq, by = "gene.name")
    cat("Control Group 1:",WT1.idx,"\n")
    cat("Control Group 2:",WT2.idx,"\n")
    gene.idx <- which(colnames(dat.annot) %in% "gene.name")
    log2FC.WT <- dat.annot[,c(WT.idx, gene.idx)]
    log2FC.WT$comp.mat <- apply(log2FC.WT[, WT.idx], 1, 
                                function(r){log2((mean(r[WT1.idx]) + 1) / 
                                                     (mean(r[WT2.idx]) +1))})
    if(shrink_lfc == FALSE){
        log2FC.KO <- logofMeans.between.A.B(dat = dat.annot, A.samples = WT.idx, 
                                            B.samples = KO.idx)
        log2FC.length <- inner_join(x = log2FC.WT[,c("gene.name","comp.mat")],
                                    y = log2FC.KO[,c("gene.name","logFC.crude",
                                                     "gene.length")], 
                                    by = "gene.name")
    }else{
        log2FC.length <- inner_join(x = log2FC.WT[,c("gene.name","comp.mat")],
                                    y = dat.annot[,c("gene.name",
                                                     "log2FoldChange", 
                                                     "gene.length")],
                                    by = "gene.name")
    }
    res <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)], 
                               comp.between1 = "(WT/WT)", 
                               comp.between2 = "(KO/WT)", 
                               bin.size, shift.size, conf_int)
    plot.margin <- unit(c(1,0.5,0.5,0.5), "cm")
    res.plot <- plot_grid(res$plot1 + coord_cartesian(ylim = c(-0.4,0.4)),
                          res$plot2 + coord_cartesian(ylim = c(0,50)), ncol = 1, 
                          align = 'v') + theme(plot.margin = plot.margin)
    return(list(res = res, plot = res.plot, log2FC.length = log2FC.length))
}

## original k-means
WTgrp_kmeans <- function(control_mat, centers = 2, iter.max = 1000){
    group <- kmeans(t(control_mat), centers, iter.max)$cluster
    idx1 <- which(group %in% 1)
    idx2 <- which(group %in% 2)
    return(list(WT.idx1 = idx1, WT.idx2 = idx2))
}

## k-means variation for equal cluster size 
## tutorial: https://elki-project.github.io/tutorial/same-size_k_means
## points are ordered by their distance to the closest center minus the distance
## to the farthest cluster. Each point is assigned to the best cluster in order.
WTgrp_kmeans_eqSize <- function(control_mat, centers = 2, iter.max = 1000){
    size <- ceiling(nrow(t(control_mat))/centers)
    group <- kmeans(t(control_mat), centers, iter.max)
    new_group <- rep(NA, nrow(t(control_mat)))
    new_centers <- lapply(1:centers, function(r){
        euc_dist <- control_mat - group$centers[r,]
        sqrt(apply(euc_dist, 2, function(x) sum(x^2)))})
    new_centers <- matrix(unlist(new_centers), ncol = centers)
    new_clust_size <- rep(0, centers)
    sample_ord <- order(apply(new_centers, 1, min) - apply(new_centers, 1, max))
    for(i in sample_ord){
        bestcl <- which.max(new_centers[i,])
        new_group[i] <- bestcl
        new_clust_size[bestcl] <- new_clust_size[bestcl] + 1
        if(new_clust_size[bestcl] >= size){
            new_centers[,bestcl] <- NA
        }
    }
    idx1 <- which(new_group %in% 1)
    idx2 <- which(new_group %in% 2)
    return(list(WT.idx1 = idx1, WT.idx2 = idx2))
}

overlap_degs_mCA_wrapper <- function(degs.dat, count.dat, refseq, WT1.idx, 
                                     WT2.idx, bin.size, shift.size, methyl.type, 
                                     degs = TRUE){
    if(degs == TRUE){
        degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
    }
    cat("Number of genes =",dim(degs.dat)[1],"\n\n")
    degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"), 
                           refseq, #%>% rownames_to_column(var = "gene.name"), 
                           by = "gene.name")
    ## overlap plot
    cat("Control Group 1:",WT1.idx,"\n")
    cat("Control Group 2:",WT2.idx,"\n")
    log2FC.WT <- data.frame(count.dat[,c(1:10, 21)], stringsAsFactors = F)
    log2FC.WT$comp.mat <- apply(count.dat[,c(WT1.idx, WT2.idx)], 1, 
                    function(r){log2((mean(r[1:length(WT1.idx)])+1) / 
                                         (mean(r[(length(WT1.idx)+1):10])+1))})
    log2FC.length <- inner_join(x = log2FC.WT[,c("gene.name","comp.mat")],
                                y = degs.dat[,c("gene.name","logFC", 
                                                "gene.length","mCA.CA")], 
                                by = "gene.name")
    message(dim(log2FC.length)[1])
    res <- overlay.mC(mat = log2FC.length[,c(2:5,1)], comp.between1 = "(WT/WT)",
                      comp.between2 = "(KO/WT)", bin.size = bin.size,
                      shift.size = shift.size, methyl.type = methyl.type)
    return(res = res)
}

## for DEGs 
# overlap_degs_wrapper <- function(degs.dat, count.dat, refseq,
#                                  WT1.idx, WT2.idx,
#                                  bin.size, shift.size){
#   degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
#   cat("Number of DEGs =",dim(degs.dat)[1],"\n\n")
#   degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"),
#                          refseq, by = "gene.name")
#   Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
#                log2FC = log2(1), comp.between = "")
#   Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
#                log2FC = log2(1.2), comp.between = "")
# 
#   ## overlap plot
#   cat("Control Group 1:",WT1.idx,"\n")
#   cat("Control Group 2:",WT2.idx,"\n")
#   log2FC.WT <- data.frame(count.dat[,c(1:10)], gene.name = rownames(count.dat),
#                           stringsAsFactors = F)
#   log2FC.WT$comp.mat <- apply(count.dat, 1,
#                               function(r){log2((mean(r[WT1.idx])+1) / 
#                                                    (mean(r[WT2.idx])+1))})
#   log2FC.length <- inner_join(x = log2FC.WT[,c("gene.name","comp.mat")],
#                               y = degs.dat[,c("gene.name","logFC", "gene.length")],
#                               by = "gene.name")
#   message(dim(log2FC.length)[1])
#   res <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)],
#                              comp.between1 = "(WT/WT)", comp.between2 = "(KO/WT)",
#                              bin.size = bin.size, shift.size = shift.size)
#   res.plot <- plot_grid(res$plot1 + coord_cartesian(ylim = c(-0.4,0.4)),
#                         res$plot2 + coord_cartesian(ylim = c(0,10)), 
#                         ncol = 1, align = 'v') +
#               theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
#   return(list(res = res, plot = res.plot, log2FC.length =log2FC.length))
# }
#
# overlap_degs_mCA_wrapper <- function(degs.dat, count.dat, refseq, 
#                                      WT1.idx, WT2.idx, bin.size, 
#                                      shift.size, methyl.type, degs = TRUE){
#   if(degs == TRUE){
#     degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
#   }
#   cat("Number of genes =",dim(degs.dat)[1],"\n\n")
#   degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"), 
#                          refseq, #%>% rownames_to_column(var = "gene.name"), 
#                          by = "gene.name")
#   ## overlap plot
#   cat("Control Group 1:",WT1.idx,"\n")
#   cat("Control Group 2:",WT2.idx,"\n")
#   log2FC.WT <- data.frame(count.dat[,c(1:10)], gene.name = rownames(count.dat), 
#                           stringsAsFactors = F)
#   log2FC.WT$comp.mat <- apply(count.dat[,c(WT1.idx,WT2.idx)], 1,
#                               function(r){log2((mean(r[1:length(WT1.idx)])+1)/(mean(r[(length(WT1.idx)+1):10])+1))})
#   log2FC.length <- inner_join(x = log2FC.WT[,c("gene.name","comp.mat")],
#                               y = degs.dat[,c("gene.name","logFC", 
#                                               "gene.length","mCA.CA")], 
#                               by = "gene.name")
#   message(dim(log2FC.length)[1])
#   res <- overlay.mC(mat = log2FC.length[,c(2:5,1)], comp.between1 = "(WT/WT)",
#                     comp.between2 = "(KO/WT)", bin.size = bin.size,
#                     shift.size = shift.size, methyl.type = methyl.type)
#   return(res = res)
# }