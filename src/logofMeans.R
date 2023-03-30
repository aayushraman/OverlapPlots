## log of means between two/three samples
logofMeans.between.A.B <- function(dat, A.samples, B.samples){
  dat$Mean.A <- apply(dat[,A.samples], 1, function(r) {(mean(r))})
  dat$Mean.B <- apply(dat[,B.samples], 1, function(r) {(mean(r))})
  dat$FC.crude <- apply(dat[,c("Mean.A", "Mean.B")], 1, function(r) {(r[2]/r[1])})
  dat$logFC.crude <- apply(dat[,c("Mean.A", "Mean.B")], 1, function(r) {log2((r[2]+1)/(r[1]+1))})
  dat <- dat[!is.na(dat$logFC.crude),]
  return(dat)
}

logofMeans.between.ABC <- function(dat, A.samples, B.samples, C.samples){
  dat$Mean.A <- apply(dat[,A.samples], 1, function(r) {(mean(r))})
  dat$Mean.B <- apply(dat[,B.samples], 1, function(r) {(mean(r))})
  dat$Mean.C <- apply(dat[,C.samples], 1, function(r) {(mean(r))})
  dat$FC.crude <- apply(dat[,c("Mean.A", "Mean.B", "Mean.C")], 1, function(r) {((r[2]-r[1])/(r[3]-r[1]))})
  dat$logFC.crude <- apply(dat[,c("Mean.A", "Mean.B", "Mean.C")], 1, function(r) {log2((r[2]-r[1]+1)/(r[3]-r[1]+1))})
  dat <- dat[!is.na(dat$logFC.crude),]
  return(dat)
}

## log2FC between genotypes
log2FCwithingenotypes <- function(dat){

  ## converting log2 values to intensity and then calculating log2FC
  log2FC.dat <- data.frame(row.names = rownames(dat))
  rownames(dat) <- dat$gene.name
  gene.length <- dat[,c("gene.length","gene.name")]
  dat <- dat[,!names(dat) %in% c("gene.name","gene.length")]
  for(j in 2:ncol(dat)){
    i = j-1
    if(j > i){
      sample1 <- colnames(dat[,c(i,j)])
      sample2 <- colnames(dat[,-c(i,j)])
      dat.mean <- logofMeans.between.A.B(dat = dat, A.samples = sample2, B.samples = sample1)
      log2FC.dat <- cbind(log2FC.dat, dat.mean$logFC.crude)
    }
  }
  colnames(log2FC.dat) <- c("comp.mat1", "comp.mat2", "comp.mat3")
  log2FC.dat <- data.frame(log2FC.dat, gene.length)
  return(log2FC.dat)
}

## log2FC between genotypes
# log2FCwithingenotypes <- function(dat, gene.length){
#   
#   ## converting log2 values to intensity and then calculating log2FC
#   log2FC.dat <- data.frame(row.names = rownames(dat))
#   rownames(dat) <- dat$gene.name
#   dat <- dat[,-1]
#   for(j in 2:ncol(dat)){
#     i = j-1
#     if(j > i){
#       dat.mean = data.frame((apply(dat, 1, function(r) 
#         log2((mean((as.numeric(as.matrix(r[c(i,j)]))))+1)/(mean((as.numeric(as.matrix(r[-c(i,j)]))))+1)))))
#       log2FC.dat <- cbind(log2FC.dat, dat.mean)
#     }     
#   }  
#   colnames(log2FC.dat) <- c("comp.mat1", "comp.mat2", "comp.mat3")
#   log2FC.dat <- data.frame(log2FC.dat, gene.length)
#   return(log2FC.dat)
# }

## merging the exprs and annot
# merge.dat.annot <- function(exprs.dat, annot.mat){
#   ## for RNA-Seq
#   if(missing(annot.mat)){
#     dat.agg <- aggregate(. ~ gene.name, data = exprs.dat %>% dplyr::select(matches("SEQC_|gene.name")), max)
#     gene.length <- aggregate(. ~ gene.name, data = unique(exprs.dat 
#                                                           %>% dplyr::select(matches("gene.name|gene.length"))), max)  
#     dat.agg <- inner_join(x = dat.agg, y = gene.length, by = "gene.name")
#   }
#   else{## for array
#       annot.mat <- annot.mat %>% dplyr::select(contains("gene.")) %>% rownames_to_column() %>% mutate(rowname=as.character(rowname))
#       exprs.dat <- data.frame(exprs.dat) %>% rownames_to_column() %>% mutate(rowname = as.character(rowname))
#       data.annot <- inner_join(x = exprs.dat, y = annot.mat, by="rowname")
#       dat.agg <- aggregate(. ~ gene.name, data = data.annot %>% dplyr::select(matches("GSM|gene.name")), mean)
#       #print(dim(dat.agg))
#       dat.agg <- inner_join(x = dat.agg, y = unique(data.annot %>% dplyr::select(matches("gene.name|gene.length"))), 
#                             by = "gene.name")
#       dups.gene <- dat.agg[duplicated(dat.agg$gene.name),"gene.name"]
#       for(i in 1:length(dups.gene)){
#         idx <- which(dat.agg$gene.name == dups.gene[i])
#         #cat(dups.gene[i],"\t",dat.agg[idx, "gene.length"],"\n")
#         max.gene.length <- max(dat.agg[idx, "gene.length"])
#         dat.agg$gene.length[idx] = max.gene.length
#       }
#       dat.agg <- unique(dat.agg)
#   }
#   dat.agg <- unique(dat.agg)
#   #print(dim(dat.agg))
#   return(dat.agg)
# }
# 
# ## illumina merge annot
# merge.dat.annot.illu <- function(exprs.dat, annot.mat){
#   annot.mat <- annot.mat %>% dplyr::select(matches("Probe|gene."))
#   exprs.dat <- data.frame(exprs.dat) %>% rownames_to_column()
#   data.annot <- inner_join(x = exprs.dat, y = annot.mat, by=c("rowname" = "Probe_Id"))
#   dat.agg <- aggregate(. ~ gene.name, data = data.annot %>% dplyr::select(matches("GSM|gene.name")), mean)
#   #print(dim(dat.agg))
#   dat.agg <- inner_join(x = dat.agg, y = unique(data.annot %>% dplyr::select(matches("gene.name|gene.length"))), 
#                         by = "gene.name")
#   dups.gene <- dat.agg[duplicated(dat.agg$gene.name),"gene.name"]
#   for(i in 1:length(dups.gene)){
#     idx <- which(dat.agg$gene.name == dups.gene[i])
#     #cat(dups.gene[i],"\t",dat.agg[idx, "gene.length"],"\n")
#     max.gene.length <- max(dat.agg[idx, "gene.length"])
#     dat.agg$gene.length[idx] = max.gene.length
#   }
#   dat.agg <- unique(dat.agg)
#   dat.agg <- unique(dat.agg)
#   #print(dim(dat.agg))
#   return(dat.agg)
# }