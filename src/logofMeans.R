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
