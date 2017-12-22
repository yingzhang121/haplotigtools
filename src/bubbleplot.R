#!/usr/bin/env Rscript

library(R.utils)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2){
    stop("usage: Rscript --vanilla bubbleplot.r input_file length_of_primary_contig",call.=FALSE)
}
infile <- args[1]
scontig <- as.integer(args[2])

plot_h1 = function(st,end) {
  eq = function(x) {
    y=(2*x-st-end)/(end-st)
    return(0.9-y^2*0.4)
  }
  curve(eq, from=st, to=end, add=T, lwd=2)
}
plot_h2 = function(st,end) {
  eq = function(x) {
    y=(2*x-st-end)/(end-st)
    return(y^2*0.4)
  }
  curve(eq, from=st, to=end, add=T, lwd=2)
}

x <- read.table(infile, header=T)
x$V8 <- x$Primary_Start/scontig
x$V9 <- x$Primary_end/scontig
pcoords <- mergeIntervals(as.matrix(x[,8:9]))
pseg <- matrix(1L,nrow=nrow(pcoords)+1, ncol=2)
for (i in 1:nrow(pcoords)) {pseg[i,2]=pcoords[i,1]; pseg[i+1,1]=pcoords[i,2]}
pseg[1,1] <- 0
pseg <- pseg[ pseg[,1]!=pseg[,2], ]

pdf("bubbleplot.pdf")
plot(0:10, rep(0.5, 11), type="n", ann=F, axes=F, ylim=c(0,1), xlim=c(0, 10))
for (i in 1:nrow(pcoords)) {plot_h1(pcoords[i,1]*10, pcoords[i,2]*10)}
for (i in 1:nrow(pseg)){segments(pseg[i,1]*10, 0.5, pseg[i,2]*10, 0.5, lwd=2)}
for (i in 1:nrow(x)) {plot_h2(x[i,8]*10,x[i,9]*10)}
dev.off()
