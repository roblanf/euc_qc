#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

inputf = args[1]
outputf = args[2]

library(reshape2)

d = read.delim(inputf, header=F)
d = d[,c(1,2,3)]
names(d) = c("sample1", "sample2", "distance")

dst = acast(d, sample1 ~ sample2)
dst = data.matrix(dst)
diag(dst) = NA # set diagonals to NA to avoid colour washout, HT: https://twitter.com/BEDecato/status/847646772453285890
dim = ncol(dst)

samples = rownames(dst)
splits = strsplit(samples, "_")
splitmat = matrix(unlist(splits),ncol=3,byrow=TRUE)
samples = as.character(as.data.frame(splitmat)[,2])

rownames(dst) = samples
colnames(dst) = samples

pdf(file=outputf)

image(1:dim, 1:dim, dst, axes = FALSE, xlab="", ylab="")
axis(1, 1:dim, rownames(dst), cex.axis = 0.5, las=3)
axis(2, 1:dim, rownames(dst), cex.axis = 0.5, las=1)
text(expand.grid(1:dim, 1:dim), sprintf("%.3f", dst), cex=0.4)

dev.off()