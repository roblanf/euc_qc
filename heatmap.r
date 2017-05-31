#!/usr/bin/env Rscript
#### plot distance matrix using Pretty Heatmap (pheatmap) ####

#inputf = "c:/devwork/NGS/Epauc/QC/mash/distances.tab"
#outputf = "c:/devwork/NGS/Epauc/QC/mash/heatmap.pdf"

args = commandArgs(trailingOnly=TRUE)

inputf = args[1]
outputf = args[2]

library(reshape2)
library(pheatmap)

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
pheatmap(dst, cluster_rows = F, cluster_cols = F, display_numbers = T, fontsize_number = 6, number_format = "%.3f")
dev.off()
