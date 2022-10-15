# run meta-analysis

library(GMMAT)

args = commandArgs(trailingOnly=TRUE)
in1 <- args[1] 
in2 <- args[2]
outfn <- args[3]

# meta
glmm.score.meta(files = c(in1, in2), outfile = outfn, SNP = c('SNP', 'SNP'), A2 = c('REF', 'REF'), A1 = c('ALT', 'ALT'))