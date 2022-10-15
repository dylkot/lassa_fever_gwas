library(GMMAT)
library(tools)
library(stringr)
library(SeqArray)

args = commandArgs(trailingOnly=TRUE)
modelfn <- args[1] # rds file saving the null model object
genotypefn <- args[2] # path to a text file, gds file, or PLINK bed file that can be used for associations. Assumes fam/bim file present as well for bed
outfn <- args[3] # output pathe
ncores <- as.numeric(args[4]) # ncores

load(modelfn) ## populates model0

## Convert bed to gds
gdsfn <- paste(genotypefn, ".gds", sep="")
SeqArray::seqVCF2GDS(genotypefn, gdsfn, parallel=TRUE, fmt.import="DS")
glmm.score(model0, infile = gdsfn, outfile=outfn, ncores=ncores, is.dosage=TRUE)

