library(GMMAT)

args = commandArgs(trailingOnly=TRUE)
modelfn <- args[1] # rds file saving the null model object
genotypefn <- args[2] # path to a text file, gds file, or PLINK bed file that can be used for associations. Assumes fam/bim file present as well for bed
outfn <- args[3] # output pathe

load(modelfn) ## populates model0

## Convert bed to gds
glmm.score(model0, infile = genotypefn, outfile=outfn, is.dosage=FALSE, infile.ncol.skip=3, infile.nrow.skip = 1)
