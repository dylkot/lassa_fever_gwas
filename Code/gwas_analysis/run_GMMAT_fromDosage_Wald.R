library(GMMAT)
library(SeqArray)
library(stringr)
library(tools)

args = commandArgs(trailingOnly=TRUE)
phenofn <- args[1] # Contains columns for all variables in formula 
infile <- args[2] # Path to VCF
relfn <- args[3] # Path to a relatedness matrix containing IIDs for rows and columns that match bedbase
formula <- as.formula(args[4]) # Formula like "PHENO ~ COVAR1 + COVAR2 + COVAR3"
outfile <- args[5] # Output file containing P-values for the score test for each variant in bedbase
ncores <- as.numeric(args[6]) # Number of cores to use for the score tests to run in parallel

## Load phenotype file. Because PLINK wants phe to be 1,2 we subtract 1 from the phenotype value
pheno <- read.table(phenofn, header = TRUE, sep = "\t")

## Relatedness
rel <- as.matrix(read.table(relfn, sep = "\t", check.names = FALSE, header=TRUE, row.names=1))

gdsfn <- paste(infile, ".gds", sep="")
SeqArray::seqVCF2GDS(infile, gdsfn, parallel=TRUE, fmt.import="DS", scenario=c("imputation"))
snps = seqGetData(gdsfn, "annotation/id")
res <- glmm.wald(formula, data = pheno, kins = rel, id = "IID", random.slope = NULL, groups = NULL,
                 family = binomial(link = "logit"), infile = gdsfn,  snps=snps, is.dosage = TRUE, verbose = FALSE,
                 missing.method="impute2mean")
          
write.table(res, sep = '\t', file = outfile)