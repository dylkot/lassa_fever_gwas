library(GMMAT)
library(SeqArray)
library(stringr)
library(tools)

args = commandArgs(trailingOnly=TRUE)
phenofn <- args[1] # Contains columns for all variables in formula 
bedbase <- args[2] # Path to plink bed files with extension
z <- str_sub(bedbase, 1, -5)
print(z)
relfn <- args[3] # Path to a relatedness matrix containing IIDs for rows and columns that match bedbase
formula <- as.formula(args[4]) # Formula like "PHENO ~ COVAR1 + COVAR2 + COVAR3"
nulloutfile <- args[5] # Output file for the R object containing the null model model0
outfile <- args[6] # Output file containing P-values for the score test for each variant in bedbase
ncores <- as.numeric(args[7]) # Number of cores to use for the score tests to run in parallel

pheno <- read.table(phenofn, header = TRUE, sep = "\t")

## Relatedness
rel <- as.matrix(read.table(relfn, sep = "\t", check.names = FALSE, header=TRUE, row.names=1))

## Fit the null model
print(formula)
print(head(pheno))
print(rel[1:5, 1:5])
model0 <- glmmkin(formula, data = pheno, kins = rel,
                  id = "IID", family = binomial(link = "logit"),
                  singular.ok = FALSE)

save(model0, file=nulloutfile)

## If bed file, we actually want to pass the basename
if (file_ext(bedbase) != "bed"){    
    print('ERROR')
}
bedbase <- str_sub(bedbase, 1, -5)

## Convert bed to gds
gdsfn <- paste(bedbase, ".gds", sep="")
SeqArray::seqBED2GDS(paste(bedbase, ".bed", sep=""),
                      paste(bedbase, ".fam", sep=""),
                      paste(bedbase, ".bim", sep=""),
                      gdsfn)

## Test each of the variants against the null model
glmm.score(model0, infile=gdsfn, outfile=outfile, ncores=ncores)