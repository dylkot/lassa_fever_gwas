suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("HIBAG"))
suppressPackageStartupMessages(library("stringr"))

parser <- ArgumentParser()
parser$add_argument("-b", "--bed", type="character",
    help="Path to genotype file")
parser$add_argument("-i", "--indexfile", type="character",
    help="Path to Rdata file containing pre-formed index")
parser$add_argument("--gene", type="character",
    help="HLA gene to use")
parser$add_argument("-o", "--out", type="character",
    help="Path to output the HLA predictions")
parser$add_argument("-p", "--proboutput", type="character",
    help="Path to output the HLA posterior probabilities")

args <- parser$parse_args()

## Load the input PLINK data
plinkbase <- str_sub(args$bed,0,-5)
geno <- hlaBED2Geno(bed.fn=paste(plinkbase, ".bed", sep=""),
                    fam.fn=paste(plinkbase, ".fam", sep=""),
                    bim.fn=paste(plinkbase, ".bim", sep=""))

hlaindex <- get(load(args$indexfile))
print(hlaindex)

model <- hlaModelFromObj(hlaindex)
hla <- hlaPredict(model, geno, allele.check=FALSE, same.strand=TRUE, type="response+prob")

write.table(hla$value, file = args$out, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(hla$postprob, file = args$proboutput, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")