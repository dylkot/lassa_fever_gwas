suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("HIBAG"))
suppressPackageStartupMessages(library("stringr"))

parser <- ArgumentParser()
parser$add_argument("-b", "--bed", type="character",
    help="Path to genotype file")
parser$add_argument("--hlafile", type="character",
    help="Path to HLA file")
parser$add_argument("--gene", type="character",
    help="HLA gene to use")
parser$add_argument("-o", "--out", type="character",
    help="Path to output the constructed model")
parser$add_argument("-f", "--flankbp", type="integer", default=500000, 
    help="Number of flanking BPs to consider [default %(default)s]",
    metavar="number")
parser$add_argument("-n", "--nclassifier", type="integer", default=100, 
    help="Number of classifiers to construct [default %(default)s]",
    metavar="number")
parser$add_argument("-c", "--cores", type="integer", default=100, 
    help="Number of cores to use [default %(default)s]",
    metavar="number")
args <- parser$parse_args()



## Load HLA data
hladat <- read.table(args$hlafile, header=TRUE, stringsAsFactors=FALSE, sep="\t")
head(hladat)

## The input is a PLINK bed file
plinkbase <- str_sub(args$bed,0,-5)
geno <- hlaBED2Geno(bed.fn=paste(plinkbase, ".bed", sep=""),
                    fam.fn=paste(plinkbase, ".fam", sep=""),
                    bim.fn=paste(plinkbase, ".bim", sep=""))

print(Sys.time())

hlaid = args$gene
A1 <- paste(hlaid, "1", sep="_")
A2 <- paste(hlaid, "2", sep="_")

ind = (hladat[[A1]]!='-') & (hladat[[A2]]!='-')
print(sum(ind))
train.HLA <- hlaAllele(hladat[ind,1], H1=hladat[ind,A1], H2=hladat[ind,A2], locus=hlaid)
snpid <- hlaFlankingSNP(geno$snp.id, geno$snp.position, hlaid, args$flankbp)
train.geno <- hlaGenoSubset(geno, snp.sel=match(snpid, geno$snp.id))
set.seed(1000)
model <- hlaParallelAttrBagging(args$cores, train.HLA, train.geno, nclassifier=args$nclassifier,verbose=TRUE)
model.obj <- hlaModelToObj(model)
save(model.obj, file=args$out)
print(Sys.time())
