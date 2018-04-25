# read-vcf-snps.R
# Read of SLCS vcf files and get subset of 
# candidate snps for my dissertation work.
# ---------------------------------------

library(VariantAnnotation)
library(data.table)

source("https://bioconductor.org/biocLite.R")
biocLite("snpStats")

# read in group of snps to subset
# from create-list-snps-for-vcf.R
# ---------------------------------------------------

load("vars.Rdata") # contains the out.1 and out.2 data frame with chr and hg pos #
str(out.2)
summary(out.2)

list.pos = split(out.2, f=out.2$Chr)
list.pos[[length(list.pos)-2]] # check that variants are separated by chr
names = names(list.pos); names


# Select subsets of data from vcf files
# ---------------------------------------
setwd("../../../Genetic_Data_Center/chile/Imputed_Data")


# Look at the info.gz file
# ----------------------------------------
info.df <- read.table(gzfile("chr1.info.gz"), header=T)
head(info.df)


# Look at vcf file structure for one file
# --------------------------------------

hdr <- scanVcfHeader("chr1.dose.vcf.gz"); hdr
info(hdr) 
geno(hdr)


# function to select snps from each chromosome file
# -------------------------------------------------

test = lapply(names, function(chr) df=list.pos[[chr]])
class(test)
test[[1]]

# loop through all relevant chromosomes and retrieve snp info from each person
dat.snp = lapply(names, function(chr) {
  
  # set up parameters for reading in and subsetting the the vcf file
    param = ScanVcfParam(
      info=c("AF", "MAF", "R2", "ER2"), # AF: Estimated Alternate Allele Frequency
#      MAF: Estimated Minor Allele Frequency
#      R2: Estimated Imputation Accuracy
#      ER2: Empirical (Leave-One-Out) R-square (available only for genotyped variants)
      geno=c("GT", "DS"), # GT: Genotype; DS: Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]
      #GP: Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 
      which=GRanges(paste0(chr), 
                    IRanges( start=list.pos[[chr]]$Pos_hg19,
                             end=list.pos[[chr]]$Pos_hg19)))

    vcf = readVcf(paste0("chr", chr, ".dose.vcf.gz"),
                         "hg19", 
                         param=param)
    
    res <- genotypeToSnpMatrix(vcf)
    snps = as(res$genotype, "numeric")
    info = info(vcf)
    
    return(list(chr, snps, res$map, info)) # first item in list is chr number, second item is the snps for each person, third element is info on the allele coding
    
    })

setwd("../../../CVDGeneNas/avonholle/ms-d1")
save(dat.snp, file="datsnp.Rdata")

head(dat.snp[[2]][[2]]) # test
