# create-list-snps-for-vcf.R
# --------------------------------

#library(biomaRt)
library(readxl)
library(plyr)

# read in list of snps from Misa
# from fine mapping paper by Zubair et al, 2016 (doi: 10.1093/hmg/ddw358)
# ---------------------------------------------------

# sample of one way to get location of rsid of interest
# -------------------------------------------------------
# snp_mart = useMart("ENSEMBL_MART_SNP", 
#                    host="grch37.ensembl.org",
#                    dataset="hsapiens_snp")
# 
# snp_ids = c("rs1042034", "rs884366")
# snp_attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end")
# snp_locations = getBM(attributes=snp_attributes, filters="snp_filter", 
#                       values=snp_ids, mart=snp_mart)
# snp_locations

# obtain rsid, chromosome number and hg19 position from excel file titled, 
# 'Copy of Full_list_lipids_known_loci_02102017-annotated.xlsx'

setwd("~/GitHub/unc-dissertation-markdown/includes/table-data")


# -------------------------------------------------
# 1) Get info to compare reference alleles for use in look-vcf-snps.Rmd
# -------------------------------------------------

# column values in Eur excel file, Copy of SOL_Lipids_Suppl_Tables_061617.xlsx,
# that match names of columns in HL excel file. Copy of Full_list_lipids_known_loci_02102017-annotated.xlsx

# X__9: Effect allele
# X__8: Pos_hg19
# X__11: Author
# X__1: rsid
# Beta: "Beta (SE)"

# NOTE: Wu had mmol/l
sheets = list('Supplementary Table S6',
                'Supplementary Table S7',
                'Supplementary Table S8')
substr(sheets[[1]],22,22) #check

head(read_excel("Copy of SOL_Lipids_Suppl_Tables_061617.xlsx", 
           sheet = 'Supplementary Table S6',
           skip=2))
# check colnames

# loop through each of the sheets for the lipids to get info
eur.dat = lapply(sheets, function(x) {

  # x='Supplementary Table S6' # for debugging purposes
  val.trait = substr(x,22,22)
  num.rows = ifelse(val.trait=='6', 128,
                    ifelse(val.trait=="7", 137,
                           ifelse(val.trait=="8", 79, NA)))
  
  vals = read_excel("Copy of SOL_Lipids_Suppl_Tables_061617.xlsx", 
                    sheet = x,
                    skip=2,
                    n_max=num.rows)
  
  col.keep =  c("X__9", "X__8", "X__11", "X__1", "X__7",
                "Beta", "TRUE/ FALSE")
  vals2 = vals[which(colnames(vals) %in% col.keep)]
  vals2 = vals2[col.keep] # make sure order is correct for renaming below

  colnames(vals2) = c("Effect allele", "Pos_hg19", "Author", "rsid",
                      "Chr", "Beta (SE)", "generalize")
  vals2$beta.units = ifelse(vals2$Author=="Wu", "mmol", "mg")
  
  vals2$Trait = ifelse(val.trait == '6', "LDL",
                       ifelse(val.trait == '7', "HDL",
                              ifelse(val.trait =='8', 'TG', NA)))
  
  return(vals2[vals2$generalize=="TRUE",]) # take out only the snps that generalize to HL population
})

# make data frame from lists
eur.dat.df = do.call(rbind.data.frame, eur.dat)
colnames(eur.dat.df)
dim(eur.dat.df) # 173 matches with number or rsids in each of those sheets above.
eur.dat.df$eur = 1 # make a marker for non-hl population snps
head(eur.dat.df)

# concatenate onto the HL-specific snps
setwd("~/GitHub/unc-dissertation-markdown/includes/table-data")
df1 = read_excel("Copy of Full_list_lipids_known_loci_02102017-annotated.xlsx", 
                                  sheet = 'GWASCAT+nonGWASCAT_lipids')
head(df1)
df1$eur=0
df.allele = rbind.fill(df1[df1$significance=="gw_sig",], 
                       eur.dat.df)
dim(df.allele) # adding hl snps now goes from 173 to 252
head(df.allele)
table(df.allele$beta.units)
table(df.allele$Author)
table(df.allele$eur)
nrow(df.allele)

# NOTE: keep in mind that I've selected out the snps but have not done any 
# transformations to make units equivalent. That part happens downstream in other
# programs like power-calcs-ind-assoc.Rmd


# -------------------------------------------------
# 2) Get rsid info for data extraction from .vcf files
# -------------------------------------------------

names(df.allele)
out.2 = unique(df.allele[c("rsid", "Pos_hg19", "Chr")]) # df.allele has 252 rows, but 
# snps are in df.allele more than once because of multiple traits and different authors, same rsid.
# so just get unique rsid
nrow(out.2) # 208 snps
nrow(unique(df.allele[c("rsid", "Pos_hg19", "Chr", "Trait", "Author")])) # 252

# Add files to output files
# -----------------------------------------------------
# NOTE: to do, need to double check if 362 snps sounds reasonable

setwd("~/GitHub/unc-dissertation-markdown-p2/includes/scripts/paper3/longleaf/")

save(df.allele, out.2, file="vars.Rdata") # has out.1 data frame with chr and hg19 position
# for variants for HL population.
# out2 has both variants from HL and from European population that generalizes
# from the European population to the HL population.
# will read this file in with variantannotation to subset the vcf files in 
# ~\GitHub\unc-dissertation-markdown-p2\includes\scripts\paper3\longleaf\read-vcf.R
# run on longleaf with vcf files in the /proj/epi/Genetic_Data_Center/chile/Imputed_Data/ 
# directory.

