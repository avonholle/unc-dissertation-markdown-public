# create-list-snps-for-vcf.R
# --------------------------------

# NOTE: this can be run on local computer. Not using any SLS data.

#library(biomaRt)
library(readxl)
library(plyr)

# -------------------------------------------------
# 1) Get info to compare reference alleles for use in look-vcf-snps.Rmd
# ----------------------------------------------------------------------
# 1a) read in list of snps from Misa
# from fine mapping paper by Zubair et al, 2016 (doi: 10.1093/hmg/ddw358)
# These are European SNPs that generalize to HL ancestral groups via SOL
# Got these from Misa, see 9/13/2017 email for info on Copy of SOL_Lipids_Suppl_Tables_061617.xlsx
# file.

# NOTE: this file states that "References [and transformations and units], all strand converted to positive strand: "
# the vcf file also uses positive strand.
# ----------------------------------------------------------------------

#setwd("~/GitHub/unc-dissertation-markdown/includes/table-data")
setwd("~/../Dropbox/unc.grad.school/dissertation/data/genetics/")


# column values in Eur excel file, Copy of SOL_Lipids_Suppl_Tables_061617.xlsx,
# that match names of columns in HL excel file. Copy of Full_list_lipids_known_loci_02102017-annotated.xlsx
# Should match 
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


# From Zubair supplement: ddw358_Supp.docx (doi: 10.1093/hmg/ddw358)

# For each trait results from the GLGC the Betas were transformed from quintile 
# normalized by multiplying the beta and standard error estimates 
# by the standard deviation of the trait.  
# The standard deviation used for HDL is 13mg/dl, for LDL is 37 mg/dl, and for lnTG is 0.53 mg/dl.
# so transform the sol beta estimates this way?

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
  
  col.keep =  c("X__9", "X__10", "X__8", 
                "X__11", "X__1", "X__7",
                "Beta", "SE", 
                "Beta__1", "SE__1", "Effect allele frequency__1",
                "TRUE/ FALSE")
  vals2 = vals[which(colnames(vals) %in% col.keep)]
  vals2 = vals2[col.keep] # make sure order is correct for renaming below

  colnames(vals2) = c("Effect allele", "Other allele", "Pos_hg19", 
                      "Author", "rsid",
                      "Chr", "Beta (SE)", "SE", 
                      "beta.sol", "se.sol", 'eaf.sol',
                      "generalize")
  vals2$beta.units = ifelse(vals2$Author=="Wu", "mmol", "mg") # I would ignore this? all values were converted to mg?
  
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
length(unique(eur.dat.df$rsid)) # 154 unique variants of 173
eur.dat.df$eur=1

eur.dat.df[eur.dat.df$rsid=="rs10903129",]

# Taken from the Zubair 2016 paper:
# The standard deviation used for HDL is 13mg/dl, for LDL is 37 mg/dl, and for lnTG is 0.53 mg/dl.
# Not sure if this is correct
names(eur.dat.df)
eur.dat.df$beta.sol2 = with(eur.dat.df, ifelse(Trait=="LDL", beta.sol*37,
                                              ifelse(Trait=="HDL", beta.sol*13,
                                                     ifelse(Trait=="TG", beta.sol*0.53))))

# ===========================================================================
# 1b) Get HL-specific SNPs and concatenate onto European SNPs that generalize
# See 5/2/2017 email for info on Copy of Full_list_lipids_known_loci_02102017.xlsx file

# ===========================================================================

# concatenate onto the HL-specific snps from another list that is a catalog of all SNPs by ancestral groups
df1 = read_excel("Full_list_lipids_known_loci_02102017.xlsx", 
                                  sheet = 'GWASCAT+nonGWASCAT_lipids')
head(df1)
df1$eur=0
dim(df1) # 981 by 19
names(df1)

df1.hispanic = df1[df1$significance=="gw_sig" &
                     df1$Ancestry=="Hispanic",] # only take Hispanic and genome-wide sig.
dim(df1.hispanic) # 38 by 19
length(unique(df1.hispanic$rsid)) # 36 unique rsid

table(df1.hispanic$rsid %in% eur.dat.df$rsid) # 26 of 36 are unique to this data frame and not in generalized variants.
df1.hispanic$eur=0

# also get SNPs from Zubair et al 2016; doi: 10.1093/hmg/ddw358
# ==============================================================
df2 = read_excel("Supp-T8-Zubair-doi-201610-1093-hmg-ddw358.xlsx", 
                 sheet = 'GWASCAT+nonGWASCAT_lipids')
dim(df2) # check, 41x19
df2$beta.units = "mg"
df2$eur=0
names(df2)
names(df1.hispanic)

df3 = rbind.fill(df1.hispanic,df2)
head(df3)

dim(df3) # 79 by 20
length(unique(df3$rsid)) # 75 unique rsid


# ==========================================================
# 1c) combine both European and HL together
# ========================================================

df.allele = rbind.fill(df3, 
                       eur.dat.df)
dim(df.allele) # adding hl snps now goes from 173 to 252

head(df.allele)
table(df.allele$beta.units)
table(df.allele$Author)
table(df.allele$eur)
nrow(df.allele) # 252
length(unique(df.allele$rsid)) # 208 unique rsid from 252 rsid

# Will output these data for use in other programs such as 
# ~\GitHub\unc-dissertation-markdown-p2\includes\scripts\power\aim3\power-calcs-ind-assoc.Rmd

# -----------------------------------------------------
# 2) Get rsid info for data extraction from .vcf files
# -----------------------------------------------------

names(df.allele)
out.2 = unique(df.allele[c("rsid", "Pos_hg19", "Chr")]) 
nrow(out.2) # 208, df.allele has 252 rows, but 
# snps are in df.allele more than once because of multiple traits and different authors, same rsid.
# so just get unique rsid
nrow(out.2) # 208 snps
nrow(unique(df.allele[c("rsid", "Pos_hg19", "Chr", "Trait", "Author")])) # 252

out.2$Pos_hg19 = as.numeric(as.character(out.2$Pos_hg19))
summary(out.2)


df.allele[df.allele$rsid=="rs964184",] # example of rsid appearing more than once.

# Add files to output files
# -----------------------------------------------------


setwd("~/GitHub/unc-dissertation-markdown-p2/includes/scripts/paper3/longleaf/")
save(df.allele, out.2, eur.dat.df , file="vars.Rdata") # has out.1 data frame with chr and hg19 position
# for variants for HL population.


list.pos = split(out.2, f=out.2$Chr)
list.pos[[length(list.pos)-2]] # check that variants are separated by chr
names = names(list.pos)
test = lapply(names, function(chr) df=list.pos[[chr]])

test
# will read this file in with variantannotation to subset the vcf files in 
# ~\GitHub\unc-dissertation-markdown-p2\includes\scripts\paper3\longleaf\read-vcf-snps.R
# run on longleaf with vcf files in the /proj/epi/Genetic_Data_Center/chile/Imputed_Data/ 
# directory.

