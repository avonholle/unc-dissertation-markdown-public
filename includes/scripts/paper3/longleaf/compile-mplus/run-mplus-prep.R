# run-mplus-prep.R
# -------------------------------------------------

# program to run rest of programs to prep data and run MplusAutomate


#install.packages("markdown", repos="http://cran.r-project.org", lib="/pine/scr/v/o/vonholle/Rlibs", dependencies=TRUE)

library("knitr")
require("markdown")

# 1) Get imputed genotype data (selected from the /proj/epi/CVDGeneNas/avonholle/ms-d1/read-vcf-snps.R file)
# --------------------------------------------------------------

# TO DO: need to go back and make sure all data is read in this folder
# see https://stackoverflow.com/questions/10646665/how-to-convert-r-markdown-to-html-i-e-what-does-knit-html-do-in-rstudio-0-9 

# NOTE: prior to running this program need to transfer over the 
# sig-snp.csv file from 
# ~\GitHub\unc-dissertation-markdown\includes\scripts\power\aim3\power-calcs-ind-assoc.Rmd
# and phenotype data
# phen.Rdata from 
# ~/GitHub/unc-dissertation-markdown-p2/includes/scripts/paper3/longleaf/compile-mplus
# made from get-phen-data.Rmd in ~/GitHub/unc-dissertation-markdown-p2/includes/scripts/paper3/longleaf/
# lastly, need the 
# IDsToMatch.csv file 
# that has cross walk between geno and pheno ids
# from ~\Dropbox\unc.grad.school\my-papers\ms-201608-1\programs\kure-analysis

knit("get-snp-data.Rmd") 
markdownToHTML('get-snp-data.md', 'get-snp-data.html') # creates html file


# 3) make Mplus .inp files
# -------------------------------------------------------------------
knit("m3-data-scripts.Rmd")
markdownToHTML("m3-data-scripts.md", 'm3-data-scripts.html')
