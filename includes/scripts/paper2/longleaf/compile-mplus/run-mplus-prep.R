# run-mplus-prep.R
# -------------------------------------------------

# program to run rest of programs to prep data and run MplusAutomate
# in longleaf use the run-files-data.sh shell script, by typing, "sbatch run-files-data.sh &"


#install.packages("markdown", repos="http://cran.r-project.org", lib="/pine/scr/v/o/vonholle/Rlibs", dependencies=TRUE)
library("knitr")
require("markdown")


# 1) make Mplus .inp files
# -------------------------------------------------------------------
knit("m2-data-scripts.Rmd") # create .md file
markdownToHTML("m2-data-scripts.md", 'm2-data-scripts.html') # knit md file to html w/o using pandoc
