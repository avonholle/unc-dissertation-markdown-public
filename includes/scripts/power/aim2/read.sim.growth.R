# read.sim.growth.R
# Read in the sim.growth.R results (run on killdevil)

require(splines)
require(nlme)
require(ggplot2)
library(gridExtra)
require(reshape2)
require(plyr)

# 1) read in data
# .............................
setwd("C:/Users/vonholle/Dropbox/unc.grad.school/applications/aha-2016/power")

#load("sim.growth.Rda") # has sg object (a list of 3 different model coefficients for 1,000 iterations).
load("sim.growth-scenario2.Rda") # has sg object (a list of 3 different model coefficients for 1,000 iterations).
class(sg)
dim(sg)

# convert from list to data frame
# See http://stackoverflow.com/questions/14376506/how-to-extract-elements-from-list-of-lists
sg.1 = do.call(rbind.data.frame, lapply(sg, "[[", 1)) #  for non-interaction analysis
summary(sg.1); dim(sg.1)

sg.2 = do.call(rbind.data.frame, lapply(sg, "[[", 2)) #  for just snp analysis
summary(sg.2)

sg.3 = do.call(rbind.data.frame, lapply(sg, "[[", 3)) #  for snp and interaction analysis
summary(sg.3)
# 
# sg.4 = do.call(rbind.data.frame, lapply(sg, "[[", 4)) #  for predictors of growth curves
# summary(sg.4)


# 2) Get power numbers
# .............................................................

# First, a refresher.
# Model: y_i = beta.0.2 + beta.1.2*size_i + beta.2.2*tempo_i + beta.3.2*velocity_i + e_i
# Hypothesis: beta.3.2 is positive and significantly associated with y_i, HDL-C on a continuous scale.
# ................................................................................

# 2a) Look for percent significant coefficients
# .............................................

# Power for estimates with larger variance
l1 = list(sg.1, sg.2, sg.3) # make list of all 3 data frames w/ different scenarios

sapply(l1, names) # get column names for each of the data frames -- note, they are all the same

sg.1$sim=1
sg.2$sim=2
sg.3$sim=3

sg.dat = rbind(sg.1, sg.2, sg.3) # append all data frames together
names(sg.dat)
head(sg.dat)

sg.dat[sg.dat$sim==2, ][1:10,] # look at labels for effects for sim 2
sg.dat[sg.dat$sim==3, ][1:20,] # look at labels for effects for sim 3

sg.dat$sig = ifelse(sg.dat[,4]<0.05, 1, 0) # make indicator variable for significant effect

get.power = function(x){
  alpha.0 = mean(x[substr(rownames(x),1,6)=="alpha0",]$sig) # alpha0
  beta.0  = mean(x[substr(rownames(x),1,5)=="beta0",]$sig) # beta0
  beta.1 = mean(x[substr(rownames(x),1,5)=="beta1",]$sig) # beta1

  snp.1.add = mean(x[substr(rownames(x),1,9)=="snp.1.add",]$sig) # additive effect meant for model w/o interaction
  snp.1.dom = mean(x[substr(rownames(x),1,9)=="snp.1.dom",]$sig) # dominant effect meant for model w/o interaction

  alpha0.snp.1.add = mean(x[substr(rownames(x),1,16)=="alpha0:snp.1.add",]$sig) # additive effect for model w/ interaction
  beta0.snp.1.add = mean(x[substr(rownames(x),1,15)=="snp.1.add:beta0",]$sig) # additive effect for model w/ interaction
  beta1.snp.1.add = mean(x[substr(rownames(x),1,15)=="snp.1.add:beta1",]$sig) # additive effect for model w/ interaction

  alpha0.snp.1.dom = mean(x[substr(rownames(x),1,16)=="alpha0:snp.1.dom",]$sig) # additive effect for model w/ interaction
  beta0.snp.1.dom = mean(x[substr(rownames(x),1,15)=="beta0:snp.1.dom",]$sig) # additive effect for model w/ interaction
  beta1.snp.1.dom = mean(x[substr(rownames(x),1,15)=="beta1:snp.1.dom",]$sig) # additive effect for model w/ interaction
  
  return(c(alpha.0, beta.0, beta.1, 
           snp.1.add, snp.1.dom, 
           alpha0.snp.1.add, beta0.snp.1.add, beta1.snp.1.add,
           alpha0.snp.1.dom, beta0.snp.1.dom, beta1.snp.1.dom))
}

power.l = by(sg.dat, sg.dat$sim, FUN=get.power) # apply a function to a data frame split by factors (sim in this case)
power.l

power.dat = do.call(rbind.data.frame, power.l) # make list into data frame
colnames(power.dat) = c("alpha0", "beta0", "beta1", "snp.1.add", "snp.1.dom", 
                        "alpha0.snp.1.add", "beta0.snp.1.add", "beta1.snp.1.add",
                        "alpha0.snp.1.dom", "beta0.snp.1.dom", "beta1.snp.1.dom")
power.dat$sim=rownames(power.dat)
power.dat

write.csv(power.dat, "power-sitar.csv")