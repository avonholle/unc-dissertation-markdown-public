# sim.growth-scenario2.R
# This is the version of sim.growth.one.R with a function to replicate
# 1,000 times.
# In scenario 2, the MAF for sample SNP will be 0.19 (for rs7241918)
# and if following Hardy-Weinberg equilibrium then genotype freq will be 0.04/0.31/0.65 (after rounding to equal 1)

#  pks = c("ggplot2", "gridExtra", "reshape")
#  install.packages(pks)

require(splines)
require(nlme)
require(ggplot2)
library(gridExtra)
require(reshape2)

# simulate infant growth data for figure in R
# some code adapted from simcurve.an.20131120.R
# ...............................................

# SITAR model 
# y_it = alpha_i + h( [t-beta_i] / exp(-gamma_i) )
# h is the cubic natural spline

# let t = 0, 1, 2, 5, 6, 7, 9, and 12 (8 fixed time points)

# make function for model similar to what was in the program to obtain parameters.
fitnlme2 <- function(age, 
                     s1,s2,s3,s4,
                     salpha0, sbeta0, sbeta1) {
  myknots=c(1,2,9) # 3 knots at 2, 6 and 9 months
  mybounds=c(0,12) # bounds at 0 to 12 months for spline
  
  nsmatrix <- as.matrix( ns( (age-sbeta0)/exp(-sbeta1), knots=myknots, Boundary.knots=mybounds) )
  as.vector(  salpha0 + 
                t(   matrix(rep(1,4),ncol=4) %*%
                       t( cbind(  s1*nsmatrix[,1], 
                                  s2*nsmatrix[,2], 
                                  s3*nsmatrix[,3], 
                                  s4*nsmatrix[,4]))
                )
  )
}

# 1) Generate the random effects 
# Simple model adapted from Beath 2008, doi: 10.1002/sim.2718
# ######################################################

power.1 = function(){
  # see http://stackoverflow.com/questions/17033577/try-or-trycatch-with-boot-r
  tryCatch({ # tryCatch function will return a missing value for iterations that don't work
    
    # Take function from Beath and use it to generate data 
    # I got these numbers from the sample plots in sim.growth.data.R
    s1=0.5
    s2=0.5
    s3=0.8
    s4=0.4
    
    # splinecoefs <- as.matrix(cbind(s1,s2,s3,s4)); splinecoefs
    # salpha0=0.5
    # sbeta0=0.2
    # sbeta1=-0.6
    
    age.sim = c(0, 1, 2, 5, 6, 7, 9, 12) # age at 0, 1, ..., 12 months
    
    
    # 2) 
    # simulate alpha0, beta0, and beta1
    # ..................................
    ns = 500
    alpha.0 = rnorm(ns, 2, 0.5)
    beta.0 = rnorm(ns, 0, 0.5)
    beta.1 = rnorm(ns, 0, 0.5)
    
    df.1 = cbind.data.frame(alpha.0, beta.0, beta.1)
    df.1$id = as.numeric(seq(1,100,1))
    head(df.1)
    
    df.2 = cbind(data.frame(t(age.sim)), id = seq(1,100,1))
    dim(df.2)
    head(df.2)
    
    df.3 = merge(df.1, df.2, by="id")
    head(df.3)
    
    df.4 = df.3
    df.4[,c(5:12)] = apply(df.3[,c(5:12)], 2, function(x) {
      (x - beta.0) / exp(beta.1)
    })
    #  head(df.4)
    
    # 3)  Now generate data frame to simulate curve.
    # #######################################
    
    # function with sbeta0 or sbeta1, meant to indicate random effects
    y.base <- function (x, alpha, beta0, beta1)
    {fitnlme2(x, s1, s2, s3, s4,
              alpha, beta0, beta1)
    }
    
    # Get outcome via ns function
    sim.dat.1 = data.frame(t(apply(df.4, 1, function(x) {y.base(age.sim, x[2], x[3], x[4])} )))
    sim.dat.1$id = rownames(sim.dat.1)
    
    
    # 4) Convert simulated data frame with simulated random effects from wide to long format
    # ######################################################
    
    # first, convert the df from wide to long for plotting
    long.sim.dat = melt(sim.dat.1, id.vars="id")
    
    long.sim.dat = within(long.sim.dat, {
      time = as.numeric(substr(variable, 2,2))
    })
    
    
    # Hypothesis: Children with less favorable weight change profiles, such as steeper weight change slopes (higher velocity), 
    #             are at greater risk of lower levels of HDL-C than children with less extreme weight change profiles.
    # Now run model to get p-value for ???
    
    # Model: y_i = beta.0.2 + beta.1.2*size_i + beta.2.2*tempo_i + beta.3.2*velocity_i + e_i
    # Hypothesis: beta.3.2 is positive and significantly associated with y_i, HDL-C on a continuous scale.
    # ................................................................................
    
    # set betas
    beta.0.2 = 40 # HDL-C level average in this group is 40 with sd=10 in sample of 667
    beta.1.2 = 0.25 # assume that 0.25 is a conservative estimate of the standard deviation of the betas (random effects for growth)
    beta.2.2 = 0.25
    beta.3.2 = 0.25
    
    beta.4.2 = 1   # snp.1.add effect
    beta.5.2 = 0  # snp.1.dom effect
    
    beta.6.2 = 0.25  # snp.1.add interaction effect with alpha0
    beta.7.2 = 0.25  # snp.1.add interaction effect with beta0
    beta.8.2 = 0.25  # snp.1.add interaction effect with beta1
    
    beta.9.2 = 0   # snp.1.dom interaction effect with alpha0
    beta.10.2 = 0  # snp.1.dom interaction effect with beta0
    beta.11.2 = 0 # snp.1.dom interaction effect with beta1
    
    beta.2 = c(beta.0.2, beta.1.2, beta.2.2, beta.3.2,
               beta.4.2, beta.5.2, beta.6.2, beta.7.2,
               beta.8.2, beta.10.2, beta.11.2)
    e.i.2 = rnorm(ns, 0, 0.25)
    #  summary(e.i.2)
    
    
    # 5) run sitar model on the simulated data set from above, sim.dat.1.
    # with growth curve specifications based on visual inspection of plots from sim.growth.data.R
    # extract out the three growth parameters to put in simulated y above.
    # .............................................................
    
    # set initial params for nlme function
    # based on simulation values...
    inits1 = s1
    inits2 = s2
    inits3 = s3
    inits4 = s4
    inits0 = 2
    
    test.reg = nlme(value ~ fitnlme2(time,s1,s2,s3,s4,alpha0,beta0,beta1),
                    data=long.sim.dat,
                    fixed = s1+s2+s3+s4+alpha0 ~ 1,
                    random = alpha0 + beta0 + beta1 ~ 1 | id,
                    start = c(inits1,inits2,inits3,inits4,inits0))
    
    rf = ranef(test.reg)
    
    # 6)  Run final model with random effects from nlme as covariates and hdl-c as response variable
    #     in simple linear regression
    
    # generate allele frequencies
    # ...............................
    # Genotype AA = 1, Aa = 2 and aa = 3
    # Gamete A =1 and a = 0
    # Uniformly distribute AA, Aa and aa
    snp.1 <- sample(1:3, ns, replace=TRUE, prob=c(0.04, 0.31, 0.65)) # genotype frequencies to match hypothetical distn with sample from Willer 2012, rs7241918
    #    table(snp.1)
    snp.1.add = ifelse(snp.1==1, 1, 
                       ifelse(snp.1==2, 0,
                              ifelse(snp.1==3, -1, NA)))
    
    snp.1.dom = ifelse(snp.1==1, -0.5, 
                       ifelse(snp.1==2, 0.5,
                              ifelse(snp.1==3, -0.5, NA)))
    
    y.2 = beta.0.2 + beta.1.2*rf$alpha0 + beta.2.2*rf$beta0 + beta.3.2*rf$beta1 + e.i.2 # no interaction, no snp
    y.3 = beta.0.2 + beta.1.2*rf$alpha0 + beta.2.2*rf$beta0 + beta.3.2*rf$beta1 + 
      beta.4.2*snp.1.add + beta.5.2*snp.1.dom +
      e.i.2 # with snp
    y.4 = beta.0.2 + beta.1.2*rf$alpha0 + beta.2.2*rf$beta0 + beta.3.2*rf$beta1 + 
      beta.4.2*snp.1.add + beta.5.2*snp.1.dom +
      beta.6.2*rf$alpha0*snp.1.add + beta.7.2*rf$beta0*snp.1.add + beta.8.2*rf$beta1*snp.1.add +
      beta.9.2*rf$alpha0*snp.1.dom + beta.10.2*rf$beta0*snp.1.dom + beta.11.2*rf$beta1*snp.1.dom +
      e.i.2 # with snp and interaction
    df.5 = data.frame(y.2, y.3, y.4, alpha0=rf$alpha0, beta0=rf$beta0, beta1=rf$beta1) # data frame for regression models
    
    # regression models
    m.1 = glm(y.2 ~ alpha0 + beta0 + beta1, data=df.5) # regression with just the three growth parameters
    
    
    m.2 = glm(y.3 ~ alpha0 + beta0 + beta1 +
                snp.1.add + snp.1.dom , data=df.5) # snps
    
    #  summary(m.2)
    m.3 = glm(y.4 ~ alpha0*snp.1.add + beta0*snp.1.add + beta1*snp.1.add + 
                alpha0*snp.1.dom + beta0*snp.1.dom + beta1*snp.1.dom, data=df.5) # snps and interaction
    
    return(list(summary(m.1)$coefficients, 
                summary(m.2)$coefficients,
                summary(m.3)$coefficients, 
                summary(m.aim1)$coefficients)) # just return the parameters and their tests
  }, 
  error = function(err) {return()} # return NULL on error
  )
}

sg = replicate(2, power.1(), simplify=F) # NOTE: when this is uploaded to killdevil I go up to 1000 replicates.

#sg = lapply(1:2, function(x) {power.1()}) 
#sg

# convert from list to data frame
# see http://stackoverflow.com/questions/14376506/how-to-extract-elements-from-list-of-lists
# sg.1 = do.call(rbind.data.frame, lapply(sg, "[[", 1)) #  for non-interaction analysis
# sg.1
# 
# sg.2 = do.call(rbind.data.frame, lapply(sg, "[[", 2)) #  for just snp analysis
# sg.2
# 
# sg.3 = do.call(rbind.data.frame, lapply(sg, "[[", 3)) #  for snp and interaction analysis
# sg.3


setwd("C:/Users/vonholle/Dropbox/unc.grad.school/applications/aha-2016/power")
save(sg, file="sim.growth-scenario2.Rda")


# only look at power for aim 1 -- characterize predictors of growth trajectory
# ..................................................

test.func = function(){
  # see http://stackoverflow.com/questions/17033577/try-or-trycatch-with-boot-r
  tryCatch({ # tryCatch function will return a missing value for iterations that don't work
    ns=500
    age = rnorm(ns, 27, 7)
    alpha0.s = 0.1*age + rnorm(ns,0,1) 
    m.aim1 = glm(alpha0.s ~ age) 
    return(list(summary(m.aim1)$coefficients)) # just return the parameters and their tests
  },
  error = function(err) {return()} # return NULL on error
  )
}

just.aim1 = replicate(1000, test.func(), simplify=F) # NOTE: when this is uploaded to killdevil I go up to 1000 replicates.
class(just.aim1)

aim.1.info = do.call(rbind.data.frame, lapply(just.aim1, "[[", 1)) #  for predictors of growth curves
summary(aim.1.info)
head(aim.1.info)


aim.1.info$sig = ifelse(aim.1.info[,4]<0.05, 1, 0) # make indicator variable for significant effect
age.power = mean(aim.1.info[substr(rownames(aim.1.info),1,3)=="age",]$sig) # power for age coefficient
