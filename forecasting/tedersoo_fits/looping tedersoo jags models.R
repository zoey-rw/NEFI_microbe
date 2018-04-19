#fitting tedersoo models in JAGS
#Note - this could be run in parallel, each dependent variable model is independent of others.
rm(list=ls())
library(runjags)
library(data.table)
library(MASS)
source('NEFI_functions/write_jags_beta.r')
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/jd_format.r')

#output path - where to save 25 model summaries and data sets.
output.path <- '/fs/data3/caverill/NEFI_microbial/prior_data/ted_prior_JAGS_models.rds'

#model path - where the jags code will be saved and loaded from.
jags.model.path <- "jags_beta_test.txt"

#load tedersoo prior data.
d <- data.table(readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/ted_all_prior_data.rds'))

#column position where dependent variables start.
start.position <- grep("Ectomycorrhizal", colnames(d))

#function to transform [0,1] to (0,1) a la Cribari-Neto & Zeileis 2010.
d[,start.position:ncol(d)] <- lapply(d[,start.position:ncol(d)],crib_fun)

#dependent variables to loop over
dependent.variables <- d[,start.position:ncol(d)]
#dependent.variables <- dependent.variables[,1:2] #subset for testing loop. 

#Specify fixed effects.
fixed.effects <- d[,.(latitude, pH, moisture, pC, cn, doy, NPP, MAT, MAP, MAT_CV, MAP_CV)]

###HERE is where you sequentially fit models - loop begins.
ted.out <- list()
for(i in 1:length(dependent.variables)){
  #specify y-values.
  y <- dependent.variables[,i, with = F]
  
  #find best model combination using stepAIC function and lm.
  #generate a model formula - linear combination of all fixed effects.
  formula <- as.formula(paste('logit(',paste(colnames(y)),') ~', paste(colnames(fixed.effects), collapse="+")))
  #fit model - need complete case data for stepAIC to work without errors.
  dat <- cbind(y,fixed.effects)
  dat <- dat[complete.cases(dat),]
  mod <- lm(formula, data = dat)
  #perform stepAIC on model, subset fixed effects to only include these predictors.
  s.AIC <- data.frame(summary(stepAIC(mod))[[4]])[-1,]
  fixed.effects.jags <- fixed.effects[,rownames(s.AIC), with = F]
  
  #write beta model as a file.
  write_jags_beta(fixed.effects.jags, jags.model.path)
  
  #get complete cases of y and predictors as a jags data object.
  jd <- jd_format(y,fixed.effects.jags)
  
  #parameters to monitor - all predictors plus intercept and dic.
  to.monitor <- c(paste("m",0:length(fixed.effects.jags),sep=""),'dic')
  
  #Run the model.
  jags.out <- run.jags(model = jags.model.path,
                       data=jd,
                       adapt = 100,
                       burnin = 1000,
                       sample = 1000,
                       n.chains=3,
                       monitor = to.monitor)
  
  #summarize output.
  sum <- data.frame(summary(jags.out))
  sum$pred.name <- c('intercept',colnames(fixed.effects.jags))
  #mod.dic <- extract(jags.out, 'dic') #model dic, if you want it. requires sampling more.
  
  #save to output list.
  pre.list <- list(sum, jd)
  names(pre.list) <- c('summary','data')
  ted.out[[i]] <- pre.list
  
  #report.
  cat(i, 'of', ncol(dependent.variables),'models fit.')
  
} #end model fitting loop.

#name and save model output list.
names(ted.out) <- colnames(dependent.variables)
saveRDS(ted.out, output.path)