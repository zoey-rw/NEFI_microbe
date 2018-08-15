#Fit dirlichet models to functional groups of fungi from Tedersoo et al. Temperate Latitude Fungi.
#PROBLEM: something in the site-specific predictors is killing our predicted % ecto. ITS THE MOISTURE VALUES.
#PROBLEM: probably - no %C in model, mismatch in how CN or pC measured in two studies, 
#No hierarchy required, as everything is observed at the site level. Each observation is a unique site.
#Missing data are allowed.
#clear environment
rm(list = ls())
library(data.table)
library(doParallel)
source('paths.r')
source('NEFI_functions/ddirch_site.level_JAGS.r')
source('NEFI_functions/ddirch_site.level_JAGS_int.only.r')
source('NEFI_functions/crib_fun.r')

#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#load tedersoo data.
d <- data.table(readRDS(ted.ITSprior_data))
d <- d[,.(Ectomycorrhizal,Saprotroph,Pathogen,Arbuscular,pC,cn,pH,moisture,NPP,map,mat,forest,conifer,relEM)]
d <- d[complete.cases(d),] #optional. This works with missing data.
#d <- d[1:35,] #for testing

#organize y data
y <- d[,.(Ectomycorrhizal,Saprotroph,Pathogen,Arbuscular)]
#make other column
y <- data.frame(lapply(y,crib_fun))
y$other <- 1 - rowSums(y)
y <- as.data.frame(y)
#reorder columns. other needs to be first.
y <- y[c('other','Ectomycorrhizal','Pathogen','Saprotroph','Arbuscular')]

#Drop in intercept, setup predictor matrix.
d$intercept <- rep(1,nrow(d))
d$map <- log(d$map)
x <- d[,.(intercept,pC,cn,pH,moisture,NPP,mat,map,forest,conifer,relEM)]

#define multiple subsets
x.clim <- d[,.(intercept,NPP,mat,map)]
x.site <- d[,.(intercept,pC,cn,pH,forest,conifer,relEM)]
x.all  <- d[,.(intercept,pC,cn,pH,NPP,mat,map,forest,conifer,relEM)]
x.list <- list(x.clim,x.site,x.all)

#fit model using function.
#This take a long time to run, probably because there is so much going on.
#fit <- site.level_dirlichet_jags(y=y,x_mu=x,adapt = 50, burnin = 50, sample = 100)
#for running production fit on remote.
output.list<-
  foreach(i = 1:length(x.list)) %dopar% {
    fit <- site.level_dirlichet_jags(y=y,x_mu=x.list[i],adapt = 200, burnin = 1000, sample = 1000, parallel = T)
    return(fit)
  }

#get intercept only fit.
output.list[[length(x.list) + 1]] <- site.level_dirlichet_intercept.only_jags(y=y, silent.jags = T)

#name the items in the list
names(output.list) <- c('climate.preds','site.preds','all.preds','int.only')

cat('Saving fit...\n')
saveRDS(output.list, ted_ITS.prior_fg_JAGSfit)
cat('Script complete. \n')

#visualize fits
#par(mfrow = c(1,ncol(fit$observed)))
#for(i in 1:ncol(fit$observed)){
#  plot(fit$observed[,i] ~ fit$predicted[,i], ylim = c(0,1))
#  rsq <- summary(betareg::betareg(fit$observed[,i] ~ fit$predicted[,i]))$pseudo.r.squared
#  txt <- paste0('R2 = ',round(rsq,2))
#  mtext(colnames(fit$predicted)[i], line = -1.5, adj = 0.05)
#  mtext(txt, line = -3.5, adj = 0.05)
#  abline(0,1, lwd = 2)
#  abline(lm(fit$observed[,i] ~ fit$predicted[,i]), lty = 2, col = 'green')
#}
