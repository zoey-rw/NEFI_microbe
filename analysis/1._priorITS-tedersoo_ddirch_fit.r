#Fit dirlichet models to functional groups of fungi from Tedersoo et al. Temperate Latitude Fungi.
#No hierarchy required, as everything is observed at the site level. Each observation is a unique site.
#Missing data are allowed.
#clear environment
rm(list = ls())
library(data.table)
source('paths.r')
source('NEFI_functions/ddirch_site.level_JAGS.r')
source('NEFI_functions/crib_fun.r')

#load tedersoo data.
d <- data.table(readRDS(ted.ITSprior_data))
d <- d[,.(Ectomycorrhizal,Saprotroph,Pathogen,Arbuscular,cn,pH,moisture,NPP,map,mat,forest,conifer,relEM)]
d <- d[complete.cases(d),] #optional. This works with missing data.

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
x <- d[,.(intercept,cn,pH,moisture,NPP,mat,map,forest,conifer,relEM)]
x$map <- log(x$map)

#fit model using function.
#This take a long time to run, probably because there is so much going on.
#fit <- site.level_dirlichet_jags(y=y,x_mu=x,adapt = 50, burnin = 50, sample = 100)
#for running production fit on remote.
fit <- site.level_dirlichet_jags(y=y,x_mu=x,adapt = 200, burnin = 1000, sample = 1000, parallel = T)

cat('Saving fit...\n')
saveRDS(fit, ted_ITS.prior_fg_JAGSfit)
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
