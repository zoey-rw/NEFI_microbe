#Fit dirlichet models to functional groups of fungi from Tedersoo et al. Temperate Latitude Fungi.
#No hierarchy required, as everything is observed at the site level. Each observation is a unique site.
#Missing data are allowed.
#clear environment
rm(list = ls())
library(data.table)
source('NEFI_functions/linear_dirlichet_jags.r')
source('NEFI_functions/ddirch_site.level_JAGS.r')
source('NEFI_functions/crib_fun.r')

#load tedersoo data.
d <- data.table(readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/ted_all_prior_data.rds'))
d <- d[,.(Ectomycorrhizal,Saprotroph,Arbuscular,cn,pH,moisture,NPP,map,mat)]
d <- d[complete.cases(d),] #optional. This works with missing data.

#organize y data
y <- d[,.(Ectomycorrhizal,Saprotroph,Arbuscular)]
#make other column
y <- data.frame(lapply(y,crib_fun))
y$other <- 1 - rowSums(y)
y <- as.data.frame(y)
#reorder columns. other needs to be first.
y <- y[c('other','Ectomycorrhizal','Saprotroph','Arbuscular')]

#Drop in intercept, setup predictor matrix.
d$intercept <- rep(1,nrow(d))
x <- d[,.(intercept,cn,pH,moisture,NPP,mat,map)]
x$map <- log(x$map)

#fit model using function.
#This take a long time to run, probably because there is so much going on.
fit <- site.level_dirlichet_jags(y=y,x_mu=x,adapt = 50, burnin = 50, sample = 100)

#visualize fits
par(mfrow = c(1,ncol(fit$observed)))
for(i in 1:ncol(fit$observed)){
  plot(fit$observed[,i] ~ fit$predicted[,i], ylim = c(0,1))
  rsq <- summary(betareg::betareg(fit$observed[,i] ~ fit$predicted[,i]))$pseudo.r.squared
  txt <- paste0('R2 = ',round(rsq,2))
  mtext(colnames(fit$predicted)[i], line = -1.5, adj = 0.05)
  mtext(txt, line = -3.5, adj = 0.05)
  abline(0,1, lwd = 2)
  abline(lm(fit$observed[,i] ~ fit$predicted[,i]), lty = 2, col = 'green')
}
