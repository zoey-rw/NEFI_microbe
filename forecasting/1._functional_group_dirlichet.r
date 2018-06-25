#Fit dirlichet models to functional groups of fungi from Tedersoo et al. Temperate Latitude Fungi.
#No hierarchy required, as everything is observed at the site level. Each observation is a unique site.
#would love to tell function to returns a vector of predicted values and residuals, r2 values.
#Key to get this to work: z-transform predictors.
#clear environment
rm(list = ls())
library(data.table)
source('NEFI_functions/linear_dirlichet_jags.r')

#load tedersoo data.
d <- data.table(readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/ted_all_prior_data.rds'))

#get complete cases.
#Would love to tell JAGS to handle missing independent variables.
d <- d[,.(Ectomycorrhizal, Saprotroph, Arbuscular, Pathogen, mat, map, mat_sd, map_sd, pC, cn, pH, doy)]
d <- d[complete.cases(d),]

#define matrix of dependent variables (y)
y <- d[,.(Ectomycorrhizal, Saprotroph, Arbuscular, Pathogen)]
y <- d[,.(Ectomycorrhizal, Saprotroph)]
y$other <- 1 - rowSums(y)

#setup predictors (x)
d$intercept <- rep(1,nrow(d))
x_mu <- d[,.(intercept,mat,map,pC,cn,pH,doy)]
x_mu <- d[,.(intercept,mat,pC,pH)]
#x_mu <- d[,.(intercept,mat)]

#Get sd matrix for predictors. function will convert to precision.
#important: need column names need to map to x_mu
#if a predictor is named 'mat' in x_mu, then its corresponding sd column must also contain 'mat' somewhere
#can be mat.sd, mat_sd, matSD, SDmat, grep will find it.
#important that other predictors do not have a 'mat' motif anywhere in the name!!!
x_sd <- d[,.(mat_sd,map_sd)]

#run the function. returns a list with the jags model, as well as a list of parameter output by species.
test <- linear_dirlichet_jags(y = y, x_mu = x_mu, x_sd = x_sd, z.trans = F, parallel = F)

#check fits
par(mfrow = c(1,ncol(test$observed)))
for(i in 1:ncol(test$observed)){
  plot(test$observed[,i] ~ test$predicted[,i], ylim = c(0,1))
  rsq <- summary(betareg::betareg(test$observed[,i] ~ test$predicted[,i]))$pseudo.r.squared
  txt <- paste0('R2 = ',round(rsq,2))
  mtext(colnames(test$predicted)[i], line = -1.5, adj = 0.05)
  mtext(txt, line = -3.5, adj = 0.05)
  abline(0,1, lwd = 2)
}
