#covariate selection on tedersoo ITSprior.
#linear combination of as many predictors as we can think of, no interactions.
#clear environment, load paths and functions.
rm(list=ls())
source('NEFI_functions/covariate_selection_JAGS.r')
source('NEFI_functions/crib_fun.r')
source('paths.r')

#load data.
d <- data.table::data.table(readRDS(ted.ITSprior_data))
d <- d[,.(Ectomycorrhizal,Saprotroph,Pathogen,Arbuscular,cn,pH,moisture,NPP,map,mat,forest,conifer,relEM)]
d <- d[complete.cases(d),] #optional. This works with missing data.
#d <- d[1:35,] #subset to test.

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
y <- as.data.frame(y)
x <- as.data.frame(x)

#run the algorithm.
send_it <- covariate_selection_JAGS(y=y,x_mu=x, n.adapt = 300, n.burnin = 1000, n.sample = 1000, parallel = T)

#save the output.
saveRDS(send_it, ITS.prior_linear_fg_cov.selection_JAGS)
