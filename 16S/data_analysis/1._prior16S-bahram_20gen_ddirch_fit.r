#Fit dirlichet models to cosmopolitan groups of bacteria/archaea from Bahram et al. Temperate Latitude only.
#No hierarchy required, as everything is observed at the site level. Each observation is a unique site.
#Missing data are allowed.
#clear environment
rm(list = ls())
library(data.table)
library(runjags)
library(doParallel)
source('paths.r')
source('NEFI_functions/ddirch_site.level_JAGS.r')
source('NEFI_functions/ddirch_site.level_JAGS_int.only.r')
source('NEFI_functions/crib_fun.r')

#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#set output path.----
output.path <- bahram_16S.prior_12gen_JAGSfit

#load tedersoo data.----
#d <- data.table(readRDS(tedersoo_ITS.prior_for_analysis.path)) #old analysis path.
#d <- data.table(readRDS(tedersoo_ITS.prior_fromSV_analysis.path))
d <- data.table(readRDS(bahram_prior_gen.path))
#d <- data.table(readRDS(tedersoo_ITS_clean_map.path))
y <- readRDS(cosmo_output_16S.path)
y <- y$abundances
d <- d[,.(Run,pC,cn,PH,moisture,NPP,map,mat,forest,conifer,relEM)]
#d <- d[complete.cases(d),] #optional. This works with missing data.
#d <- d[1:3,] #for testing
y <- y[rownames(y) %in% d$Run,]
if(!sum(rownames(y) == d$Run) == nrow(y)){
  cat('Warning. x and y covariates not in the same order!')
}

#Get relative counts by adding 1 to all observations (can't handle zeros).----
y <- y + 1
y <- y/rowSums(y)
y <- as.data.frame(y)

#Drop in intercept, setup predictor matrix.
x <- d
rownames(x) <- x$Run
x$Run <- NULL
intercept <- rep(1, nrow(x))
x <- cbind(intercept, x)

#IMPORTANT: LOG TRANSFORM MAP.
#log transform map, magnitudes in 100s-1000s break JAGS code.
x$map <- log(x$map)

#define multiple subsets
x.clim <- x[,.(intercept,NPP,mat,map)]
x.site <- x[,.(intercept,pC,cn,PH,forest,conifer,relEM)]
x.all  <- x[,.(intercept,pC,cn,PH,NPP,mat,map,forest,conifer,relEM)]
x.list <- list(x.clim,x.site,x.all)

#fit model using function.
#This take a long time to run, probably because there is so much going on.
#fit <- site.level_dirlichet_jags(y=y,x_mu=x,adapt = 50, burnin = 50, sample = 100)
#for running production fit on remote.
output.list<-
  foreach(i = 1:length(x.list)) %dopar% {
    fit <- site.level_dirlichet_jags(y=y,x_mu=x.list[i],adapt = 200, burnin = 1000, sample = 1000, parallel = F)
    return(fit)
  }

#get intercept only fit.
output.list[[length(x.list) + 1]] <- site.level_dirlichet_intercept.only_jags(y=y, silent.jags = T)

#name the items in the list
names(output.list) <- c('climate.preds','site.preds','all.preds','int.only')

cat('Saving fit...\n')
saveRDS(output.list, output.path)
cat('Script complete. \n')