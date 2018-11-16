#Fit MULTINOMIAL dirlichet models to functional groups of fungi from Tedersoo et al. Temperate Latitude Fungi.
#No hierarchy required, as everything is observed at the site level. Each observation is a unique site.
#Missing data are allowed.
#clear environment
rm(list = ls())
library(data.table)
library(doParallel)
source('paths.r')
source('NEFI_functions/dmulti-ddirch_site.level_JAGS.r')
source('NEFI_functions/crib_fun.r')

#detect and register cores.----
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#set output path.----
output.path <- ted_ITS.prior_dmulti.ddirch_fg_JAGSfit

#load tedersoo data.----
d <- data.table(readRDS(tedersoo_ITS_clean_map.path))
y <- readRDS(tedersoo_ITS_fg_list.path)
seq.depth <- y$seq_total
y <- y$abundances
#d <- d[,.(Ectomycorrhizal,Saprotroph,Pathogen,Arbuscular,pC,cn,pH,moisture,NPP,map,mat,forest,conifer,relEM)]
d <- d[,.(SRR.id,pC,cn,pH,moisture,NPP,map,mat,forest,conifer,relEM)]
d <- d[complete.cases(d),] #optional. This works with missing data.
y <- y[rownames(y) %in% d$SRR.id,]
#y$Arbuscular <- NULL
y <- as.data.frame(y)
testing = F
if(testing == T){
  d <- d[1:35,] #for testing
  y <- y[rownames(y) %in% d$SRR.id,]
}

#organize y data----
#y <- d[,.(Ectomycorrhizal,Saprotroph,Pathogen)]
#y$other <- 1 - rowSums(y)
#y <- data.frame(lapply(y,crib_fun, N = nrow(y) * ncol(y)))
#in the rare case where one column actually needs to be a zero for a row to prevent to summing over 1...
#for(i in 1:nrow(y)){
#  if(rowSums(y[i,]) > 1){
#    y[i,] <- y[i,] / rowSums(y[i,])
#  }
#}

#Drop in intercept, setup predictor matrix.----
#IMPORTANT: LOG TRANSFORM MAP.
d$intercept <- rep(1,nrow(d))
d$map <- log(d$map)
x <- d[,.(intercept,pC,cn,pH,moisture,NPP,mat,map,forest,conifer,relEM)]

#define multiple subsets
x.clim <- d[,.(intercept,NPP,mat,map)]
x.site <- d[,.(intercept,pC,cn,pH,forest,conifer,relEM)]
x.all  <- d[,.(intercept,pC,cn,pH,NPP,mat,map,forest,conifer,relEM)]
x.list <- list(x.clim,x.site,x.all)

#fit model using function.----
#This take a long time to run, probably because there is so much going on.
#fit <- site.level_dirlichet_jags(y=y,x_mu=x,adapt = 50, burnin = 50, sample = 100)
#for running production fit on remote.
output.list<-
  foreach(i = 1:length(x.list)) %dopar% {
    fit <- site.level_multi.dirich_jags(y=y,seq.depth=seq.depth,x_mu=x.list[i],adapt = 200, burnin = 16000, sample = 5000, parallel = T)
    return(fit)
  }

#get intercept only fit.
#output.list[[length(x.list) + 1]] <- site.level_dirlichet_intercept.only_jags(y=y, silent.jags = T)

#name the items in the list, save output.----
#names(output.list) <- c('climate.preds','site.preds','all.preds','int.only')
names(output.list) <- c('climate.preds','site.preds','all.preds')

cat('Saving fit...\n')
saveRDS(output.list, output.path)
cat('Script complete. \n')
