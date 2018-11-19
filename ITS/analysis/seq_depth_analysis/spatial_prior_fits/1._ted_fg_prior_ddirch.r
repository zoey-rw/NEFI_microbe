#Fit dirlichet models to functional groups of fungi from Tedersoo et al. Temperate Latitude Fungi.
#rarefying to different depths to understand how this affects prediction.
#request at least 18 processors for parallel run.
#clear environment
rm(list = ls())
library(data.table)
library(doParallel)
source('paths.r')
source('NEFI_functions/ddirch_site.level_JAGS.r')

#detect and register cores.----
n.cores <- detectCores()
n.cores <- 18
registerDoParallel(cores=n.cores)

#set output path.----
output.path <- ted_fg_ddirch_fit_seqdepth.path

#load tedersoo data.----
d <- data.table(readRDS(tedersoo_ITS_clean_map.path))
#d <- d[1:35,] #for testing
y <- readRDS(tedersoo_ITS_fg_list.path)
y <- y$abundances
d <- d[,.(SRR.id,pC,cn,pH,moisture,NPP,map,mat,forest,conifer,relEM)]
d <- d[complete.cases(d),] #optional. This works with missing data.
y <- y[rownames(y) %in% d$SRR.id,]

#rarefy species table to multiple depths.
y.frame <- list()
seq.depth <- c(100,200,300,500,1000,2000)
for(i in 1:length(seq.depth)){
  r.y <- vegan::rrarefy(y, seq.depth[i]) + 1 #dd 1, dirichlet cannot handle zeros.
  r.y <- r.y / rowSums(r.y)
  y.frame[[i]] <- as.data.frame(r.y)
}
names(y.frame) <- paste0('d_',seq.depth)


#Drop in intercept, setup predictor matrix.----
#IMPORTANT: LOG TRANSFORM MAP.
d$intercept <- rep(1,nrow(d))
d$map <- log(d$map)
x <- d[,.(intercept,pC,cn,pH,NPP,mat,map,forest,conifer,relEM)]

#fit model using function.----
#This take a long time to run, probably because there is so much going on.
#fit <- site.level_dirlichet_jags(y=y.frame[[1]],x_mu=x,adapt = 50, burnin = 50, sample = 100)
#for running production fit on remote.
output.list<-
  foreach(i = 1:length(y.frame)) %dopar% {
    fit <- site.level_dirlichet_jags(y=y.frame[[i]],x_mu=x,adapt = 200, burnin = 1000, sample = 1000, parallel = F)
    return(fit)
  }

#name the items in the list, save output.----
names(output.list) <- paste0('d_',seq.depth)

cat('Saving fit...\n')
saveRDS(output.list, output.path)
cat('Script complete. \n')
