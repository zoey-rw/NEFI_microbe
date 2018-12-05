#Fit dirlichet models to functional groups of bacteria/archaea from Bahram et al. Temperate Latitude.
#No hierarchy required, as everything is observed at the site level. Each observation is a unique site.
#Missing data are allowed.
#clear environment
rm(list = ls())
library(data.table)
library(doParallel)
source('paths.r')
source('NEFI_functions/ddirch_site.level_JAGS.r')
library(RCurl)
# source function from colins github
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/ddirch_site.level_JAGS.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

source('NEFI_functions/ddirch_site.level_JAGS_int.only.r')
source('NEFI_functions/tic_toc.r')
source('NEFI_functions/crib_fun.r')

#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#set output path.----
output.path <- bahram_16S_prior_phylo.group_JAGSfits

#load bahram data.----
d <- data.table(readRDS(bahram_metadata.path))
d <- d[,.(Run,moisture,pC,cn,PH,NPP,map,mat,forest,conifer,relEM)]
d <- d[complete.cases(d),] #optional. This works with missing data.

y <- readRDS(bahram_16S_common_phylo_groups_list.path)

#Drop in intercept, setup predictor matrix.
x <- d
rownames(x) <- x$Run
x$Run <- NULL
intercept <- rep(1, nrow(x))
x <- cbind(intercept, x)
#IMPORTANT: LOG TRANSFORM MAP.
#log transform map, magnitudes in 100s-1000s break JAGS code.
x$map <- log(x$map)


#fit model using function.
#for running production fit on remote.
cat('Begin model fitting loop...\n')
tic()
output.list<-
  foreach(i = 1:length(y)) %dopar% {
    y.group <- y[[i]]
    y.group <- y.group$abundances
    y.group <- y.group[rownames(y.group) %in% d$Run,]
    y.group <- y.group + 1
    y.group <- y.group/rowSums(y.group)
    d <- d[order(d$Run,rownames(y.group)),]
    if(!sum(rownames(y.group) == d$Run) == nrow(y.group)){
      cat('Warning. x and y covariates not in the same order!')
    }
    fit <- site.level_dirlichet_jags(y=y.group,x_mu=x,
                                     adapt = 200, burnin = 2000, sample = 1000, 
                                     parallel = T, parallel_method = 'parallel') #setting parallel rather than rjparallel. 
    return(fit)                                                                  #allows nested loop to work.
  }
cat('Model fitting loop complete! ')
toc()

#get intercept only fit.
#output.list[[length(x.list) + 1]] <- site.level_dirlichet_intercept.only_jags(y=y, silent.jags = T)

#name the items in the list
names(output.list) <- names(y)

cat('Saving fit...\n')
saveRDS(output.list, output.path)
cat('Script complete. \n')