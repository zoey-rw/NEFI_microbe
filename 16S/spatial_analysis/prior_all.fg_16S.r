# Fit dirlichet models to all functional groups of bacteria/archaea from Bahram et al. 
# Model fits separately for each group, since they are not mutually-exclusive.
# No hierarchy required, as everything is observed at the site level. Each observation is a unique site.

#clear environment
rm(list = ls())
library(data.table)
library(runjags)
library(doParallel)
library(future)
library(furrr) 
library(magrittr)
plan(multiprocess)
source('paths.r')
source('NEFI_functions/ddirch_site.level_JAGS_int.only.r')
source('NEFI_functions/crib_fun.r')
#source('NEFI_functions/ddirch_site.level_JAGS.r')
library(RCurl)
# source function from colins github
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/ddirch_site.level_JAGS.r", ssl.verifypeer = FALSE)
eval(parse(text = script))
source('NEFI_functions/tic_toc.r')

#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#set output path.----
output.path <- prior_16S_all.fg.groups_JAGSfits.path

#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#load Bahram metadata.----
m <- data.table(readRDS(bahram_metadata.path))

# load covariate selection data.
covs <- readRDS(bahram_16S_prior_fg_cov.selection_JAGS)

# load Bahram functional group abundances
a1 <- readRDS(prior_N_cyclers_abundances.path)
a2 <- readRDS(prior_C_cyclers_abundances.path)
a3 <- list(readRDS(prior_cop_olig_abundances.path))

# combine these three lists of lists - just get abundances
a <- do.call(c, list(a1, a2, a3))
a <- sapply(a, "[[", 1)

# combine list of covariates
covs <- do.call(c, covs)
names(covs) <- group_names

# set pathway names
group_names <- list()
for (i in 1:11) {
  group_names[[i]] <- colnames(a[[i]])[2]
}
group_names[[12]] <- "Cop_olig" #Cop_olig has one more column than the other 11 

# subset covariate dataset
d <- m[,.(Run,pC,cn,pH,Ca,Mg,P,K,pN,moisture,NPP,map,mat,forest,conifer,relEM)] # all nutrients, no moisture
d <- d[complete.cases(d),] #optional. This works with missing data.

# create vector for model output
all_fg_models <- vector("list", 12)

# loop through each functional group
for (i in 1:length(a)) {
  y <- a[[i]]
  y <- y[rownames(y) %in% d$Run,]
  d <- d[d$Run %in% rownames(y),]
  
  #order abundance table to match the metadata file
  y <- y[match(d$Run, rownames(y)),]
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
  
  # define multiple subsets
  # get covariates from model selection output
  covariates <- c("intercept",rownames(covs[[i]])) # N_cyclers are first item; each of 7 pathways has different covariates
  cols <- which(colnames(x) %in% covariates)
  x.cov_select <- x[,cols, with=FALSE]
  x.all  <- x[,.(intercept,pC,cn,pH,Ca,Mg,P,K,pN,moisture,NPP,mat,map,forest,conifer,relEM)] # all nutrients + moisture
  x.list <- list(x.cov_select,x.all)
  
  #fit model using function.
  #This take a long time to run, probably because there is so much going on.
  #fit <- site.level_dirlichet_jags(y=y,x_mu=x,adapt = 50, burnin = 50, sample = 100)
  #for running production fit on remote.
  output.list<-
    (1:length(x.list)) %>%	
    future_map(function(i){
      fit <- site.level_dirlichet_jags(y=y,x_mu=x.list[i],adapt = 200, burnin = 5000, sample = 1000, parallel = T)
      return(fit)
    })
  
  #get intercept only fit.
  output.list[[length(x.list) + 1]] <- site.level_dirlichet_intercept.only_jags(y=y, silent.jags = T)
  
  #name the items in the list
  names(output.list) <- c('cov_select','all.preds','int.only')
  all_fg_models[[i]] <- output.list
  cat(paste0('Prior fit completed for ', pathway_names[[i]]))
}

cat('Saving fit...\n')
saveRDS(all_fg_models, output.path)
cat('Script complete. \n')
