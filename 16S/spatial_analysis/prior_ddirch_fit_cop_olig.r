#Fit dirlichet models to oligotrophic/copiotrophic groups of bacteria/archaea from Bahram et al. 
#No hierarchy required, as everything is observed at the site level. Each observation is a unique site.

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
source('NEFI_functions/ddirch_site.level_JAGS.r')
source('NEFI_functions/ddirch_site.level_JAGS_int.only.r')
source('NEFI_functions/crib_fun.r')


#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

# model these groups individually?
indiv <- F

#set output path.----
output.path <- bahram_16S.prior_cop_olig_JAGSfit 
#output.path <- bahram_16S.prior_cop_olig_indiv_JAGSfit
  
#load tedersoo data.----
d <- data.table(readRDS(bahram_metadata.path))
d <- d[,.(Run,pC,cn,PH,Ca,Mg,P,K,pN,moisture,NPP,map,mat,forest,conifer,relEM)] # all nutrients
d <- d[complete.cases(d),] #optional. This works with missing data.

# load covariate selection data.
covs <- readRDS(bahram_16S_prior_fg_cov.selection_JAGS)

if (indiv == T) {
  a <- readRDS(prior_cop_olig_abundances_indiv.path) 
  group_names <- c("copiotroph","oligotroph")
  all_cop_olig_models <- vector("list", 2)
  
   } else {
     y <- readRDS(prior_cop_olig_abundances.path)
       a <- list()
       a[[1]] <- y
       group_names <- c("cop_olig")
       all_cop_olig_models <- vector("list", 1)
     }

for (i in 1:length(a)) {
  y <- a[[i]]

y <- y$abundances

y <- y[rownames(y) %in% d$Run,]
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
# get covariates from model selection output, + intercept column
covariates <- c("intercept", rownames(covs[[3]][[i]])) # Cop_olig model is third item; one set of covariates
cols <- which(colnames(x) %in% covariates)
x.cov_select <- x[,cols, with=FALSE]
x.all  <- x[,.(intercept,pC,cn,PH,Ca,Mg,P,K,pN,moisture,NPP,mat,map,forest,conifer,relEM)]
x.list <- list(x.cov_select,x.all)

#fit model using function.
#This take a long time to run, probably because there is so much going on.
#fit <- site.level_dirlichet_jags(y=y,x_mu=x,adapt = 50, burnin = 50, sample = 100)
#for running production fit on remote.
output.list<-
(1:length(x.list)) %>%	
future_map(function(i){
    fit <- site.level_dirlichet_jags(y=y,x_mu=x.list[i],adapt = 200, burnin = 1000, sample = 1000, parallel = T)
    return(fit)
  })

#get intercept only fit.
output.list[[length(x.list) + 1]] <- site.level_dirlichet_intercept.only_jags(y=y, silent.jags = T)

#name the items in the list
names(output.list) <- c('cov_select','all.preds','int.only')

all_cop_olig_models[[i]] <- output.list
cat(paste0('Cop-olig fit completed for ', group_names[[i]]))
}


cat('Saving fit...\n')
saveRDS(all_cop_olig_models, output.path)
cat('Script complete. \n')