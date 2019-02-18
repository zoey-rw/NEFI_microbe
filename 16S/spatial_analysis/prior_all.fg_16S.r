#Fit dirlichet models to all functional groups at once. 
#No hierarchy required, as everything is observed at the site level. Each observation is a unique site.
#clear environment
rm(list = ls())
library(data.table)
library(future)
library(furrr) 
library(magrittr)
library(doParallel)
source('NEFI_functions/tic_toc.r')
source('paths.r')
#source('NEFI_functions/ddirch_site.level_JAGS.r')
library(RCurl)
# source function from colins github
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/ddirch_site.level_JAGS.r", ssl.verifypeer = FALSE)
eval(parse(text = script))
# source paths.r from colins github
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/paths.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#set output path.----
output.path <- prior_16S_all.fg.groups_JAGSfits.path

#load Bahram metadata.----
d <- data.table(readRDS(bahram_metadata.path))

# load covariate selection data.
# covs <- readRDS(bahram_16S_prior_fg_cov.selection_JAGS)

# load Bahram functional group abundances
a1 <- readRDS(prior_N_cyclers_abundances.path)
a2 <- readRDS(prior_C_cyclers_abundances.path)
a3 <- list(readRDS(prior_cop_olig_abundances.path))

# combine these three lists of lists - just get abundances
a <- do.call(c, list(a1, a2, a3))
y <- sapply(a, "[[", 1)

# set pathway names
group_names <- list()
for (i in 1:12) {
  group_names[[i]] <- colnames(a[[i]])[2]
}
group_names[[12]] <- "Cop_olig" #Cop_olig has one more column than the other 11 

#subset to predictors of interest, complete case the thing.
d <- d[,.(Run,pC,cn,pH,moisture,NPP,map,mat,forest,conifer,relEM, Ca, Mg, P, K)] #with micronutrients.
#d <- d[,.(Run,pC,cn,pH,moisture,NPP,map,mat,forest,conifer,relEM)]
d <- d[complete.cases(d),] #optional. This works with missing data.
d <- d[d$Run %in% rownames(y[[1]]),]

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
#covs <- do.call(c, covs)
#names(covs) <- group_names
# covs <- sapply(covs, "[[", 1)
# covariates <- c("intercept",rownames(covs[[i]])) # N_cyclers are first item; each of 7 pathways has different covariates
# cols <- which(colnames(x) %in% covariates)
# x.cov_select <- x[,cols, with=FALSE]
# 
# covs.no.nutr <- sapply(covs, "[[", 2)
# covariates <- c("intercept",rownames(covs.no.nutr[[i]])) # N_cyclers are first item; each of 7 pathways has different covariates
# cols <- which(colnames(x) %in% covariates)
# x.cov_select.no.nutr <- x[,cols, with=FALSE]

#define multiple subsets
x.no.nutr <- x[,.(intercept,pC,cn,pH,pN,moisture,NPP,mat,map,forest,conifer,relEM)]
x.all  <- x[,.(intercept,pC,cn,pH,Ca,Mg,P,K,pN,moisture,NPP,mat,map,forest,conifer,relEM)] # all nutrients + moisture
x.list <- list(
  #x.cov_select,x.cov_select.no.nutr,
  x.no.nutr,
  x.all)

#fit model using function.
#for running production fit on remote.
cat('Begin model fitting loop...\n')
tic()
output.list<-
  foreach(i = 1:length(y)) %dopar% {
    y.group <- y[[i]]
    y.group <- y.group[rownames(y.group) %in% d$Run,]
    y.group <- y.group[match(d$Run, rownames(y.group)),] #order abundance table to match the metadata file
    y.group <- y.group + 1
    y.group <- as.data.frame(y.group)
    y.group <- y.group/rowSums(y.group)
    if(!sum(rownames(y.group) == d$Run) == nrow(y.group)){
      cat('Warning. x and y covariates not in the same order!')
    }
    
    output <-
      (1:length(x.list)) %>%	
      future_map(function(i){
        fit <- site.level_dirlichet_jags(y=y.group,x_mu=x.list[i],adapt = 200, burnin = 5000, sample = 1000, parallel = T)
        return(fit)
      })
    names(output) <- c("no.nutr.preds","all.preds")
    return(output)                                                  #allows nested loop to work.
  }
cat('Model fitting loop complete! ')
toc()

#name the items in the list
names(output.list) <- group_names

cat('Saving fit...\n')
saveRDS(output.list, output.path)
cat('Script complete. \n')