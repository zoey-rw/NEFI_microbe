#covariate selection on bahram 16Sprior functional groups.
#linear combination of as many predictors as we can think of, no interactions.
#clear environment, load paths and functions.
rm(list=ls())
library(runjags)
library(doParallel)
library(future)
library(furrr) 
library(magrittr)
plan(multiprocess)
source('NEFI_functions/covariate_selection_JAGS.r')
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/tic_toc.r')
source('paths.r')

#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#load data.
d <- data.table::data.table(readRDS(bahram_metadata.path))
d <- d[,.(Run,pC,cn,pH,Ca,Mg,P,K,pN,moisture,NPP,map,mat,forest,conifer,relEM)]
d <- d[complete.cases(d),] #optional. This works with missing data.

output.path <- bahram_16S_prior_fg_cov.selection_JAGS
all_fg_output <- list()

# load Bahram functional group abundances
a1 <- readRDS(prior_N_cyclers_abundances.path)
a2 <- readRDS(prior_C_cyclers_abundances.path)
a3 <- list(readRDS(prior_cop_olig_abundances.path))

# combine these three lists of lists - just get abundances
a <- do.call(c, list(a1, a2, a3))
a <- sapply(a, "[[", 1)

tic()
#for (p in 1:length(a)) {
  # loop over models
  output <- list()
  output <-
    foreach(p = 1:length(a)) %dopar% { # loop through each functional group
  
#organize y data
y <- a[[p]]
y <- y[rownames(y) %in% d$Run,]

#order abundance table to match the metadata file and vice versa
d <- d[d$Run %in% rownames(y),]
y <- y[match(d$Run, rownames(y)),]
if(!sum(rownames(y) == d$Run) == nrow(y)){
  cat('Warning. x and y covariates not in the same order!')
}

#Get relative counts by adding 1 to all observations (can't handle zeros).----
y <- y + 1
y <- y/rowSums(y)
y <- as.data.frame(y)

#Drop in intercept, setup predictor matrix.
d$intercept <- rep(1,nrow(d))
d$map <- log(d$map)

x <- as.data.frame(d[,.(intercept,pC,cn,pH,Ca,Mg,P,K,moisture,NPP,map,mat,forest,conifer,relEM)])
x.no.nutr <- as.data.frame(d[,.(intercept,pC,cn,pH,moisture,NPP,map,mat,forest,conifer,relEM)])

#run the algorithm.
covs <- covariate_selection_JAGS(y=y,x_mu=x, n.adapt = 300, n.burnin = 1000, 
                                 n.sample = 1000,parallel = F)
cat(paste0("First covariate selection complete for group: ", colnames(y)[2]))
covs.no.nutr <- covariate_selection_JAGS(y=y,x_mu=x.no.nutr, n.adapt = 300, n.burnin = 1000, 
                                         n.sample = 1000, parallel = F)
cat(paste0("Second (no.nutr) covariate selection complete for group: ", colnames(y)[2]))

cov_output <- list(covs, covs.no.nutr)
names(cov_output) <- c(paste0("Cov_select ", colnames(y)[2]), paste0("No.nutr cov_select ", colnames(y)[2]))
cat(paste0("All covariate selection complete for group: ", colnames(y)[2]))
return(cov_output)
} # end functional group loop
toc()
# change name of last item 
names(output)[[12]] <- "Cop_olig"

# save the output
saveRDS(output, output.path)
