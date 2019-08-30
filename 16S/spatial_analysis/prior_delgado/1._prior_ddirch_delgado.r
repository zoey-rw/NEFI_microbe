# fit prior models for bacterial taxa from Delgado-Baquerizo et al. 2018
rm(list = ls())
library(data.table)
library(doParallel)
source('NEFI_functions/tic_toc.r')
source('NEFI_functions/crib_fun.r')
source('paths.r')
#source('NEFI_functions/ddirch_site.level_JAGS.r')

library(RCurl)
# source function from colins github
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/ddirch_site.level_JAGS.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)


y.all <- readRDS("/projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/prior_abundance_mapping/Delgado/delgado_16S_common_phylo_fg_abun.rds")
d <- readRDS("/projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/prior_abundance_mapping/Delgado/delgado_metadata_spatial.rds")

# set predictors of interest
preds <- c("pC","cn","pH","NPP","forest","map","mat","ndep.glob","relEM")

# format data
x <- d
x <- x[,colnames(x) %in% preds]
x$map <- x$map/1000
intercept <- rep(1, nrow(x))
x <- cbind(intercept, x)
y.all$Species <- NULL

#fit model using function.
#for running production fit on remote.
allfits <- list()
cat('Begin model fitting loop...\n')
tic()
#for (i in 1:length(y.all)){
for (i in 1:length(y.all)){
    y <- y.all[[i]]$rel.abundances
y <- y[,colnames(y)!="other", drop=FALSE]
y <- as.data.frame(y)
y$other <- 1-rowSums(y)
y <- as.matrix(y)
y <- crib_fun(y)
y <- as.data.frame(y)

x <- x[complete.cases(x),]
y <- y[complete.cases(y),]
x <- x[rownames(x) %in% rownames(y),]
y <- y[rownames(y) %in% rownames(x),]
y <- y[match(rownames(x), rownames(y)),] #order abundance table to match the metadata file
if(!sum(rownames(y) == rownames(x)) == nrow(y)){
  cat('Warning. x and y covariates not in the same order!')
}

fit <- site.level_dirlichet_jags(y=y,x_mu=x,
                                 adapt = 200, burnin = 500, sample = 200, 
                                 parallel = T, parallel_method="parallel")
allfits[[i]] <- fit
cat(paste("Model fit for", names(y.all)[i], "\n"))
cat('Model fitting loop complete! ')
toc()
}

saveRDS(allfits, "/projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/JAGS_output/prior_delgado/dir_delgado_8-29-19.rds")
