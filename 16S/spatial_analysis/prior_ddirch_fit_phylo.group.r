#Fit dirlichet models to phylogenetic groups of bacteria/archaea from Bahram et al. Temperate Latitude.
#No hierarchy required, as everything is observed at the site level. Each observation is a unique site.
#Missing data are allowed.
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

#source('NEFI_functions/ddirch_site.level_JAGS_int.only.r')
source('NEFI_functions/tic_toc.r')

#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#set output path.----
#output.path <- bahram_16S_prior_phylo.group_JAGSfits
output.path <- "/projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/JAGS_output/prior_phylo_JAGSfit_phylumtest.rds"

#load bahram abundance data.----
y <- readRDS(bahram_16S_common_phylo_groups_list.path)

#load bahram metadata and format----
d <- data.table(readRDS(bahram_metadata.path))
#subset to predictors of interest, complete case the thing.
d <- d[,.(Run,pC,cn,pH,NPP,map,mat,forest,conifer,relEM, Ca, Mg, P, K)]
d <- d[complete.cases(d),] #optional. This works with missing data.
d <- d[d$Run %in% rownames(y[[1]]$abundances),]

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
x.no.nutr <- x[,.(intercept,pC,cn,pH,NPP,mat,map,forest,conifer,relEM)]
x.all  <- x[,.(intercept,pC,cn,pH,Ca,Mg,P,K,NPP,mat,map,forest,conifer,relEM)] # all nutrients + moisture
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
    y.group <- y[[i]]$abundances
    y.group <- y.group[rownames(y.group) %in% d$Run,]
    y.group <- y.group + 1
    y.group <- y.group/rowSums(y.group)
    y.group <- y.group[order(rownames(y.group),d$Run),]
    if(!sum(rownames(y.group) == d$Run) == nrow(y.group)){
      cat('Warning. x and y covariates not in the same order!')
    }
    
    output <- list()
    for(k in 1:length(x.list)){
      fit <- site.level_dirlichet_jags(y=y.group,x_mu=x.list[[k]],
                                       adapt = 200, burnin = 8000, sample = 1000, 
                                       parallel = T, parallel_method="parallel")
      output[[k]] <- fit
    }
    names(output) <- c("no.nutr.preds","all.preds")
    cat(paste("Models fit for", names(y)[i], "\n"))
    return(output) #allows nested loop to work.
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