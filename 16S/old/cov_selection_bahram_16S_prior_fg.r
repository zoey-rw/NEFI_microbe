#covariate selection on bahram 16Sprior functional groups.
#linear combination of as many predictors as we can think of, no interactions.
#clear environment, load paths and functions.
rm(list=ls())
source('NEFI_functions/covariate_selection_JAGS.r')
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/tic_toc.r')
source('paths.r')
library(RCurl)
# source function from colins github
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/ddirch_site.level_JAGS.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#load data.
d <- data.table::data.table(readRDS(bahram_metadata.path))
d <- d[,.(Run,pC,cn,pH,Ca,Mg,P,K,NPP,map,mat,forest,conifer,relEM)]
d <- d[complete.cases(d),] #optional. This works with missing data.

output.path <- bahram_16S_prior_fg_cov.selection_JAGS

# load Bahram functional group abundances
a1 <- readRDS(prior_N_cyclers_abundances.path)
a2 <- readRDS(prior_C_cyclers_abundances.path)
a3 <- list(readRDS(prior_cop_olig_abundances.path))

# combine these three lists of lists - just get abundances
a <- do.call(c, list(a1, a2, a3))
fg_all <- sapply(a, "[[", 1)


#Drop in intercept, setup predictor matrix.
d <- d[d$Run %in% rownames(a[[1]][[1]])]
d$intercept <- rep(1,nrow(d))
d$map <- log(d$map)
x <- as.data.frame(d[,.(intercept,pC,cn,pH,Ca,Mg,P,K,NPP,map,mat,forest,conifer,relEM)])
x.no.nutr <- as.data.frame(d[,.(intercept,pC,cn,pH,NPP,map,mat,forest,conifer,relEM)])
x.list <- list(x,x.no.nutr)
#fg_all <- c("N_cyclers", "C_cyclers", "Cop_olig")
all_fg_output <- list()

all_fg_output<-
  foreach(f = 1:length(fg_all)) %dopar% {

  #organize y data
  y <- fg_all[f]
  y <- y[[1]]
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
  
  fg_output <- list()
  for(k in 1:length(x.list)){
    fg_output[k] <- covariate_selection_JAGS(y=y,x_mu=x.list[[k]], n.adapt = 300, n.burnin = 5000, 
                                   n.sample = 1000,parallel = F)
  cat(paste0("Covariate selection loop #",k, " for group: ", colnames(y)[2]))
  }
  names(fg_output) <- c(paste0("Cov_select ", colnames(y)[2]), 
                        paste0("No.nutr cov_select ", colnames(y)[2])
                        )
  cat(paste0("All covariate selection complete for group: ", colnames(y)[2]))
  #all_fg_output[[f]] <- fg_output
  return(fg_output)
  
} # end functional group loop


# set pathway names
group_names <- list()
for (i in 1:12) {
  group_names[[i]] <- colnames(fg_all[[i]])[2]
}
group_names[[12]] <- "Cop_olig" #Cop_olig has one more column than the other 11 
# add names to each item
names(all_fg_output) <- group_names

# save the output
saveRDS(all_fg_output, bahram_16S_prior_fg_cov.selection_JAGS)
