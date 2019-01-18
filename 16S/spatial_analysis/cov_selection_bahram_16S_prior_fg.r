#covariate selection on bahram 16Sprior functional groups.
#linear combination of as many predictors as we can think of, no interactions.
#clear environment, load paths and functions.
rm(list=ls())
source('NEFI_functions/covariate_selection_JAGS.r')
source('NEFI_functions/crib_fun.r')
source('paths.r')

#load data.
d <- data.table::data.table(readRDS(bahram_metadata.path))
d <- d[,.(Run,pC,cn,PH,Ca,Mg,P,K,pN,moisture,NPP,map,mat,forest,conifer,relEM)]
d <- d[complete.cases(d),] #optional. This works with missing data.

fg_all <- c("N_cyclers", "C_cyclers", "Cop_olig")
all_fg_output <- list()

for (f in 1:length(fg_all)) {
fg_output <- list()
fg <- fg_all[f]
abundances <- switch(fg,
                     "N_cyclers" = readRDS(prior_N_cyclers_abundances.path),
                     "C_cyclers" = readRDS(prior_C_cyclers_abundances.path),
                     "Cop_olig" = readRDS(prior_cop_olig_abundances.path)
)

if (fg == "Cop_olig") {
abun <- list()
abun[[1]] <- abundances
abundances <- abun
}

for (p in 1:length(abundances)) {

#organize y data
y <- abundances[[p]]
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

#Drop in intercept, setup predictor matrix.
d$intercept <- rep(1,nrow(d))
x <- d[,.(intercept,pC,cn,PH,Ca,Mg,P,K,pN,moisture,NPP,map,mat,forest,conifer,relEM)]
x$map <- log(x$map)
#y <- as.data.frame(y)
x <- as.data.frame(x)

#run the algorithm.
send_it <- covariate_selection_JAGS(y=y,x_mu=x, n.adapt = 300, n.burnin = 1000, n.sample = 1000, parallel = F)

fg_output[p] <- send_it

} # end pathway loop

all_fg_output[[f]] <- fg_output

} # end functional group loop

#save the output.
saveRDS(all_fg_output, bahram_16S_prior_fg_cov.selection_JAGS)