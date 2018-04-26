#Fit dirlichet models to functional groups of fungi from Tedersoo et al. Temperate Latitude Fungi.
#Key to get this to work: z-transform predictors.
# still working on how to z-transform precision matrix.
# then need to write in ability to run chains in parallel.
#clear environment
rm(list = ls())
library(data.table)
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/z_transform.r')
source('NEFI_functions/mean_center.r')
source('NEFI_functions/precision_matrix_match.r')

#load tedersoo data.
d <- data.table(readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/ted_all_prior_data.rds'))

#get complete cases.
#Would love to tell JAGS to handle missing independent variables.
d <- d[,.(Ectomycorrhizal, Saprotroph, Arbuscular, Pathogen, mat, map, mat_sd, map_sd, cn, pH, doy)]
d <- d[complete.cases(d),]

#define matrix of dependent variables (y)
y <- d[,.(Ectomycorrhizal, Saprotroph, Arbuscular, Pathogen)]
y <- data.frame(lapply(y, crib_fun)) #this deals with the zero relative abundances that dirlichet doesn't like.

#setup predictors (x)
d$intercept <- rep(1,nrow(d))
x_mu <- d[,.(intercept,mat,map,cn,pH,doy)]

#Get precision matrix for predictors.
#important: need same column names as x_mu
x_precision <- d[,.(mat_sd,map_sd)]
x_precision <- precision_matrix_match(x_mu,x_precision)

#z transform?
z.trans = T
if(z.trans == T){
  sd.out <- data.frame(lapply(x_mu, sd))
    x_mu <- data.frame(lapply(x_mu, z_transform))
  #z-transform SD values in precision matrix.
    for(i in 1:ncol(x_precision)){
      if(mean(x_precision[,i]) != 10000){
        x_precision[,i] <- x_precision[,i] / sd.out[i]
      }
    }
}

#convert SD to precision values
for(i in 1:ncol(x_precision)){
  if(mean(x_precision[,i]) != 10000){
    x_precision[,i] <- 1/(x_precision[,i]^2)
  }
}

jags.model = "
model {
    #setup parameter priors for each species * predictor combination.
    for(j in 1:N.spp){
      for(k in 1:N.preds){
        #m[k,j] ~ dgamma(1.0E-3, 1.0E-3)
        m[k,j] ~ dnorm(0, 1.0E-4)
      }
    }

    #account for uncertainty in predictors
    for(j in 1:N.preds){
      for(i in 1:N){
        x[i,j] ~ dnorm(x_mu[i,j], x_precision[i,j])
      }
    }

    #go ahead and fit means of species abundances as a linear combination of predictors and parameters.
    for(i in 1:N){
        for(j in 1:N.spp){
             log(a0[i,j]) <- inprod(m[,j] , x[i,])
           }
    y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp]) 
    }

} #close model loop.
"

#specify JAGS data object.
jags.data <- list(y = as.matrix(y), x_mu = as.matrix(x_mu), x_precision = as.matrix(x_precision),
                  N.spp = ncol(y), N.preds = ncol(x_mu), N = nrow(y))

#fit runjags model.
jags.out <- runjags::run.jags(model = jags.model,
                               data = jags.data,
                              adapt = 100,
                             burnin = 300,
                             sample = 800,
                           n.chains = 3,
                            monitor = c('m'))
out <- summary(jags.out)

#grab parmeters by species, make a list of species-parameter dataframes
output.list <- list()
for(i in 1:ncol(y)){
  z <- out[grep(paste0(',',i,']'),rownames(out)),]
  z <- cbind(names(x_mu),z)
  colnames(z)[1] <- 'predictor'
  output.list[[i]] <- data.frame(z)
}
#Species names are the list entry names.
names(output.list) <- colnames(y)

#did we get relative abundances right intercept only case?
#colMeans(y)
#a.m <- exp(out[,2])
#a.m / sum(a.m)
