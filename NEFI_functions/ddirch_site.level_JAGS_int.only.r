#' site.level_dirlichet_intercept.only_jags()
#' Fits a site.level dirlichet model with x-values ONLY OBSERVED AT SITE LEVEL.
#' Accounts for missing data.
#'
#' @param y species matrix as a dataframe. Column names are species names.
#' @param adapt   number of adaptive iterations for JAGS simulation.
#' @param burnin  number of burnin iterations for JAGS simulation.
#' @param sample  number of sample iterations for JAGS simulation.
#' @param n.chains number of chains for JAGS simulation.
#' @param parallel whether or not to run JAGS chains in parallel.
#'
#' @return returns a list with the fitted model, a list of species by parameter tables, and matrices of predicted, observed and residual values on the observation scale (0,1)
#' @export
#'
#' @examples
site.level_dirlichet_intercept.only_jags<-function(y,
                                          adapt = 500, burnin = 1000, sample = 2000, n.chains = 3, 
                                          parallel = F, silent.jags = F){
  #Load some important dependencies.
  source('NEFI_functions/crib_fun.r')
  source('NEFI_functions/sd_to_precision.r')
  source('NEFI_functions/precision_matrix_match.r')
  
  #Some checks before we get started.
  y    <- as.data.frame(y)
  y.names <- colnames(y)

  ###massage your data together.
  #deal with zero relative abundances.
  y <- data.frame(y)

  #make sure every else is a matrix.
  y <- as.matrix(y)
  x <- as.matrix(rep(1, length(y)))

  ###setup jags data object.
  jags.data <- list(N = nrow(y), N.spp = ncol(y), #number of observations and number of species
                    y = y, x=x)                   #species matrix, y, x is vector of 1's.
  
  ###specify JAGS model.
  jags.model = "
  model {
  #parameter priors for each species.
  for(j in 1:N.spp){
    x.m[j] ~ dnorm(0, 1.0E-3)
  }

  #fit species abundances as a linear combination of predictors and parameters.
  for(i in 1:N){
    for(j in 1:N.spp){
      log(a0[i,j]) <- x.m[j] * x[i,]
    }
  y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp]) 
  }
  
  } #close model loop.
  "
  
  ###Fit JAGS model.
  #parallel or not parallel.
  run.method <- ifelse(parallel == F,'rjags','rjparallel')
  #run jags model.
  jags.out <- runjags::run.jags(   model = jags.model,
                                   data = jags.data,
                                   adapt = adapt,
                                   burnin = burnin,
                                   sample = sample,
                                   n.chains = n.chains,
                                   method = run.method,
                                   silent.jags = silent.jags,
                                   monitor = c('x.m','deviance'))
  #summarize output
  out <- summary(jags.out)
  
  #grab parmeters by species, make a list of species-parameter dataframes
  output.list <- list()
  for(i in 1:ncol(y)){
    z <- out[grep("^x\\.m\\[",rownames(out), value = T),]
    z <-   z[i,]
    names <- c('intercept')
    z <- c(names,z)
    names(z)[1] <- 'predictor'
    output.list[[i]] <- data.frame(t(z))
  }
  names(output.list) <- y.names
  
  #get deviance score.
  deviance <- out[grep('deviance',rownames(out)),]
  
  #make a super output that also returns model
  super.list <- list(jags.out, output.list,y,deviance)
  names(super.list) <- c('jags_model','species_parameter_output','observed','deviance')
  
  #return model and output
  return(super.list)
  
} #end function.
