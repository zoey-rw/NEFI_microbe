#' linear_dirlichet_jags
#' Fits a multivariate dirlichet model to species relative abundances as a linear combination of predictors and parameters.
#' Models are fit in JAGS.
#' requires runjags package.
#' requires multiple other functinos in the NEFI_microbe/NEFI_functions directory.
#'
#' @param x_mu data.frame or matrix of predictors, with column names.
#' @param x_sd data.frame or matrix of predictor sd values. Column names should map from x_mu to x_sd as mat -> mat_sd or mat.sd or matSD (or similar)
#' @param y data.frame or matrix of species, with column names.
#' @param n.chains number of chains to run
#' @param adapt number of adaptive iterations
#' @param burnin number of burnin iterations
#' @param sample number of samples to draw post burnin
#'
#' @return returns a list containing the jags model and a list of parameter tables by species.
#' @export
#'
#' @examples
#' #load tedersoo data.
#' library(data.table)
#' d <- data.table(readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/ted_all_prior_data.rds'))
#' 
#' d <- d[,.(Ectomycorrhizal, Saprotroph, Arbuscular, Pathogen, mat, map, mat_sd, map_sd, cn, pH, doy)]
#' d <- d[complete.cases(d),]
#' 
#' #define matrix of dependent variables (y)
#' y <- d[,.(Ectomycorrhizal, Saprotroph)]
#' 
#' #setup predictors (x)
#' d$intercept <- rep(1,nrow(d))
#' x_mu <- d[,.(intercept,map)]
#' 
#' #setup x_sd
#' x_sd <- d[,.(mat_sd,map_sd)]
#' 
#' #Run function.
#' fit_dirlichet_jags(x_mu=x_mu, x_sd=x_sd, y=y, z.trans = T, sample = 800)

fit_dirlichet_jags <- function(x_mu, x_sd, y, z.trans = F, n.chains = 3, adapt = 100, burnin = 200, sample = 400){
  
  #source some other functions.
  source('/home/caverill/NEFI_microbe/NEFI_functions/crib_fun.r')
  source('/home/caverill/NEFI_microbe/NEFI_functions/z_transform.r')
  source('/home/caverill/NEFI_microbe/NEFI_functions/precision_matrix_match.r')
  
  #round out y data with other category.
  #Requires first dealing with zero relative abundances dirlichet does not like.
  y <- data.frame(lapply(y, crib_fun))

  #Match up predictors and their SD. if no SD supplied we assign ~perfect precision.
  x_sd <- precision_matrix_match(x_mu,x_sd)
  
  #z-transform loop. transforms both x obs and associated sd.
  if(z.trans == T){
    sd.out <- data.frame(lapply(x_mu, sd))
    x_mu <- data.frame(lapply(x_mu, z_transform))
    #z-transform SD values in precision matrix.
    for(i in 1:ncol(x_sd)){
      if(mean(x_sd[,i]) != 10000){
        x_sd[,i] <- x_sd[,i] / sd.out[i]
      }
    }
  }
  
  #convert SD to precision values
  for(i in 1:ncol(x_sd)){
    if(mean(x_sd[,i]) != 10000){
      x_sd[,i] <- 1/(x_sd[,i]^2)
    }
  }
  x_precision <- x_sd
  
  #specify JAGS model. This is general for any linear combination of x's.
  jags.model = "
  model {
    #setup parameter priors for each species * predictor combination.
    for(j in 1:N.spp){
      for(k in 1:N.preds){
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
  
  #setup JAGS data object.
  jags.data <- list(y = as.matrix(y), x_mu = as.matrix(x_mu), x_precision = as.matrix(x_precision),
                    N.spp = ncol(y), N.preds = ncol(x_mu), N = nrow(y))
  
  #run jags model.
  jags.out <- runjags::run.jags(   model = jags.model,
                                    data = jags.data,
                                   adapt = adapt,
                                  burnin = burnin,
                                  sample = sample,
                                n.chains = n.chains,
                                 monitor = c('m'))
  #summarize output
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
  
  #make a super output that also returns model
  super.list <- list(jags.out, output.list)
  names(super.list) <- c('jags_model','species_parameter_output')
  
  #return model and output
  return(super.list)
  
} ##end function.
