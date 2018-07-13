#' return_pD.r
#' Gets site-level ddirch model pD value.
#'
#' @param y            #matrix of observed relative abundances used to fit the original model.
#' @param x_mu         #matrix of predictors used to fit original model. No missing values allowed!
#' @param x.mm         #parameter matrix returned from original model.
#' @param alpha        #alpha parameter of original model.
#' @param silent.jags  #if TRUE then JAGS model fitting will shut up.
#'
#' @return             #returns the pD value.
#' @export
#'
#' @examples
#' y=y;x=x_mu; x.mm=fit$x.mm;alpha=fit$alpha;silent.jags=F
return_pD <- function(y, x_mu, x.mm, alpha, silent.jags = F){
  #Specify the JAGS model. This is all it is. can add in variation in x-predictors, return a pD probability distribution, but its really narrow.
  y <- as.matrix(y)
  x_mu <- as.matrix(x_mu)
  
  fix.model = "
  model{
  #mean center all predictors (except intercept).
  for(i in 1:N){
    x.center[i,1] <- 1
      for(j in 2:N.preds){
        x.center[i,j] <- x[i,j] - mean(x[,j])
    }
  }

  #fit species abundances as a linear combination of predictors and parameters.
  for(i in 1:N){
    for(j in 1:N.spp){
      log(a0[i,j]) <- alpha + inprod(x.mm[,j], x.center[i,])
    }
    y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp]) 
  }
  
  } #close model loop.
  "
  
  #setup new JAGS data object with parameters fixed at means.
  new.data <- list(y=y, x=x_mu, N=nrow(y), N.spp = ncol(y), N.preds = ncol(x_mu), x.mm=x.mm, alpha = alpha)
  
  #fit the JAGS model.
  jags.out <- runjags::run.jags(   model = fix.model,
                                   data = new.data,
                                   adapt = 0,
                                   burnin = 0,
                                   sample = 100,
                                   n.chains = 3,
                                   method = 'rjags',
                                   silent.jags = silent.jags,
                                   monitor = c('deviance'))
  jags.sum <- summary(jags.out)
  pD <- jags.sum[1,3]
  return(pD)
}

