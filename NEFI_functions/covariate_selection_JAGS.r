#' covariate_selection_JAGS
#' algorithm to perform backwards covariate selection based on JAGS model deviance scores.
#' Given a model and a set of predictors the algorithm proceeds this way:
#' Step 1. calculate model deviance with all predictors.
#' Step 2. Rerun the model without a particular predictor, calculate deviance. Do this for all predictors.
#' Step 3. Drop the predictor that results in the largest improvement in deviance (apparently this is the greedy part of a greedy algorithm.)
#' Step 4. Using our new best model, repeat Steps 2 and 3 until:
#'         (a.) deviance no longer improves,
#'         (b.) deviance does not improve beyond a threshold change (say 2 units),
#'         (c.) there are no more predictors to remove.
#'         
#' @param y           #matrix of species relative abundances.
#' @param x_mu        #matrix of covaraites.
#' @param n.adapt     #number of adaptive iterations.
#' @param n.burnin    #number of burnin iterations.
#' @param n.sample    #number of sampling iterations.
#'
#' @return            #returns a list of covariates retained, deviance scores of the last comparison, and the JAGS model that is the best.
#'                    #JAGS model in list is also a list with the full model, model summary, species X parameter tables and predicted/observed/residual matrices.
#' @export
#'
#' @examples
#' #### generate pseudo data to test ####
#'n.obs <- 25
#'x1 <- runif(n.obs,20,40);x2 <- runif(n.obs, 2, 4);x3 <- runif(n.obs,100,110);intercept <- rep(1,n.obs);y <- x1*0.5 + x2*3 + rnorm(n.obs)
#'y1 <- y;y2 <- rep(10, length(y1));spp.y <- data.frame(cbind(y1,y2));spp.y <- spp.y / rowSums(spp.y)
#'x_mu <- data.frame(intercept,x1,x2,x3)
#'test <- covariate_selection_JAGS(y=spp.y, x_mu = x_mu)
#'

#Function depends on these things.
library(runjags)
library(doParallel)
source('NEFI_functions/ddirch_site.level_JAGS.r')
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/tic_toc.r')

covariate_selection_JAGS <- function(y, x_mu, n.adapt = 100, n.burnin = 100, n.sample = 200, parallel = F, silent.jags = T){
  #### Fit first JAGS model to all covariates. ####
  #register parallel cores
  registerDoParallel(cores = detectCores())
  
  #start a loop counter and a done indicator for the while loop.
  donezo <- 0
  loop.count <- 1
  runmode <- ifelse(parallel == F, F, T)
  
  #fit first model, calculate deviance.
  cat('Fitting initial model...\n')
  tic()
  fit <- site.level_dirlichet_jags(y = y, x_mu = x_mu, adapt = n.adapt, burnin = n.burnin, sample = n.sample, parallel = runmode)
  dev <- fit$deviance #grab deviance
  cat('Initial model fit complete.')
  toc()
  while(donezo < 1){
    #sequentially leave out one predictor, except intercept.
    tic()
    fit.list <- list() 
    fit.list <-
      foreach(i = 2:ncol(x_mu)) %dopar% {
        x_new <- x_mu[,-i]
        mod.fit <- site.level_dirlichet_jags(y=y, x_mu = x_new, adapt = n.adapt, burnin = n.burnin, sample = n.sample, parallel = runmode, silent.jags = silent.jags)
        fit.list[[i-1]] <- mod.fit
      }
    #grab deviance statistics as data.frame.
    dev.list <- list()
    for(i in 1:length(fit.list)){
      dev.list[[i]] <- fit.list[[i]]$deviance
    }
    dev.list <-data.frame(do.call(rbind,dev.list))
    rownames(dev.list) <- colnames(x_mu)[2:ncol(x_mu)]
    
    #calculate difference in deviance scores.
    dev.list$diff <- dev[4] - dev.list$Mean
    ###Here I think I should compare deviance distributions, if 95% CI overlaps with reference model deviance, then set difference to zero.
    #removing predictor that results in the largest improvement in deviance score.
    max.diff <- max(dev.list$diff)
    to.drop <- rownames(dev.list[which.max(dev.list$diff),])
    
    #report how long current selection loop took to complete.
    cat('Model selection loop',loop.count,'complete. ')
    toc()
    
    #if the difference is less than 2 or negative you are done.
    if(max.diff < 2){
      cat('*Best* model found! results are in output.list. \n')
      to_keep <- colnames(x_mu)
      output.list <- list(dev.list, to_keep, fit)
      names(output.list) <- c('deviance_table','variables','jags_model')
      donezo <- 1
      next
    }
    
    #If removing the next predictor brings you to a single predictor remaining, you are done.
    if(ncol(x_mu) == 3){
      #Drop the predictor.
      x_mu <- x_mu[,-which(colnames(x_mu) %in% to.drop)]
      #Grab the model thats the best.
      pos <- which(rownames(dev.list) == to.drop)
      fit <- fit.list[[pos]]
      #wrap the restof the output.
      to_keep <- colnames(x_mu)
      output.list <- list(dev.list, to_keep, fit)
      names(output.list) <- c('deviance_table','variables','jags_model')
      #tell R you are done and report.
      donezo <- 1
      cat('No more predictors to remove. Returning *best* model. Results in output.list. \n')
      next
    }
    
    #If dropping the variables improves deviance by more than 2, awesome. Drop it and repeat the procedure!
    x_mu <- x_mu[,-which(colnames(x_mu) %in% to.drop)]
    #get new reference model to compare things to.
    pos <- which(rownames(dev.list) == to.drop)
    fit <- fit.list[[pos]]
    dev <- fit$deviance
    loop.count <- loop.count + 1
    #tell yourself whats going on.
    cat(to.drop,'removed from model. Refitting...\n')
  } #end while loop
  
  return(output.list)
} #end function
