#Model comparison and variable selection.
#We want to carve down the number of predictors we use to avoid overfitting our prior.
#To do this we propose a few model structures.
#1. intercept only.
#2. linear combination of all predictors.
#3. 'Expert' opinion (Colin's opinion...). Incorporating interactions based on biological hypotheses. We may test multiple models like this.
#We then perform backward selection on each of these model structures (well, nothing to select in the intercept only case.)
#Within a model structure the bakwards variable selection proceeds like this:
#Step 1. calculate model deviance.
#Step 2. Rerun the model without a particular predictor, calculate deviance. Do this for all predictors.
#Step 3. Drop the predictor that results in the largest improvement in deviance (apparently this is the greedy part of a greedy algorithm.)
#Step 4. Using our new best model, repeat Steps 2 and 3 until either:
#         (a.) deviance no longer improves, or 
#         (b.) deviance does not improve beyond a threshold change (say 2 units).


#testing with pseudo data.
rm(list = ls())
library(runjags)
library(doParallel)
source('NEFI_functions/ddirch_site.level_JAGS.r')
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/tic_toc.r')

#### generate pseudo data ####
n.obs <- 60
x1 <- runif(n.obs,20,40)
x2 <- runif(n.obs, 2, 4)
x3 <- runif(n.obs,100,110)
intercept <- rep(1,n.obs)
y <- x1*0.5 + x2*3 + rnorm(n.obs)
#get two species abundances
y1 <- y
y2 <- rep(10, length(y1))
spp.y <- data.frame(cbind(y1,y2))
spp.y <- spp.y / rowSums(spp.y)
x_mu <- data.frame(intercept,x1,x2,x3)

#### Fit first JAGS model to all covariates. ####
#register parallel cores
registerDoParallel(cores = detectCores())

#fit first model, calculate deviance.
tic()
fit <- site.level_dirlichet_jags(y = spp.y, x_mu = x_mu, adapt = 100, burnin = 100, sample = 200)
dev <- fit$deviance
donezo <- 0
loop.count <- 1
cat('Initial model fit complete.')
toc()
while(donezo < 1){
  #sequentially leave out one predictor, except intercept.
  tic()
  fit.list <- list() 
  fit.list <-
  foreach(i = 2:ncol(x_mu)) %dopar% {
    x_new <- x_mu[,-i]
    mod.fit <- site.level_dirlichet_jags(y=spp.y, x_mu = x_new, adapt = 100, burnin = 100, sample = 200)
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
}
