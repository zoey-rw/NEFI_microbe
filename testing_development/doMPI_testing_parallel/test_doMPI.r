#clear environment, load packages.
rm(list=ls())
library(doMPI)
#clock functions
tic = function() {assign("timer", Sys.time(), envir=.GlobalEnv)}
toc = function() print(Sys.time()-timer)

#generate huge vectors to fit a model to.
n.sims <- 1000

#fit model without parallelization.
cat('Fitting model without parallelization...n')
output.list <- list()
tic()
for(i in 1:n.sims){
  set.seed(i)
  #x <- rnorm(1000000)
  #y <- x*0.5 + 1 + rnorm(length(x))
  #mod <- lm(y ~ x)
  #output.list[[i]] <- mod
  mean(rnorm(1e7))
}
toc()

#Use foreach and instructions to fit in parallel.
cat('Fitting model with parallelization...n')
cl <- startMPIcluster(count=30)
registerDoMPI(cl)
tic()
parallel.output <-
  foreach(i = 1:n.sims) %dopar% {
    set.seed(i)
    #x <- rnorm(1000000)
    #y <- x*0.5 + 1 + rnorm(length(x))
    #mod <- lm(y ~ x)
    #return(mod)
    return(mean(rnorm(1e7)))
  }
toc()

#close connection thing.
closeCluster(cl)
mpi.quit()