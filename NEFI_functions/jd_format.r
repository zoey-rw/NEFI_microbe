#This function formats a dataframe into a jags data object.
#N in your model must be the number of observations.

jd_format <- function(y,preds){
  #make sure y is not a data.table
  if(is.data.table(y) == T){ y <- as.data.frame(y)}
  #subset to complete cases.
  all <- cbind(y,preds)
  all <- all[complete.cases(all),]
  y <- all[,1]
  y <- as.vector(y)
  preds <- all[,-1]
  #assign fixed effects to objects x1 -> xn
  colnames(preds)  <- paste0('x', 1:(ncol(preds)))
  jd <- as.list((preds))
  
  #Get jags data object together.
  pre.list <- list(y=y, N = length(y))
  jd <- c(pre.list, as.list(preds))
  return(jd)
}