#' sd_to_precision()
#'only converts numeric columns.
#'
#' @param x_sd table of standard deviation values.
#'
#' @return return a matrix or dataframe of precision values.
#' @export
#'
#' @examples
sd_to_precision <- function(x_sd){
  num <- unlist(lapply(x_sd, is.numeric))  
  z_sd <- data.frame(x_sd[,num])
  for(i in 1:ncol(z_sd)){
      z_sd[,i] <- 1/(z_sd[,i]^2)
  }
  x_sd[,num] <- z_sd
  x_sd <- as.matrix(x_sd) #JAGS likes things served as matrices.
  return(x_sd)
}