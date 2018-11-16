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
  to_flip <- list()
  for(i in 1:ncol(x_sd)){
    to_flip[[i]] <- is.numeric(x_sd[,i])
  }
  to_flip <- unlist(to_flip)
  z_sd <- data.frame(x_sd[,to_flip])
  for(i in 1:ncol(z_sd)){
      z_sd[,i] <- 1/(z_sd[,i]^2)
  }
  #x_sd[,to_flip] <- z_sd
  for(i in 1:ncol(x_sd)){
    if(to_flip[i] == T){
      x_sd[,i] <- z_sd[,i]
    }
  }
  x_sd <- as.matrix(x_sd) #JAGS likes things served as matrices.
  return(x_sd)
}