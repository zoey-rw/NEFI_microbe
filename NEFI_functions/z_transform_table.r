#' z_transform_table()
#' Gotta tell this to only work on numeric columns
#'
#' @param x_mu table of x values
#' @param x_sd associated table of sd values
#'
#' @return a list of length 2. First item is z-transformed x values. second item is z-transformed sd table.
#' @export
#'
#' @examples
source('NEFI_functions/z_transform.r')
z_transform_table <- function(x_mu, x_sd){
  x_mu_nums <- unlist(lapply(x_mu, is.numeric))  
  x_sd_nums <- unlist(lapply(x_sd, is.numeric))  
  sd.out <- data.frame(lapply(x_mu[,x_mu_nums], sd, na.rm= T))
  x_mu[,x_mu_nums] <- data.frame(lapply(x_mu[,x_mu_nums], z_transform))
  #z-transform SD values in precision matrix.
  z_sd <- x_sd[,x_sd_nums]
  for(i in 1:ncol(z_sd)){
    z_sd[,i] <- z_sd[,i] / sd.out[i]
  }
  x_sd[,x_sd_nums] <- z_sd
 return(list(x_mu,x_sd)) 
}



