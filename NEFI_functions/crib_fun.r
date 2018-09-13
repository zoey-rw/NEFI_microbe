#' crib_fun
#' converts a vector of [0,1] values to (0,1) a la Cribari-Neto & Zeileis 2010
#' @param x a vector of values on the interval [0,1]
#' @param N alternative sample size. This is useful when tranforming a matrix in the dirchlet case, rather than just a vector as in the beta case.
#'
#' @return  a vector of values on the interval (0,1)
#' @export
#'
#' @examples
crib_fun <- function(x,N = NA){
  #default use length of vector.
  if( is.na(N)){
    out <- (x * (length(x) - 1) + 0.5) / length(x)  
  }
  #custom- useful when I am in multivariate case.
  if(!is.na(N)){
    out <- (x * (N - 1) + 0.5) / N
  }
  return(out)
}