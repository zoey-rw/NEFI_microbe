#' Title pro_norm_otu.r
#' This proportionally normalizes an OTU table of counts to an OTU table of relative abundances.
#'
#' @param otu #an otu table of counts
#'
#' @return    #an otu table of relative abundances.
#' @export
#'
#' @examples
pro_norm_otu <- function(otu){
  for(i in 1:ncol(otu)){
    otu[,i] <- otu[,i] / sum(otu[,i])
  }
  return(otu)
}
