#' cn_uncertainty_neon.r
#' Calculates observation sd of C:N based on analytical replicates of NEON data.
#' sd scales with magnitude of C:N observation (heteroscedastic).
#'
#' @param x value or vector of C:N observations to estimate observation uncertainty of.
#'
#' @return vector of observation sd.
#' @export
#'
#' @examples
pC_uncertainty_neon.r <- function(x){
  #load model fit.
  source('paths.r')
  mod <- readRDS(cn_NEON_uncertainty_model.path)
  #grab k parameters that predict standard deviation.
  pars <- mod[grep('k',rownames(mod)),][,4]
  sd <- pars[1] + x*pars[2]
  #return standard deviation.
  return(sd)
}