#' mean_center
#' centers a vector by subtracting the mean.
#' if there no variance in a column, just returns the column untransformed (i.e. a vector of 1s for intercept).
#'
#' @param x matrix or dataframe
#'
#' @return martix or dataframe, mean centered.
#' @export
#'
#' @examples
mean_center <- function(x) {
  ifelse(sd(x) == 0,
         return(x),
         return(x - mean(x, na.rm=T))
  )
}