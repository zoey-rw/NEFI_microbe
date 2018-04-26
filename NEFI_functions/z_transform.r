#' z_transform
#' z transforms all columns of a matrix or data farme by subtracting each vectors mean, dividing by sd.
#' if there no variance in a column, just returns the column untransformed.
#'
#' @param x matrix or dataframe
#'
#' @return martix or dataframe, z transformed.
#' @export
#'
#' @examples
z_transform <- function(x) {
  ifelse(sd(x) == 0,
         return(x),
         return((x - mean(x, na.rm=T)) / sd(x, na.rm=T))
  )
}