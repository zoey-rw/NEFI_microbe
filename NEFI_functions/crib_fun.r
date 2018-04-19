#converts a vector of [0,1] values to (0,1) a la Cribari-Neto & Zeileis 2010

crib_fun <- function(x){(x * (length(x) - 1) + 0.5) / length(x)}