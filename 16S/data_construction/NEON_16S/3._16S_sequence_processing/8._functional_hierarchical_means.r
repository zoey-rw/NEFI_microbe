# Get plot- and site-level hierarchical means of functional groups from rel abundances at NEON cores

rm(list=ls())
library(runjags)
source('paths.r')
#source('NEFI_functions/hierarch_ddirch_means.r')

# source hierarch means function
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/hierarch_ddirch_means.r", ssl.verifypeer = FALSE)
eval(parse(text = script))
# colin's paths.R
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/paths.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

