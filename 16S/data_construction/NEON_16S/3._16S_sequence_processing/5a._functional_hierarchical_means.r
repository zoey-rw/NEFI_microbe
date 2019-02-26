# Get plot- and site-level hierarchical means of functional groups from rel abundances at NEON cores

rm(list=ls())
library(runjags)
library(foreach)
library(doParallel)
library(parallel)
source('paths.r')
source('NEFI_functions/tic_toc.r')
#source('NEFI_functions/hierarch_ddirch_means.r')

# source hierarch means function
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/hierarch_ddirch_means.r", ssl.verifypeer = FALSE)
eval(parse(text = script))
# colin's paths.R
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/paths.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

#set output.path----
output.path <- NEON_all.fg_plot.site_obs_16S.path

#load data and format.----
d <- readRDS(prior_fg_abundances_16S.path)

#register parallel environment.----
n.cores <- detectCores()
registerDoParallel(n.cores)

tic()
# loop over models
output <- list()
output <-
  foreach(i = 1:length(d)) %dopar% { # loop through each functional group
    
    y <- d[[i]]

    #Get y multivariate matrix.
    abundances <- y$abundances
    y <- as.matrix((abundances + 1) / rowSums(abundances + 1))
    
    #get core_plot and plot_site indexing.
    core_plot <- substr(rownames(y), 1, 8)
    core_site <- substr(rownames(y), 1, 4)
    plot_site <- unique(core_plot)
    plot_site <- substr(plot_site, 1, 4)
    
    #fit the hierarchical means.
    fit <- hierarch_ddirch_means(y=y, core_plot = core_plot, plot_site = plot_site, jags.method = 'parallel')
    
    #add row and column names - plot level (you should put this in the function).
    for(j in 1:length(fit$plot.fit)){
      rownames(fit$plot.fit[[j]]) <- unique(core_plot)
      rownames(fit$site.fit[[j]]) <- unique(plot_site)
      colnames(fit$plot.fit[[j]]) <- colnames(y)
      colnames(fit$site.fit[[j]]) <- colnames(y)
    }
    fit$core.fit <- y #add in the core-level data!
    cat(paste0("Models fit for set ",i,"/12. \n"))
    return(fit)
  }
names(output) <- names(d)

#save matrix lists.----
saveRDS(output, output.path)
cat('Script complete. ');toc()