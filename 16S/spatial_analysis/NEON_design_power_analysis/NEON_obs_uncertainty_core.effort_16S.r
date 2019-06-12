#Simulating the effect of observation uncertainty on overall predictability at core, plot and site scales.
#manipulating number of cores per plot.
#clear environment, load packages, source functions.
rm(list=ls())
library(DirichletReg)
source('paths.r')
#source('NEFI_functions/ddirch_obs_uncertainty.r')
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/ddirch_obs_uncertainty.r", ssl.verifypeer = FALSE)
eval(parse(text = script))
source('NEFI_functions/tic_toc.r')

#set output path.----
output.path <- NEON_observation_uncertainty_core_16S.path

#SEt number of cores per plot.5 levels.
cores.per.plot <- c(3,6,9,12,15)

#Generate site values, set parameters.----
set.seed(1234)
n.site <- 13
y1 <- runif(n.site, 0.1, 0.9)    #abundant taxon, varies alot.
y2 <- runif(n.site, 0.001, 0.01) #rare taxon.
y3 <- runif(n.site, 0.1, 0.2)    #abudant taxon, varies comparatively less.
y4 <- 1 - (y1 + y2 + y3)              #'other' so we sum to 1.
site_mu <- as.matrix(data.frame(y1,y2,y3,y4))

#assign variances, higher number means lower variance.
#These are based on Harvard Forest and Dopheide 2018 intra-core uncertainty.
core.var <- 300 
plot.var <- 200
site.var <- 60

#run the simulation.----
#This takes a real long time. Function automatically detects number of available processors.
cat('Running simulation...\n');tic()
plot_output <- list()
for(i in 1:length(cores.per.plot)){
  run <- ddirch_obs_uncertainty(site_mu, site.var = site.var, plot.var = plot.var, core.var = core.var, 
                                n.sim = 100, n.core = cores.per.plot[i])
  plot_output[[i]] <- run
  cat(i,'of',length(cores.per.plot),'simulations complete. ');toc()
}
names(plot_output) <- cores.per.plot

cat('All simulations complete.\n');toc()

#Save output.----
saveRDS(plot_output, output.path)

