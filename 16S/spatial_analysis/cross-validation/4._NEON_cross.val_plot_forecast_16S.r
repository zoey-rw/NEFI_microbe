#Forecast using multinomial-dirichlet fit - NEON plot-level 50/50 cross validation.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
#source('NEFI_functions/dmulti_ddirch_forecast.r')
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/dmulti_ddirch_forecast.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

#set output path.----
output.path <- plot.CV_NEON_fcast_16S.path

#load model and NEON site predictors..----
all.mod <- readRDS(plot.CV_NEON_dmulti.ddirch_16S.path)
dat <- readRDS(plot.CV_NEON_cal.val_data_16S.path) #calibration/validation dat core-level NEON.

#Define x_mu and x_sd values.----
plot.preds <- dat$val$x_mu.val
plot.sd    <- dat$val$x_sd.val

#run forecast over all phylo/functional levels.----
all.output <- list()
for(i in 1:length(all.mod)){
  mod <- all.mod[[i]]
  plot.fit <- dmulti_ddirch_forecast(mod, cov_mu = plot.preds, cov_sd = plot.sd, names = plot.preds$plotID)
  #store output as a list and save.----
  output <- list(plot.fit,plot.preds,plot.sd)
  names(output) <- c('plot.fit','plot.preds','plot.sd')
  all.output[[i]] <- output
  cat(names(all.mod)[i],'forecast complete.',i,'of',length(all.mod),'forecasts completed.\n')
}
names(all.output) <- names(all.mod)

#Save output.----
saveRDS(all.output, output.path)
