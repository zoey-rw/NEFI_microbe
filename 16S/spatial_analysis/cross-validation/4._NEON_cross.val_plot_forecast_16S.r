#Forecast using multinomial-dirichlet fit - NEON plot-level 50/50 cross validation.
#clear environment, source paths, packages and functions.
rm(list=ls())
library(data.table)
source('paths.r')
source('paths_fall2019.r')
source('NEFI_functions/ddirch_forecast_noLogMap.r')

#set output path.----
output.path <- plot.CV_NEON_fcast_16S.path

#load old forecasts to use as backup
old_fcast <- readRDS(paste0(scc_gen_16S_dir, "JAGS_output/cross_val/plot.CV_NEON_fcast_16S_Sep19.rds"))

#load model and NEON site predictors..----
all.mod <- readRDS(plot.CV_NEON_ddirch_16S_JAGSfit)
dat <- readRDS(plot.CV_NEON_cal.val_data_16S.path) #calibration/validation dat core-level NEON.

#Define x_mu and x_sd values.----
plot.preds <- dat$val$x_mu.val
plot.sd    <- dat$val$x_sd.val

#run forecast over all phylo/functional levels.----
all.output <- list()
for(i in 1:length(all.mod)){
  if (length(all.mod[[i]]) < 8) { # only selects the models that didn't converge
    all.output[[i]] <- old_fcast[[i]] # and uses older forecasts for those models.
    for(m in seq_along(all.output[[i]]$plot.fit)) { # make sure group names are lowercase
      colnames(all.output[[i]]$plot.fit[[m]]) <- tolower(colnames(all.output[[i]]$plot.fit[[m]]))
      }
  } else {
  mod <- all.mod[[i]]
  plot.fit <- ddirch_forecast_noLogMap(mod, cov_mu = plot.preds, cov_sd = plot.sd, names = plot.preds$plotID)
  #store output as a list and save.----
  output <- list(plot.fit,plot.preds,plot.sd)
  names(output) <- c('plot.fit','plot.preds','plot.sd')
  all.output[[i]] <- output
  cat(names(all.mod)[i],'forecast complete.',i,'of',length(all.mod),'forecasts completed.\n')
  }
}
names(all.output) <- names(all.mod)

#Save output.----
saveRDS(all.output, output.path)
