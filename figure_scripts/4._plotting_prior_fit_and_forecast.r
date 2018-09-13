#Generating predicted vs. observed data for prior fungal w/ 95% credible interval and 95% predictive interval.
#STart with just Ectomycrrhizal fungi, then generalize to all functional groups.
#clear environment, load packages, functions and paths.
rm(list=ls())
source('paths.r')
source('NEFI_functions/ddirch_forecast_site.level.r')
source('NEFI_functions/forecast_plot.r')

#output path
figure.path <- ecto_prior_prediction.fig.path

#data prep for figure.----
#load model results.
#mod 1 is data from maps.
#mod 2 is site-specific data, no maps.
#mod 3 is all covariates.
mod <- readRDS(ted_ITS.prior_fg_JAGSfit)

#load tedersoo prior data.
d <- readRDS(tedersoo_ITS.prior_for_analysis.path)

#get covariates as complete cases.
preds <- as.character(mod[[3]]$species_parameter_output[[1]]$predictor)
d.covs <- d[,colnames(d) %in% preds]
d <- d[complete.cases(d.covs),]
d.covs <- d.covs[complete.cases(d.covs),]
#always log map.
d.covs$map <- log(d.covs$map)

#Use forecast function to sample from parameter space.
f.cast <- ddirch_forecast_site.level(mod = mod,site_covs = d.covs, n.samp = 2000)

#isolate observed.
obs <- d$Ectomycorrhizal
#get obs, predicted and intervals.
of_interest <- 'Ectomycorrhizal'
for(i in 1:length(f.cast$all.preds)){
  to_add <- f.cast$all.preds[[i]][,of_interest]
  if(i == 1){out <- to_add}
  if(i  > 1){out <- cbind(out,to_add)}
}
colnames(out) <- names(f.cast$all.preds)

#bind together observed and predicted, order by predicted mean.
to.plot <- data.frame(obs,out)
to.plot <- to.plot[order(to.plot$mean),]
to.plot <- to.plot*100

#how many of your observations fall w/ in the 95% predictive interval?

#begin plot.----
png(filename=figure.path,width=6,height=6,units='in',res=300)
forecast_plot(to.plot,ylab = '% Ectomycorrhizal Fungi')
#end plot.----
dev.off()