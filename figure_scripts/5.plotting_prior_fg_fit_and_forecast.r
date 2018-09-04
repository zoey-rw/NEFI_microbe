#plotting all FG with same forecast_plot.
rm(list=ls())
source('paths.r')
source('NEFI_functions/ddirch_forecast_site.level.r')
source('NEFI_functions/forecast_plot.r')

#output path
figure.path <- 'figures/fg_prior_prediction.png'

#data prep.----
#load models and tedersoo prior data.
mod <- readRDS(ted_ITS.prior_fg_JAGSfit)
d <- readRDS(tedersoo_ITS.prior_for_analysis.path)

#get covariates as complete cases.
preds <- as.character(mod[[3]]$species_parameter_output[[1]]$predictor)
d.covs <- d[,colnames(d) %in% preds]
d <- d[complete.cases(d.covs),]
d.covs <- d.covs[complete.cases(d.covs),]
#always log map.
d.covs$map <- log(d.covs$map)

#isolate observed functional groups.
obs <- d[,colnames(d) %in% c('Ectomycorrhizal','Arbuscular','Saprotroph','Pathogen')]

#get forecast
f.cast <- ddirch_forecast_site.level(mod = mod,site_covs = d.covs)



#getting plotting dataframes together.----
plot.list <- list()
for(i in 1:ncol(obs)){
  of_interest <- colnames(obs)[i]
  for(k in 1:length(f.cast$all.preds)){
    to_add <- f.cast$all.preds[[k]][,of_interest]
    if(k == 1){out <- to_add}
    if(k  > 1){out <- cbind(out,to_add)}
  }
  colnames(out) <- names(f.cast$all.preds)
        to.plot <- data.frame(obs[,colnames(obs) %in% of_interest],out)
        to.plot <- to.plot[order(to.plot$mean),]
        to.plot <- to.plot*100
        colnames(to.plot)[1] <- 'obs'
        plot.list[[i]] <- to.plot
}
names(plot.list) <- colnames(obs)

#begin plotting.----
png(filename=figure.path,width=8.5,height=8,units='in',res=300)
par(mfrow = c(2,2))
for(i in 1:length(plot.list)){
  forecast_plot(plot.list[[i]],
                ylab = paste0('% ',names(plot.list)[i],' Fungi'))
}
dev.off()