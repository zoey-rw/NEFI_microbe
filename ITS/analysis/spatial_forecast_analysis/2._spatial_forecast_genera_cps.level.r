#Making a spatial forecast based on the prior data to NEON sites at core, plot and site levels.
#downstream this will log transform map, which prevents this from generalizing beyond the dirichlet example.
#This script depends on the following packages: DirichletReg.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/precision_matrix_match.r')
source('NEFI_functions/ddirch_forecast.r')

#set output path.----
output.path <- NEON_site_fcast_genera.path

#load model results.----
#mod 1 is data from maps.
#mod 2 is site-specific data, no maps.
#mod 3 is all covariates.
mod <- readRDS(ted_ITS.prior_20gen_JAGSfit)
mod <- mod[[3]] #just the all predictor case.

#get core-level covariate means and sd.----
dat <- readRDS(hierarch_filled.path)
core_mu <- dat$core.core.mu
plot_mu <- dat$plot.plot.mu
site_mu <- dat$site.site.mu
#merge together.
plot_mu$siteID <- NULL
core.preds <- merge(core_mu   , plot_mu)
core.preds <- merge(core.preds, site_mu)
core.preds$relEM <- NULL
names(core.preds)[names(core.preds)=="b.relEM"] <- "relEM"

#get core-level SD.
core_sd <- dat$core.core.sd
plot_sd <- dat$plot.plot.sd
site_sd <- dat$site.site.sd
#merge together.
plot_sd$siteID <- NULL
core.sd <- merge(core_sd   , plot_sd)
core.sd <- merge(core.sd, site_sd)
core.sd$relEM <- NULL
names(core.sd)[names(core.sd)=="b.relEM"] <- "relEM"


#get plot-level covariate means and sd.----
core_mu <- dat$core.plot.mu
plot_mu <- dat$plot.plot.mu
site_mu <- dat$site.site.mu
#merge together, .
plot_mu$siteID <- NULL
plot.preds <- merge(core_mu,plot_mu)
plot.preds <- merge(plot.preds,site_mu)
plot.preds$relEM <- NULL
names(plot.preds)[names(plot.preds)=="b.relEM"] <- "relEM"

#get plot-level SD.
core_sd <- dat$core.plot.sd
plot_sd <- dat$plot.plot.sd
site_sd <- dat$site.site.sd
#merge together
plot_sd$siteID <- NULL
plot.sd <- merge(core_sd,plot_sd)
plot.sd <- merge(plot.sd,site_sd)
plot.sd$relEM <- NULL
names(plot.sd)[names(plot.sd)=='b.relEM'] <- "relEM"

#get site-level covariate means and sd.----
core_mu <- dat$core.site.mu
plot_mu <- dat$plot.site.mu
site_mu <- dat$site.site.mu
#merge together
site.preds <- merge(core_mu, plot_mu)
site.preds <- merge(site.preds,site_mu)
names(site.preds)[names(site.preds)=='b.relEM'] <- "relEM"

#get site-level SD.
core_sd <- dat$core.site.sd
plot_sd <- dat$plot.site.sd
site_sd <- dat$site.site.sd
#merge together.
site.sd <- merge(core_sd,plot_sd)
site.sd <- merge(site.sd,site_sd)
names(site.sd)[names(site.sd)=='b.relEM'] <- "relEM"

#Get forecasts from ddirch_forecast.----
core.fit <- ddirch_forecast(mod=mod, cov_mu=core.preds, cov_sd=core.sd, names=core.preds$sampleID)
plot.fit <- ddirch_forecast(mod=mod, cov_mu=plot.preds, cov_sd=plot.sd, names=plot.preds$plotID)
site.fit <- ddirch_forecast(mod=mod, cov_mu=site.preds, cov_sd=site.sd, names=site.preds$siteID)

#store output as a list and save.----
output <- list(core.fit,plot.fit,site.fit,core.preds,plot.preds,site.preds,core.sd,plot.sd,site.sd)
names(output) <- c('core.fit','plot.fit','site.fit',
                   'core.preds','plot.preds','site.preds',
                   'core.sd','plot.sd','site.sd')
saveRDS(output, output.path)

#validate against observed data by plotting.----
par(mfrow = c(3,1))
trans <- 0.3
limy <- c(0,1)
i = 2 #ECM fungi

#core.level.----
#organize data.
fcast <- core.fit
truth <- readRDS(NEON_cosmo_genera.path) #UPDATE THIS FOR GENERA.
truth <- truth$rel.abundances
rownames(truth) <- gsub('-GEN','',truth$geneticSampleID)
truth <- truth[rownames(truth) %in% rownames(core.fit$mean),]
for(k in 1:length(fcast)){
  fcast[[k]] <- fcast[[k]][rownames(fcast[[k]]) %in% rownames(truth),]
}
truth <- as.matrix(truth)
mu <- fcast$mean[,i][order(fcast$mean[,i])]
ci_0.975 <- fcast$ci_0.975[,i][order(match(names(fcast$ci_0.975[,i]),names(mu)))]
ci_0.025 <- fcast$ci_0.025[,i][order(match(names(fcast$ci_0.025[,i]),names(mu)))]
pi_0.975 <- fcast$pi_0.975[,i][order(match(names(fcast$pi_0.975[,i]),names(mu)))]
pi_0.025 <- fcast$pi_0.025[,i][order(match(names(fcast$pi_0.025[,i]),names(mu)))]
fungi_name <- colnames(fcast$mean)[i]
obs.mu   <- truth[,i][order(match(names(truth[,c(fungi_name)]),names(mu)))]
#obs.lo95 <- fg$lo95[,i][order(match(names(fg$lo95[,i]),names(mu)))]
#obs.hi95 <- fg$hi95[,i][order(match(names(fg$hi95[,i]),names(mu)))]
#plot
plot(obs.mu ~ mu, cex = 0.7, ylim=limy, main = 'core-level')
rsq <- round(summary(lm(obs.mu ~mu))$r.squared,2)
mtext(paste0('R2=',rsq), side = 3, line = -2.7, adj = 0.03)
#add confidence interval.
range <- mu
polygon(c(range, rev(range)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
polygon(c(range, rev(range)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
#fraction within 95% predictive interval.
in_it <- round(sum(obs.mu < pi_0.975 & obs.mu > pi_0.025) / length(obs.mu),2) * 100
state <- paste0(in_it,'% of observations within 95% prediction interval.')
mtext(state,side = 3, line = -1.3, adj = 0.05)
abline(0,1,lwd=2)


#plot.level----
#organize data.
fcast <- plot.fit
truth <- readRDS(NEON_plot.level_fg_obs.path)
for(k in 1:length(truth)){
  rownames(truth[[k]]) <- gsub('.','_',rownames(truth[[k]]), fixed = T)
  truth[[k]] <- truth[[k]][rownames(truth[[k]]) %in% rownames(fcast$mean),]
}
for(k in 1:length(fcast)){
  fcast[[k]] <- fcast[[k]][rownames(fcast[[k]]) %in% rownames(truth$mean),]
}
mu <- fcast$mean[,i][order(fcast$mean[,i])]
ci_0.975 <- fcast$ci_0.975[,i][order(match(names(fcast$ci_0.975[,i]),names(mu)))]
ci_0.025 <- fcast$ci_0.025[,i][order(match(names(fcast$ci_0.025[,i]),names(mu)))]
pi_0.975 <- fcast$pi_0.975[,i][order(match(names(fcast$pi_0.975[,i]),names(mu)))]
pi_0.025 <- fcast$pi_0.025[,i][order(match(names(fcast$pi_0.025[,i]),names(mu)))]
fungi_name <- colnames(fcast$mean)[i]
obs.mu   <- truth$mean[,fungi_name][order(match(names(truth$mean[,fungi_name]),names(mu)))]
obs.lo95 <- truth$lo95[,fungi_name][order(match(names(truth$lo95[,fungi_name]),names(mu)))]
obs.hi95 <- truth$hi95[,fungi_name][order(match(names(truth$hi95[,fungi_name]),names(mu)))]
#plot
plot(obs.mu ~ mu, cex = 0.7, ylim=limy, main = 'plot-level')
arrows(c(mu), obs.lo95, c(mu), obs.hi95, length=0.05, angle=90, code=3)
rsq <- round(summary(lm(obs.mu ~mu))$r.squared,2)
mtext(paste0('R2=',rsq), side = 3, line = -2.7, adj = 0.03)
#1:1 line
abline(0,1, lwd = 2)
#add confidence interval.
range <- mu
polygon(c(range, rev(range)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
polygon(c(range, rev(range)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
#fraction within 95% predictive interval.
in_it <- round(sum(obs.mu < pi_0.975 & obs.mu > pi_0.025) / length(obs.mu),2) *100
state <- paste0(in_it,'% of observations within 95% prediction interval.')
mtext(state,side = 3, line = -1.3, adj = 0.05)


#site.level----
#organize data.
fcast <- site.fit
truth <- readRDS(NEON_site.level_fg_obs.path)
for(k in 1:length(truth)){
  rownames(truth[[k]]) <- gsub('.','_',rownames(truth[[k]]), fixed = T)
  truth[[k]] <- truth[[k]][rownames(truth[[k]]) %in% rownames(fcast$mean),]
}
for(k in 1:length(fcast)){
  fcast[[k]] <- fcast[[k]][rownames(fcast[[k]]) %in% rownames(truth$mean),]
}
mu <- fcast$mean[,i][order(fcast$mean[,i])]
ci_0.975 <- fcast$ci_0.975[,i][order(match(names(fcast$ci_0.975[,i]),names(mu)))]
ci_0.025 <- fcast$ci_0.025[,i][order(match(names(fcast$ci_0.025[,i]),names(mu)))]
pi_0.975 <- fcast$pi_0.975[,i][order(match(names(fcast$pi_0.975[,i]),names(mu)))]
pi_0.025 <- fcast$pi_0.025[,i][order(match(names(fcast$pi_0.025[,i]),names(mu)))]
fungi_name <- colnames(fcast$mean)[i]
obs.mu   <- truth$mean[,fungi_name][order(match(names(truth$mean[,fungi_name]),names(mu)))]
obs.lo95 <- truth$lo95[,fungi_name][order(match(names(truth$lo95[,fungi_name]),names(mu)))]
obs.hi95 <- truth$hi95[,fungi_name][order(match(names(truth$hi95[,fungi_name]),names(mu)))]
#plot
plot(obs.mu ~ mu, cex = 0.7, ylim=limy, main='site-level')
arrows(c(mu), obs.lo95, c(mu), obs.hi95, length=0.05, angle=90, code=3)
rsq <- round(summary(lm(obs.mu ~mu))$r.squared,2)
mtext(paste0('R2=',rsq), side = 3, line = -2.7, adj = 0.03)
#1:1 line
abline(0,1, lwd = 2)
#add confidence interval.
range <- mu
polygon(c(range, rev(range)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
polygon(c(range, rev(range)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
#fraction within 95% predictive interval.
in_it <- round(sum(obs.mu < pi_0.975 & obs.mu > pi_0.025) / length(obs.mu),2) *100
state <- paste0(in_it,'% of observations within 95% prediction interval.')
mtext(state,side = 3, line = -1.3, adj = 0.05)

