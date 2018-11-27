#Making a spatial forecast based on the prior data to NEON sites at core, plot and site levels.
#downstream this will log transform map, which prevents this from generalizing beyond the dirichlet example.
#This script depends on the following packages: DirichletReg.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/precision_matrix_match.r')
source('NEFI_functions/ddirch_forecast.r')
library(data.table)

library(RCurl)
# source paths.r from colins github
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/paths.r", ssl.verifypeer = FALSE)
eval(parse(text = script))
# source ddirch_forecast
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/ddirch_forecast.r", ssl.verifypeer = FALSE)
eval(parse(text = script))


#set output path.----
output.path <- NEON_site_fcast_genera_16S.path

#load model results.----
#mod 1 is data from maps.
#mod 2 is site-specific data, no maps.
#mod 3 is all covariates.
mod <- readRDS(bahram_16S.prior_12gen_JAGSfit_no_moisture)
mod <- mod[[3]] #just the all predictor case.

#get core-level covariate means and sd.----
dat <- readRDS(hierarch_filled_16S.path)
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

# change predictor from 'pH' to 'PH' to match prior - fix upstream if prior model is re-run
setnames(core.preds, "pH", "PH")
setnames(core.sd, "pH", "PH")
setnames(plot.preds, "pH", "PH")
setnames(plot.sd, "pH", "PH")
setnames(site.preds, "pH", "PH")
setnames(site.sd, "pH", "PH")

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




### plot observed data onto forecast ###

# load forecast 
output <- readRDS(NEON_site_fcast_genera_16S.path)
core.fit <- output$core.fit

# read in obs table to get deprecatedVialID
map<-readRDS("/fs/data3/caverill/NEFI_data/16S/pecan_gen/NEON_data_aggregation/obs.table_16S.rds")
#map$deprecatedVialID

# read in relative abundances of NEON genera
truth <- readRDS(NEON_cosmo_abundances.path) 
truth <- truth$rel.abundances

# check what fraction of rows are empty for each genus.
colSums(truth==0, na.rm=T)/colSums(!is.na(truth), na.rm = T)

truth$deprecatedVialID <- rownames(truth)
truth1 <- merge(truth, map, by = "deprecatedVialID")
truth1 <- truth1[!duplicated(truth1$geneticSampleID),]
rownames(truth1) <- gsub('-GEN','',truth1$geneticSampleID)


#validate against observed data by plotting.----
par(mfrow = c(2,2))
trans <- 0.3
limy <- c(0,.1)
i = 2 

#core.level.----
pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/NEON_core.fcast_cosmo_16S.pdf")
#organize data.
fcast <- core.fit
truth <- readRDS(NEON_cosmo_abundances.path) 
truth <- truth$rel.abundances
truth$deprecatedVialID <- rownames(truth)
truth1 <- merge(truth, map[,c("deprecatedVialID", "geneticSampleID")], by = "deprecatedVialID")
truth1 <- truth1[!duplicated(truth1$geneticSampleID),]
rownames(truth1) <- gsub('-GEN','',truth1$geneticSampleID)
truth <- truth1
truth <- truth[rownames(truth) %in% rownames(fcast$mean),]
test <- truth[-grep('DSNY',rownames(truth)),]

for (i in 1:20) {

for(k in 1:length(fcast)){
  fcast[[k]] <- fcast[[k]][rownames(fcast[[k]]) %in% rownames(truth),]
}
truth <- as.matrix(truth)
mu <- fcast$mean[,i][order(fcast$mean[,i])]
ci_0.975 <- fcast$ci_0.975[,i][order(match(names(fcast$ci_0.975[,i]),names(mu)))]
ci_0.025 <- fcast$ci_0.025[,i][order(match(names(fcast$ci_0.025[,i]),names(mu)))]
pi_0.975 <- fcast$pi_0.975[,i][order(match(names(fcast$pi_0.975[,i]),names(mu)))]
pi_0.025 <- fcast$pi_0.025[,i][order(match(names(fcast$pi_0.025[,i]),names(mu)))]
bac_name <- colnames(fcast$mean)[i]
obs.mu   <- truth[,c(bac_name)][order(match(names(truth[,c(bac_name)]),names(mu)))]
#obs.lo95 <- fglo95[,i][order(match(names(fglo95[,i]),names(mu)))]
#obs.hi95 <- fghi95[,i][order(match(names(fghi95[,i]),names(mu)))]
#plot
plot(obs.mu ~ mu, cex = 0.7, ylim=limy, main = paste0('core-level ', bac_name))
rsq <- round(summary(lm(obs.mu ~mu))$r.squared,2)
mtext(paste0('R2=',rsq), side = 3, cex = 0.7, line = -2.7, adj = 0.03)
#add confidence interval.
range <- mu
polygon(c(range, rev(range)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
polygon(c(range, rev(range)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
#fraction within 95% predictive interval.
in_it <- round(sum(as.numeric(obs.mu) < pi_0.975 & as.numeric(obs.mu) > pi_0.025, na.rm = TRUE) / length(obs.mu),2) * 100
state <- paste0(in_it,'% of observations within 95% prediction interval.')
mtext(state,side = 3, cex = 0.7, line = -1.3, adj = 0.05)
abline(0,1,lwd=2)
}
dev.off()



# plot.level
pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/NEON_plot.fcast_cosmo_16S.pdf")

fcast <- output$plot.fit	
truth <- readRDS(NEON_plot.level_genera_obs_16S.path)	

for (i in 1:20) {
  
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
bac_name <- colnames(fcast$mean)[i]	
obs.mu   <- truth$mean[,bac_name][order(match(names(truth$mean[,bac_name]),names(mu)))]	
#obs.lo95 <- truth$lo95[,bac_name][order(match(names(truth$lo95[,bac_name]),names(mu)))]	
#obs.hi95 <- truth$hi95[,bac_name][order(match(names(truth$hi95[,bac_name]),names(mu)))]	
#plot	
plot(obs.mu ~ mu, cex = 0.7, ylim=limy, main = paste0('plot-level ', bac_name))	
arrows(c(mu), obs.lo95, c(mu), obs.hi95, length=0.05, angle=90, code=3)	
rsq <- round(summary(lm(obs.mu ~mu))$r.squared,2)	
mtext(paste0('R2=',rsq), side = 3, cex = .7, line = -2.7, adj = 0.03)	
#1:1 line	
abline(0,1, lwd = 2)	
#add confidence interval.	
range <- mu	
polygon(c(range, rev(range)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)	
polygon(c(range, rev(range)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)	
#fraction within 95% predictive interval.	
in_it <- round(sum(obs.mu < pi_0.975 & obs.mu > pi_0.025) / length(obs.mu),2) *100	
state <- paste0(in_it,'% of observations within 95% prediction interval.')	
mtext(state,side = 3, cex = .7, line = -1.3, adj = 0.05)
}
dev.off()





#site.level----
pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/NEON_site.fcast_cosmo_16S.pdf")
#organize data.
fcast <- output$site.fit
truth <- readRDS(NEON_site.level_genera_obs_16S.path)

for(i in 1:20) {
  
for(k in 1:length(truth)){
  rownames(truth[[k]]) <- gsub('.','_',rownames(truth[[k]]), fixed = T)
  truth[[k]] <- truth[[k]][rownames(truth[[k]]) %in% rownames(fcast$mean),]
}
#drop DSNY
for(k in 1:length(truth)){
  truth[[k]] <- truth[[k]][-(grep('DSNY',rownames(truth[[k]]))),]
}
for(k in 1:length(fcast)){
  fcast[[k]] <- fcast[[k]][rownames(fcast[[k]]) %in% rownames(truth$mean),]
}
mu <- fcast$mean[,i][order(fcast$mean[,i])]
ci_0.975 <- fcast$ci_0.975[,i][order(match(names(fcast$ci_0.975[,i]),names(mu)))]
ci_0.025 <- fcast$ci_0.025[,i][order(match(names(fcast$ci_0.025[,i]),names(mu)))]
pi_0.975 <- fcast$pi_0.975[,i][order(match(names(fcast$pi_0.975[,i]),names(mu)))]
pi_0.025 <- fcast$pi_0.025[,i][order(match(names(fcast$pi_0.025[,i]),names(mu)))]
bac_name <- colnames(fcast$mean)[i]
obs.mu   <- truth$mean[,bac_name][order(match(names(truth$mean[,bac_name]),names(mu)))]
obs.lo95 <- truth$lo95[,bac_name][order(match(names(truth$lo95[,bac_name]),names(mu)))]
obs.hi95 <- truth$hi95[,bac_name][order(match(names(truth$hi95[,bac_name]),names(mu)))]
#plot
plot(obs.mu ~ mu, cex = 0.7, ylim=limy, main=paste0('site-level ', bac_name))
arrows(c(mu), obs.lo95, c(mu), obs.hi95, length=0.05, angle=90, code=3)
rsq <- round(summary(lm(obs.mu ~mu))$r.squared,2)
mtext(paste0('R2=',rsq), cex = .7, side = 3, line = -2.7, adj = 0.03)
#1:1 line
abline(0,1, lwd = 2)
#add confidence interval.
range <- mu
polygon(c(range, rev(range)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
polygon(c(range, rev(range)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
#fraction within 95% predictive interval.
in_it <- round(sum(obs.mu < pi_0.975 & obs.mu > pi_0.025) / length(obs.mu),2) *100
state <- paste0(in_it,'% of observations within 95% prediction interval.')
mtext(state,side = 3,  cex = .7, line = -1.3, adj = 0.05)
}
dev.off()
