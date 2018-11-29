# Validation plots of forecasts to NEON core/plot/sites, 
# plotted by genera, each with c/p/s.

rm(list=ls())
source('paths.r')
source('NEFI_functions/crib_fun.r')
library(data.table)


# load forecast 
output <- readRDS(NEON_site_fcast_genera_16S.path)
# read in obs table
map <- readRDS(obs.table_16S.path)

# get rsq from priors
# read in prior fit.
fit <- readRDS(bahram_16S.prior_12gen_JAGSfit)
fit <- fit$all.preds
prior_rsq <- list()
for(i in 1:ncol(fit$predicted)){
  mod <- betareg::betareg(crib_fun(fit$observed[,i]/rowSums(fit$observed)) ~ crib_fun(fit$predicted[,i]))
  rsq <-round(summary(mod)$pseudo.r.squared, 2)
  prior_rsq[i] <- as.character(rsq)
}
prior_rsq <- as.numeric(prior_rsq)
names(prior_rsq) <- colnames(fit$predicted)

# get % of NEON cores a genus is present in.
truth <- readRDS(NEON_cosmo_abundances_16S.path) 
truth <- truth$rel.abundances
present_percent <- colSums(truth!=0, na.rm=T)/colSums(!is.na(truth), na.rm = T)


#validate against observed data by plotting.----
trans <- 0.3
limy <- c(0,.1)
#i = 2 

pdf(NEON_cps.fcast_cosmo_16S.path)
par(mfrow = c(3,1))
par(mar = c(2,2,2,2))
par(oma = c(0,0,2,0))

for (i in 2:21) {
  
#core.level.----
#organize data.
fcast <- output$core.fit
truth <- readRDS(NEON_cosmo_abundances_16S.path) 
truth <- truth$rel.abundances
truth$deprecatedVialID <- rownames(truth)
truth1 <- merge(truth, map[,c("deprecatedVialID", "geneticSampleID")], by = "deprecatedVialID")
truth1 <- truth1[!duplicated(truth1$geneticSampleID),]
rownames(truth1) <- gsub('-GEN','',truth1$geneticSampleID)
truth <- truth1
truth <- truth[rownames(truth) %in% rownames(fcast$mean),]
test <- truth[-grep('DSNY',rownames(truth)),]

  
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
  
  # plot
  # get ylim
  obs_limit <- max(obs.mu, na.rm = T)
  if(max(pi_0.975) > as.numeric(obs_limit)){obs_limit <- max(pi_0.975)}
  limy <- as.numeric(obs_limit)*1.05
  if(limy > 0.95){limy <- 1}
  
  plot(obs.mu ~ mu, cex = 0.7, ylim=c(0,limy), main = paste0('core-level ', bac_name))
  rsq <- round(summary(lm(obs.mu ~mu))$r.squared,2)
  mtext(paste0('R2=',rsq), side = 3, cex = 0.7, line = -2.7, adj = 0.03)
  #1-to-1 line.
  abline(0,1,lwd=2)
  # best-fit line.
  abline(lm(obs.mu ~mu), lty = 2, col = 'purple')
  #add confidence interval.
  range <- mu
  polygon(c(range, rev(range)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
  polygon(c(range, rev(range)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
  #fraction within 95% predictive interval.
  in_it <- round(sum(as.numeric(obs.mu) < pi_0.975 & as.numeric(obs.mu) > pi_0.025, na.rm = TRUE) / length(obs.mu),2) * 100
  state <- paste0(in_it,'% of observations within 95% prediction interval.')
  mtext(state,side = 3, cex = 0.7, line = -1.3, adj = 0.05)



# plot.level
fcast <- output$plot.fit	
truth <- readRDS(NEON_plot.level_genera_obs_16S.path)	

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
  obs.lo95 <- truth$lo95[,bac_name][order(match(names(truth$lo95[,bac_name]),names(mu)))]	
  obs.hi95 <- truth$hi95[,bac_name][order(match(names(truth$hi95[,bac_name]),names(mu)))]	
  #plot	
  plot(obs.mu ~ mu, cex = 0.7, ylim=c(0,limy), main = paste0('plot-level ', bac_name))	
  arrows(c(mu), obs.lo95, c(mu), obs.hi95, length=0.05, angle=90, code=3)	
  rsq <- round(summary(lm(obs.mu ~mu))$r.squared,2)	
  mtext(paste0('R2=',rsq), side = 3, cex = .7, line = -2.7, adj = 0.03)	
  #1:1 line	
  abline(0,1, lwd = 2)	
  # best-fit line.
  abline(lm(obs.mu ~mu), lty = 2, col = 'purple')
  #add confidence interval.	
  range <- mu	
  polygon(c(range, rev(range)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)	
  polygon(c(range, rev(range)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)	
  #fraction within 95% predictive interval.	
  in_it <- round(sum(obs.mu < pi_0.975 & obs.mu > pi_0.025) / length(obs.mu),2) *100	
  state <- paste0(in_it,'% of observations within 95% prediction interval.')	
  mtext(state,side = 3, cex = .7, line = -1.3, adj = 0.05)






#site.level----

#organize data.
fcast <- output$site.fit
truth <- readRDS(NEON_site.level_genera_obs_16S.path)

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

#for(i in 2:21) {
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
  plot(obs.mu ~ mu, cex = 0.7, ylim=c(0,limy), main=paste0('site-level ', bac_name))
  arrows(c(mu), obs.lo95, c(mu), obs.hi95, length=0.05, angle=90, code=3)
  rsq <- round(summary(lm(obs.mu ~mu))$r.squared,2)
  mtext(paste0('R2=',rsq), cex = .7, side = 3, line = -2.7, adj = 0.03)
  #1:1 line
  abline(0,1, lwd = 2)
  # best-fit line.
  abline(lm(obs.mu ~mu), lty = 2, col = 'purple')
  #add confidence interval.
  range <- mu
  polygon(c(range, rev(range)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
  polygon(c(range, rev(range)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
  #fraction within 95% predictive interval.
  in_it <- round(sum(obs.mu < pi_0.975 & obs.mu > pi_0.025) / length(obs.mu),2) *100
  state <- paste0(in_it,'% of observations within 95% prediction interval.')
  mtext(state,side = 3,  cex = .7, line = -1.3, adj = 0.05)
#}
  mtext(paste0("Rsq for prior fit:", prior_rsq[i]), cex = .7, side = 3, line = 1, outer = TRUE)
  mtext(paste0("Present in ", present_percent[i], " of NEON cores."), cex = .7, side = 3, line = 0, outer = TRUE)
  
}
dev.off()
