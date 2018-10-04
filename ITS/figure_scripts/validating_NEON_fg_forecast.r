#validate the NEON site level functional group forecast!
rm(list=ls())
source('paths.r')
source('NEFI_functions/')

#load and prep data.----
#load forecast.
fcast <- readRDS(NEON_site_fcast_fg.path)
fcast <- fcast$all.preds
#load site-level fg observations.
fg    <- readRDS(NEON_site.level_fg_obs.path)


#order observed and predicted. drop observations not in predicted.
for(i in 1:length(fcast)){
  fcast[[i]] <- fcast[[i]][order(rownames(fcast[[i]])),]
  #fcast[[i]] <- fcast[[i]][!(rownames(fcast[[i]]) == 'DSNY'),] #DSNY is the outlier.
}
for(i in 1:length(fg)){
  fg[[i]] <- fg[[i]][rownames(fg[[i]]) %in% rownames(fcast$mean),]
  fg[[i]] <- fg[[i]][order(rownames(fg[[i]])),]
  fg[[i]] <- fg[[i]][,order(match(colnames(fg[[i]]), colnames(fcast$mean)))]
}


#plot!----
#currently only plotting Ecto fungi
par(mfrow = c(2,2))
limy <- c(0,1)
trans <- 0.3


#isolate predicted mean and interval of interest as well as observed.----
for(i in 1:ncol(fcast$mean)){
  lab <- rownames(fcast$mean)
  mu <- fcast$mean[,i][order(fcast$mean[,i])]
  ci_0.975 <- fcast$ci_0.975[,i][order(match(names(fcast$ci_0.975[,i]),names(mu)))]
  ci_0.025 <- fcast$ci_0.025[,i][order(match(names(fcast$ci_0.025[,i]),names(mu)))]
  pi_0.975 <- fcast$pi_0.975[,i][order(match(names(fcast$pi_0.975[,i]),names(mu)))]
  pi_0.025 <- fcast$pi_0.025[,i][order(match(names(fcast$pi_0.025[,i]),names(mu)))]
  obs.mu   <- fg$mean[,i][order(match(names(fg$mean[,i]),names(mu)))]
  obs.lo95 <- fg$lo95[,i][order(match(names(fg$lo95[,i]),names(mu)))]
  obs.hi95 <- fg$hi95[,i][order(match(names(fg$hi95[,i]),names(mu)))]
  fungi_name <- colnames(fcast$mean)[i]
  
  #calculate how many observations within 95% rpedictive interval.----
  win <- round(length(obs.mu[obs.mu < pi_0.975 & obs.mu > pi_0.025]) / length(obs.mu),2) * 100
  
  #plot predictions.----
  lab <-names(mu)
  plot(mu, ylim = limy, xlab = NA, xaxt='n', bty = 'n', ylab = NA, cex=0)
  axis(1, at=1:length(mu), labels=lab, las = 2)
  lines(smooth.spline(mu, spar = 0), lwd = 2)
  range <- c(1:length(mu))
  polygon(c(range, rev(range)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
  polygon(c(range, rev(range)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
  #plot observed values.
  par(new = T)
  plot(obs.mu, ylim = limy, xlab = NA, xaxt='n', yaxt='n', bty = 'n', ylab =NA, pch = 16)
  arrows(c(1:length(obs.mu)), obs.lo95, c(1:length(obs.mu)), obs.hi95, length=0.05, angle=90, code=3)
  rsq <- round(summary(lm(obs.mu ~ mu))$r.squared, 2)
  mtext(paste0(fungi_name,' fungi forecast validation: NEON'), side = 3, line = -1)
  mtext(paste0('R2 = ',rsq), side = 3, adj = 0.025, line = -4.2)
  mtext(paste0(win,'% of sites within 95% predictive interval.'), side = 3, line = -3, adj = 0.05)
}


#try again with an x axis that is values. This definitely looks shittier.----
#lab <-names(mu)
#plot(obs.mu ~ mu, ylim = limy, xlab = NA, bty = 'n', ylab = NA, cex=0)
#axis(1, at=mu, labels=lab, las = 2)
#abline(0,1, lwd = 2)
#polygon(c(mu, rev(mu)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
#polygon(c(mu, rev(mu)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
#plot observed values.
#par(new = T)
#plot(obs.mu ~ mu, ylim = limy, xlab = NA, xaxt='n', yaxt='n', bty = 'n', ylab =NA, pch = 16)
#arrows(c(mu), obs.lo95, c(mu), obs.hi95, length=0.05, angle=90, code=3)
#rsq <- round(summary(lm(obs.mu ~ mu))$r.squared, 2)
#mtext('Ectomycorrhizal forecast validation: NEON', side = 3, line = -1)
#mtext(paste0(win,'% of sites within 95% predictive interval.'), side = 3, line = -3, adj = 0.05)
#mtext(paste0('R2 = ',rsq), side = 3, adj = 0.025, line = -4.2)


