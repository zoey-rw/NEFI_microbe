#plotting functional group forecasts vs. validation data.
#clear environment, source paths, functions and packages.
rm(list=ls())
source('paths.r')

#set output path.----
output.path <- NEON_cps_genera_forecast_figure.path
#output.path <- 'test.png'

#load data.----
d <- readRDS(NEON_site_fcast_genera.path)
core.truth <- readRDS(NEON_cosmo_genera.path)
core.truth <- core.truth$rel.abundances
rownames(core.truth) <- gsub('-GEN','',core.truth$geneticSampleID)
plot.truth <- readRDS(NEON_plot.level_genera_obs.path)
site.truth <- readRDS(NEON_site.level_genera_obs.path)


#global plot settings.----
trans <- 0.3
limy <- c(0,1)
core.cex <- 0.7
plot.cex <- 1.0
site.cex <- 1.5
outer.cex <- 2
glob.pch <- 16
names <- colnames(d$core.fit$mean)
names <- names[c(2:(length(names)),1)]
names <- names[!names %in% c('other')]

#setup output spec.----
png(filename=output.path,width=12,height=12,units='in',res=300)
par(mfrow = c(7,3),
    mai = c(0.3,0.3,0.3,0.3),
    oma = c(4,6,3,1))

#loop over functional groups.----
for(i in 1:length(names)){
  #core.level.----
  #organize data.
  fcast <- d$core.fit
  obs <- core.truth
  obs <- obs[rownames(obs) %in% rownames(fcast$mean),]
  for(k in 1:length(fcast)){
    fcast[[k]] <- fcast[[k]][rownames(fcast[[k]]) %in% rownames(obs),]
  }
  obs <- as.matrix(obs)
  pos <- which(colnames(fcast$mean) == names[i]) #position that matches the fungal type we are plotting.
  mu <- fcast$mean[,pos][order(fcast$mean[,pos])]
  ci_0.975 <- fcast$ci_0.975[,pos][order(match(names(fcast$ci_0.975[,pos]),names(mu)))]
  ci_0.025 <- fcast$ci_0.025[,pos][order(match(names(fcast$ci_0.025[,pos]),names(mu)))]
  pi_0.975 <- fcast$pi_0.975[,pos][order(match(names(fcast$pi_0.975[,pos]),names(mu)))]
  pi_0.025 <- fcast$pi_0.025[,pos][order(match(names(fcast$pi_0.025[,pos]),names(mu)))]
  fungi_name <- colnames(fcast$mean)[pos]
  obs.pos <- which(colnames(obs) == names[i])
  obs.mu   <- obs[,obs.pos][order(match(names(obs[,c(fungi_name)]),names(mu)))]
  
  #get y-limit.
  obs_limit <- max(obs.mu)
  if(max(pi_0.975) > as.numeric(obs_limit)){obs_limit <- max(pi_0.975)}
  y_max <- as.numeric(obs_limit)*1.05
  if(y_max > 0.95){y_max <- 1}
  
  #plot
  plot(obs.mu ~ mu, cex = core.cex, pch=glob.pch, ylim=c(0,y_max), ylab=NA, xlab = NA)
  rsq <- round(summary(lm(obs.mu ~mu))$r.squared,2)
  mtext(paste0('R2=',rsq), side = 3, line = -2.7, adj = 0.03)
  mtext(fungi_name, side = 2, line = 2.5, cex = 1.5)
  #add confidence interval.
  range <- mu
  polygon(c(range, rev(range)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
  polygon(c(range, rev(range)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
  #fraction within 95% predictive interval.
  in_it <- round(sum(as.numeric(obs.mu) < pi_0.975 & as.numeric(obs.mu) > pi_0.025) / length(obs.mu),2) * 100
  state <- paste0(in_it,'% of obs. within interval.')
  mtext(state,side = 3, line = -1.3, adj = 0.05)
  abline(0,1,lwd=2)
  
  
  #plot.level.----
  #organize data.
  fcast <- d$plot.fit
  obs <- plot.truth
  for(k in 1:length(obs)){
    rownames(obs[[k]]) <- gsub('.','_',rownames(obs[[k]]), fixed = T)
    obs[[k]] <- obs[[k]][rownames(obs[[k]]) %in% rownames(fcast$mean),]
  }
  for(k in 1:length(fcast)){
    fcast[[k]] <- fcast[[k]][rownames(fcast[[k]]) %in% rownames(obs$mean),]
  }
  pos <- which(colnames(fcast$mean) == names[i]) #position that matches the fungal type we are plotting.
  mu <- fcast$mean[,pos][order(fcast$mean[,pos])]
  ci_0.975 <- fcast$ci_0.975[,pos][order(match(names(fcast$ci_0.975[,pos]),names(mu)))]
  ci_0.025 <- fcast$ci_0.025[,pos][order(match(names(fcast$ci_0.025[,pos]),names(mu)))]
  pi_0.975 <- fcast$pi_0.975[,pos][order(match(names(fcast$pi_0.975[,pos]),names(mu)))]
  pi_0.025 <- fcast$pi_0.025[,pos][order(match(names(fcast$pi_0.025[,pos]),names(mu)))]
  fungi_name <- colnames(fcast$mean)[pos]
  obs.mu   <- obs$mean[,fungi_name][order(match(names(obs$mean[,fungi_name]),names(mu)))]
  obs.lo95 <- obs$lo95[,fungi_name][order(match(names(obs$lo95[,fungi_name]),names(mu)))]
  obs.hi95 <- obs$hi95[,fungi_name][order(match(names(obs$hi95[,fungi_name]),names(mu)))]
  
  #get y-limit.
  obs_limit <- max(obs.hi95)
  if(max(pi_0.975) > obs_limit){obs_limit <- max(pi_0.975)}
  y_max <- as.numeric(obs_limit)*1.05
  if(y_max > 0.95){y_max <- 1}
  
  #plot
  plot(obs.mu ~ mu, cex = plot.cex, pch=glob.pch, ylim=c(0,y_max), ylab=NA, xlab = NA)
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
  state <- paste0(in_it,'% of obs. within interval.')
  mtext(state,side = 3, line = -1.3, adj = 0.05)
  
  
  #site.level.----
  #organize data.
  fcast <- d$site.fit
  obs <- site.truth
  for(k in 1:length(obs)){
    rownames(obs[[k]]) <- gsub('.','_',rownames(obs[[k]]), fixed = T)
    obs[[k]] <- obs[[k]][rownames(obs[[k]]) %in% rownames(fcast$mean),]
  }
  for(k in 1:length(fcast)){
    fcast[[k]] <- fcast[[k]][rownames(fcast[[k]]) %in% rownames(obs$mean),]
  }
  pos <- which(colnames(fcast$mean) == names[i]) #position that matches the fungal type we are plotting.
  mu <- fcast$mean[,pos][order(fcast$mean[,pos])]
  ci_0.975 <- fcast$ci_0.975[,pos][order(match(names(fcast$ci_0.975[,pos]),names(mu)))]
  ci_0.025 <- fcast$ci_0.025[,pos][order(match(names(fcast$ci_0.025[,pos]),names(mu)))]
  pi_0.975 <- fcast$pi_0.975[,pos][order(match(names(fcast$pi_0.975[,pos]),names(mu)))]
  pi_0.025 <- fcast$pi_0.025[,pos][order(match(names(fcast$pi_0.025[,pos]),names(mu)))]
  fungi_name <- colnames(fcast$mean)[pos]
  obs.mu   <- obs$mean[,fungi_name][order(match(names(obs$mean[,fungi_name]),names(mu)))]
  obs.lo95 <- obs$lo95[,fungi_name][order(match(names(obs$lo95[,fungi_name]),names(mu)))]
  obs.hi95 <- obs$hi95[,fungi_name][order(match(names(obs$hi95[,fungi_name]),names(mu)))]
  
  #get y-limit.
  obs_limit <- max(obs.hi95)
  if(max(pi_0.975) > obs_limit){obs_limit <- max(pi_0.975)}
  y_max <- as.numeric(obs_limit)*1.05
  if(y_max > 0.95){y_max <- 1}
  
  #plot
  plot(obs.mu ~ mu, cex = site.cex, pch=glob.pch, ylim=c(0,y_max), ylab=NA, xlab = NA)
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
  state <- paste0(in_it,'% of obs. within interval.')
  mtext(state,side = 3, line = -1.3, adj = 0.05)
}

#outer labels.----
mtext('core-level', side = 3, line = -1.8, adj = 0.12, cex = outer.cex, outer = T)
mtext('plot-level', side = 3, line = -1.8, adj = 0.50, cex = outer.cex, outer = T)
mtext('site-level', side = 3, line = -1.8, adj = 0.88, cex = outer.cex, outer = T)
mtext( 'observed relative abundance', side = 2, line = 3, cex = outer.cex, outer = T)
mtext('predicted relative abundance', side = 1, line = 2, cex = outer.cex, outer = T)

#end plot.----
dev.off()
