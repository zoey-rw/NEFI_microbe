# Validation plots of functional group forecasts to NEON core/plot/sites, 
# plotted by group, each with c/p/s.

rm(list=ls())
source('paths.r')
source('NEFI_functions/crib_fun.r')
library(data.table)


# load forecast 
all_fcasts <- readRDS(NEON_cps_fcast_N.cycler_16S.path)
all_fcasts <- all_fcasts[1:5]

#testing at phylum level.
all_fcasts <- readRDS(paste0(pecan_gen_16S_dir,"/NEON_forecast_data/NEON_cps_fcast_all_phylo_16S.rds"))
#all_fits <- readRDS(paste0(scc_gen_16S_dir,"/JAGS_output/prior_phylo_JAGSfit_phylum.rds"))

all_fits <- readRDS(paste0(scc_gen_16S_dir, "JAGS_output/prior_phylo_JAGSfit_fewer_taxa.rds"))
phylum.mod <- readRDS(paste0(scc_gen_16S_dir,"JAGS_output/prior_phylo_JAGSfit_phylum_fewer_taxa_more_burnin.rds"))
all_fits$phylum <- phylum.mod$phylum

# read in obs table that links deprecatedVialID and geneticSampleID
#map <- readRDS(obs.table_16S.path)
map <- readRDS(core_obs_16S.path)

# read in prior fit.
#all_fits <- readRDS("/fs/data3/caverill/NEFI_data/16S/scc_gen/JAGS_output/bahram_16S.prior_phylo_new_test.rds")

#validate against observed data by plotting.----

#i = 2 

#pdf(NEON_cps.fcast_phyla_16S.path)


#pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/figures/NEON_cps.fcast_phylum_16S_new_pH.pdf")
#pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/figures/NEON_cps.fcast_phylum_16S_moreC.pdf")
filename <- "/fs/data3/caverill/NEFI_data/16S/pecan_gen/figures/NEON_cps.fcast_16S_lessmissingdata.pdf"
pdf(filename,onefile=T)

#global plot settings.----
par(mfrow = c(3,1))
par(mar = c(2,2,2,2))
par(oma = c(0,0,2,0))
trans <- 0.3
limy <- c(0,1)
core.cex <- 0.7
plot.cex <- 1.0
site.cex <- 1.5
outer.cex <- 2
glob.pch <- 16

for (p in 1:length(all_fcasts)) {
   #p <- 1
  output <- all_fcasts[[p]]  
  fit <- all_fits[[p]]
  fit <- fit$no.nutr.preds
  
  # read in observed values
  raw.truth <- readRDS(NEON_all.phylo.levels_plot.site_obs_16S.path)
  raw.truth <- raw.truth[[p]]
  
  for (i in 2:ncol(fit$observed)){
    
    
    #core.level.----
    # read in core-level forecast
    fcast <- output$core.fit
    # read in core-level observed values
    truth <- raw.truth$core.fit
    #organize data.
    truth <- as.data.frame(truth)
    truth$deprecatedVialID <- rownames(truth)
    truth1 <- merge(truth, map[,c("deprecatedVialID", "geneticSampleID")], by = "deprecatedVialID")
    truth1 <- truth1[!duplicated(truth1$geneticSampleID),]
    rownames(truth1) <- gsub('-GEN','',truth1$geneticSampleID)
    truth <- truth1
    truth <- truth[rownames(truth) %in% rownames(fcast$mean),]
    
    for(k in 1:length(fcast)){
      fcast[[k]] <- fcast[[k]][rownames(fcast[[k]]) %in% rownames(truth),]
    }
    #truth <- as.matrix(truth)
    mu <- fcast$mean[,i][order(fcast$mean[,i])]
    ci_0.975 <- fcast$ci_0.975[,i][order(match(names(fcast$ci_0.975[,i]),names(mu)))]
    ci_0.025 <- fcast$ci_0.025[,i][order(match(names(fcast$ci_0.025[,i]),names(mu)))]
    pi_0.975 <- fcast$pi_0.975[,i][order(match(names(fcast$pi_0.975[,i]),names(mu)))]
    pi_0.025 <- fcast$pi_0.025[,i][order(match(names(fcast$pi_0.025[,i]),names(mu)))]
    group_name <- colnames(fcast$mean)[i]
    if (!group_name %in% colnames(truth)) next()
    obs.mu   <- truth[,c(group_name),drop=FALSE][order(match(names(truth[,c(group_name),drop=FALSE]),names(mu)))]
    obs.mu <- obs.mu[,c(group_name)]
    # get prior rsq.
    mod <- betareg::betareg(crib_fun(fit$observed[,i]/rowSums(fit$observed)) ~ crib_fun(fit$predicted[,i]))
    prior_rsq <-round(summary(mod)$pseudo.r.squared, 2)
    
    # create plot
    # get ylim (used for plot/site level as well)
    obs_limit <- max(obs.mu, na.rm = T)
    if(max(pi_0.975) > as.numeric(obs_limit)){obs_limit <- max(pi_0.975)}
    limy <- as.numeric(obs_limit)*1.05
    if(limy > 0.95){limy <- 1}
    
    plot(obs.mu ~ mu, cex = 0.7, ylim=c(0,limy), main = paste0('core-level ', group_name))
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
    
    
    
    # plot.level----
    # read in plot-level forecast 
    fcast <- output$plot.fit	
    # read in plot-level observed values 
    truth <- raw.truth$plot.fit
    
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
    group_name <- colnames(fcast$mean)[i]	
    obs.mu   <- truth$mean[,group_name][order(match(names(truth$mean[,group_name]),names(mu)))]	
    obs.lo95 <- truth$lo95[,group_name][order(match(names(truth$lo95[,group_name]),names(mu)))]	
    obs.hi95 <- truth$hi95[,group_name][order(match(names(truth$hi95[,group_name]),names(mu)))]	
    
    # create plot
    
    plot(obs.mu ~ mu, cex = 0.7, ylim=c(0,limy), main = paste0('plot-level ', group_name))	
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
    
    #read in site-level forecast
    fcast <- output$site.fit
    # read in site-level observed values
    truth <- raw.truth$site.fit
    #organize data.
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
    group_name <- colnames(fcast$mean)[i]
    obs.mu   <- truth$mean[,group_name][order(match(names(truth$mean[,group_name]),names(mu)))]
    obs.lo95 <- truth$lo95[,group_name][order(match(names(truth$lo95[,group_name]),names(mu)))]
    obs.hi95 <- truth$hi95[,group_name][order(match(names(truth$hi95[,group_name]),names(mu)))]
    
    #create plot
    
    plot(obs.mu ~ mu, cex = 0.7, ylim=c(0,limy), main=paste0('site-level ', group_name))
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
    mtext(paste0("Rsq for prior fit:", prior_rsq), cex = .7, side = 3, line = 1, outer = TRUE)
    #mtext(paste0("Present in ", present_percent[i], " of NEON cores."), cex = .7, side = 3, line = 0, outer = TRUE)
  }
}
dev.off()
