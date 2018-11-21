#getting plot and site level means of NEON cosmopolitan genera for forecast validation.
#This would be worth throwing on the cluster. Runs in like 20-30 on pecan2.
rm(list=ls())
library(runjags)
source('paths.r')
source('NEFI_functions/hierarch_ddirch_means.r')

#set output.path----
plot.output.path <- NEON_plot.level_genera_obs_fastq.path
site.output.path <- NEON_site.level_genera_obs_fastq.path

#load data and format.----
d <- readRDS(NEON_ITS_fastq_cosmo_genera.path)
fg.c <- d$abundances
seq.depth <- d$seq_total

#get y dependent matrix.----
y <- fg.c
y$geneticSampleID <- NULL
y <- as.matrix(y)
y <- y+1 #dirichlet doesn't like hard zeros because log-link.
y <- y/d$seq_total

#get core_plot, core_site, plot_site
core_plot <- substr(fg.c$geneticSampleID,1,8)
core_site <- substr(fg.c$geneticSampleID,1,4)
plot_site <- unique(core_plot)
plot_site <- substr(plot_site,1,4)

#Get means using hierarch_ddirch_means function.----
fit <- hierarch_ddirch_means(y = y, core_plot = core_plot, plot_site = plot_site)
plot.fit <- fit$plot.fit
site.fit <- fit$site.fit

#add row and column names - plot level.
for(i in 1:length(plot.fit)){
  rownames(plot.fit[[i]]) <- unique(core_plot)
  colnames(plot.fit[[i]]) <- colnames(y)
}

#add row and column names - site level.
for(i in 1:length(site.fit)){
  rownames(site.fit[[i]]) <- unique(plot_site)
  colnames(site.fit[[i]]) <- colnames(y)
}

#save matrix lists.----
saveRDS(site.fit, site.output.path)
saveRDS(plot.fit, plot.output.path)

#visual check
plot_check <- F
if(plot_check == T){
  #DO THESE MODELS BALL-PARK FIT STRAIGHTFORWARD MEANS?
  #Fit doesn't need to be perfect and can be biased, just need to be ~linear.
  #check if this is similar to just aggregating by site.
  par(mfrow = c(1,2))
  #site level.
  test.mu <- data.frame(core_plot,core_site, y)
  k <- aggregate(. ~ core_site, FUN=mean, data = test.mu[,colnames(y)])
  to_plot <- as.vector(as.matrix(k[,2:ncol(k)]))
  plot(to_plot ~ as.vector(site.fit$mean));abline(0,1, lwd =2)
  mtext('Site Level', side = 3)
  #plot level.
  test.mu <- data.frame(core_plot,core_site, y)
  k <- aggregate(. ~ core_plot, FUN=mean, data = test.mu[,colnames(y)])
  to_plot <- as.vector(as.matrix(k[,2:ncol(k)]))
  plot(to_plot ~ as.vector(plot.fit$mean));abline(0,1, lwd =2)
  mtext('Plot Level', side = 3)
}
