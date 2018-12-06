# Get plot- and site-level hierarchical means of cosmopolitan genera from rel abundances at NEON cores
# top 15 phyla
# clear workspace, load libraries
rm(list=ls())
library(runjags)
source('paths.r')
#source('NEFI_functions/hierarch_ddirch_means.r')

# source hierarch means function
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/hierarch_ddirch_means.r", ssl.verifypeer = FALSE)
eval(parse(text = script))
# colin's paths.R
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/paths.r", ssl.verifypeer = FALSE)
eval(parse(text = script))


#set output.path----
plot.output.path <- NEON_plot.level_phyla_obs_16S.path
site.output.path <- NEON_site.level_phyla_obs_16S.path

#load data and format.----
d <- readRDS(NEON_phyla_abundances_16S.path)
fg.c <- d$abundances
fg.c$other <- NULL
# drop anything with under 1000 sequences.
seq_total <- d$seq_total[d$seq_total>1000]
fg.c <- fg.c[d$seq_total>1000,]

#get y dependent matrix.----
y <- fg.c
rownames(y) <- NULL
y <- as.matrix(y)
y <- y+1 #dirichlet doesn't like hard zeros because log-link.
y <- y/seq_total
other <- 1- rowSums(y)
y <- cbind(other,y)


#get core_plot, core_site, plot_site
deprecatedVialID <- rownames(fg.c)
#deprecatedVialID <- str_replace_all(deprecatedVialID, c("."= "_"))
core_plot <- substr(deprecatedVialID,1,8)
core_site <- substr(deprecatedVialID,1,4)
plot_site <- unique(core_plot)
plot_site <- substr(plot_site,1,4)


#--------Get means using hierarch_ddirch_means function.-------#
# takes ~1hr on geo.

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


#--------- visual check --------#

site.fit <- readRDS(site.output.path)
site.fit <- readRDS("/fs/data3/caverill/NEFI_data/16S/scc_gen/NEON_processed/NEON_site.level_phyla_obs_16S.rds")
plot.fit <- readRDS(plot.output.path)
plot.fit <- readRDS("/fs/data3/caverill/NEFI_data/16S/scc_gen/NEON_processed/NEON_plot.level_phyla_obs_16S.rds")


plot_check <- T
if(plot_check == T){
  #DO THESE MODELS BALL-PARK FIT STRAIGHTFORWARD MEANS?
  #Fit doesn't need to be perfect and can be biased, just need to be ~linear.
  #check if this is similar to just aggregating by site.
  par(mfrow = c(1,2))
  test.mu <- data.frame(core_plot,core_site, y)
  #setnames(test.mu, old = "Candidatus.Solibacter", new = "Candidatus Solibacter")
  
  #site level.
  k <- aggregate(. ~ core_site, FUN=mean, data = test.mu[,colnames(test.mu) %in% colnames(y)])
  to_plot <- as.vector(as.matrix(k[,2:ncol(k)]))
  plot(to_plot ~ as.vector(site.fit$mean));abline(0,1, lwd =2)
  mtext('Site Level', side = 3)
  #plot level.
  k <- aggregate(. ~ core_plot, FUN=mean, data = test.mu[,colnames(test.mu) %in% colnames(y)])
  to_plot <- as.vector(as.matrix(k[,2:ncol(k)]))
  plot(to_plot ~ as.vector(plot.fit$mean));abline(0,1, lwd =2)
  mtext('Plot Level', side = 3)
}




