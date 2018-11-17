#getting site level means of NEON cosmopolitan genera of interest for forecast validation.
#Model takes 5-10 minutes to run and summarize results.
rm(list=ls())
library(runjags)
source('paths.r')

library(RCurl)
# source paths.r from colins github
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/paths.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

#set output.path----
plot.output.path <- NEON_plot.level_genera_obs.path
site.output.path <- NEON_site.level_genera_obs.path

#load data and format.----
d <- readRDS(NEON_cosmo_abundances.path)
fg.c <- d$abundances
fg.c <- fg.c[,!colnames(fg.c) %in% c('deprecatedVialID','other')]

#get y dependent matrix.----
y <- fg.c
rownames(y) <- NULL
y$geneticSampleID <- NULL
y <- as.matrix(y)
y <- y+1 #dirichlet doesn't like hard zeros because log-link.
# let's add 20 to the seq_total, since we're adding 20 sequences.
d$seq_total <- d$seq_total + 21
y <- y/d$seq_total
other <- 1- rowSums(y)
y <- cbind(other,y)
# 
# # many cores have less than 20 counts total, which creates a problematic values when we divide y+1 by the seq_total. 
# # let's remove any negative abundances we end up with.
# this creates a "Slicer stuck at value with infinite density" error...
# y.df <- as.data.frame(y)
# has.neg <- apply(y.df, 1, function(row) any(row < 0)) 
# y.df <- y.df[-which(has.neg),]
# fg.c <- fg.c[-which(has.neg),]
# y <- as.matrix(y.df)

#get core_plot, core_site, plot_site
deprecatedVialID <- rownames(fg.c)
#deprecatedVialID <- str_replace_all(deprecatedVialID, c("."= "_"))
core_plot <- substr(deprecatedVialID,1,8)
core_site <- substr(deprecatedVialID,1,4)
plot_site <- unique(core_plot)
plot_site <- substr(plot_site,1,4)

#jags data object.----
jd <- list(y=as.matrix(y), N.site = length(unique(plot_site)), N.plot =  length(plot_site), N.core = nrow(y), N.spp = ncol(y), 
           core_plot = as.factor(core_plot), plot_site=as.factor(plot_site))

#jags model.----
jags.model = "
model {
# Observations (single set per core):
for(i in 1:N.core){
y[i,1:N.spp] ~ ddirch(plot_mu[core_plot[i],1:N.spp] * core_alpha) 
}
# Plot means:
for(i in 1:N.plot){
plot_mu[i,1:N.spp] ~ ddirch(site_mu[plot_site[i],1:N.spp] * plot_alpha)
}
# Site means:
for(i in 1:N.site){
site_mu[i,1:N.spp] ~ ddirch(priors[1:N.spp] * site_alpha)
}
# Priors:
for(s in 1:N.spp){
priors[s] <- 1
}
site_alpha ~ dgamma(0.01, 0.01)
plot_alpha ~ dgamma(0.01, 0.01)
core_alpha ~ dgamma(0.01, 0.01) 
}" #end jags model.

#run JAGS model in runjags.----
test <- run.jags(model = jags.model,
                 data = jd,
                 n.chains = 3,
                 monitor = c('core_alpha','plot_alpha','site_alpha','plot_mu','site_mu'),
                 adapt = 1000,
                 burnin = 3500,
                 sample = 5000)

#get matrices of mean, upper and lower confidence intervals at plot and site level.----
#site level
site.out <- summary(test, var='site_mu')
mean.matrix <- matrix(site.out[,4], nrow = length(unique(plot_site)), ncol = ncol(y))
lo95.matrix <- matrix(site.out[,1], nrow = length(unique(plot_site)), ncol = ncol(y))
hi95.matrix <- matrix(site.out[,3], nrow = length(unique(plot_site)), ncol = ncol(y))

site.output.list <- list(mean.matrix, lo95.matrix, hi95.matrix)
names(site.output.list) <- c('mean','lo95','hi95')
#add row and column names
for(i in 1:length(site.output.list)){
  rownames(site.output.list[[i]]) <- unique(plot_site)
  colnames(site.output.list[[i]]) <- colnames(y)
}

#plot level
plot.out <- summary(test, var = 'plot_mu')
mean.matrix <- matrix(plot.out[,4], nrow = length(plot_site), ncol = ncol(y))
lo95.matrix <- matrix(plot.out[,1], nrow = length(plot_site), ncol = ncol(y))
hi95.matrix <- matrix(plot.out[,3], nrow = length(plot_site), ncol = ncol(y))

plot.output.list <- list(mean.matrix, lo95.matrix, hi95.matrix)
names(plot.output.list) <- c('mean','lo95','hi95')
#add row and column names
for(i in 1:length(plot.output.list)){
  rownames(plot.output.list[[i]]) <- unique(core_plot)
  colnames(plot.output.list[[i]]) <- colnames(y)
}


#check if this is similar to just aggregating by site.
test.mu <- data.frame(core_plot,core_site, y)
k <- aggregate(. ~ core_site, FUN=mean, data = test.mu[,colnames(y)])
plot(k$Russula ~ site.output.list$mean[,3]);abline(0,1)

#check if this is similar to just aggregating by plot.
k <- aggregate(. ~ core_plot, FUN=mean, data = test.mu[,colnames(y)])
plot(k$Russula ~ plot.output.list$mean[,3]);abline(0,1)
 