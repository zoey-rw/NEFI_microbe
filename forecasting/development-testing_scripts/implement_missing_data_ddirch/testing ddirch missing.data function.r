#testing ddirch with missing data.
#clear environment, generate data.
rm(list = ls())
library(data.table)
source('NEFI_functions/crib_fun.r')

#### generate site-level and hierarchical multivariate data.####
n <- 10*4
mat <- c(rep(2,n/4),rep(4,n/4), rep(6,n/4), rep(8,n/4))
map <- c(rep(400,n/2), rep(200,n/2))
pH <- runif(n, 3,7)
core <- (data.frame(mat,map,pH))
core$y <- mat*0.5 + map*0.04 +pH*1 + rnorm(n,2)
core$plot <- c(rep('a1',n/4),rep('a2',n/4), rep('b1',n/4), rep('b2',n/4))
core$site <- c(rep('a',n/2), rep('b',n/2))
core$intercept <- rep(1, nrow(core))
core_plot <- as.factor(core$plot)
core_site <- as.factor(core$site)
plot_plot <- as.factor(unique(core$plot))
plot_site <- as.factor(substr(plot_plot,1,1))
site_site <- as.factor(unique(core$site))
y <- as.matrix(core$y)
#get two species abundances
y1 <- y
y2 <- rep(10, length(y1))
spp.y <- data.frame(cbind(y1,y2))
spp.y <- spp.y / rowSums(spp.y)

core_mu <- data.frame(core$intercept, core$pH)
site_mu <- data.frame(c(400,200))
plot_mu <- data.frame(c(2,4,6,8))
colnames(core_mu) <- c('intercept','pH')
colnames(site_mu) <- c('map')
colnames(plot_mu) <- c('mat')

#convert to just "site-level" for downstream testing.
core <- cbind(core_mu, core_plot)
plot <- cbind(plot_mu, plot_site)
plot <- cbind(plot, plot_plot)
site <- cbind(site_mu, site_site)
all.test <- merge(core,plot, by.x = 'core_plot',by.y ='plot_plot')
all.test <- merge(all.test,site, by.x = 'plot_site', by.y = 'site_site')
test.x <- all.test[,3:ncol(all.test)]




#### testing hierarchical dirichlet w/ pseudo data ####
source('forecasting/development-testing_scripts/implement_missing_data_ddirch/ddirch_missing.data_function.r')

#go ahead and use the function to fit the model to data, summarize. 
test <- hierarchical_dirlichet_jags(y = spp.y,
                                    core_mu = core_mu, plot_mu = plot_mu, site_mu = site_mu,
                                    core_plot = core_plot, core_site = core_site, plot_site = plot_site, 
                                    plot_plot = plot_plot, site_site = site_site,
                                    adapt = 100, burnin = 200, sample = 200)

#Plot predicted vs. observed, 1:1 line, best fit line and R2 values 
par(mfrow = c(1,2))
m1 <- lm(test$observed[,1] ~ test$predicted[,1])
plot(test$observed[,1] ~ test$predicted[,1], main = 'y1 obs vs. pred'); abline(0,1, lwd = 2);abline(m1, lty = 2, col = 'green')
mtext(paste0('R2 = ',round(summary(m1)$r.squared, 2)), side = 3, adj = 0.05, line = -2)
m2 <- lm(test$observed[,2] ~ test$predicted[,2])
plot(test$observed[,2] ~ test$predicted[,2], main = 'y2 obs vs. pred'); abline(0,1, lwd = 2);abline(m2, lty = 2, col = 'green')
mtext(paste0('R2 = ',round(summary(m1)$r.squared, 2)), side = 3, adj = 0.05, line = -2)





#### testing site level dirichlet w/ pseudo data ####
source('forecasting/development-testing_scripts/implement_missing_data_ddirch/ddirch_missing.data_site.level.only_function.r')
test <- site.level_dirlichet_jags(y=spp.y,x_mu=test.x,adapt = 100, burnin = 200, sample = 200)

#plot predicted vs. observed, 1:1 line, vest fit line and R2 values.
par(mfrow = c(1,2))
plot(test$observed[,1] ~ test$predicted[,1])
mod <- lm(test$observed[,1] ~ test$predicted[,1])
abline(0,1,lwd=2);abline(mod, lty = 2, col = 'green')
mtext(paste0('R2 = ',round(summary(mod)$r.squared,2)),side = 3, line = -2, adj= 0.05)
plot(test$observed[,2] ~ test$predicted[,2])
mod <- lm(test$observed[,2] ~ test$predicted[,2])
abline(0,1,lwd=2);abline(mod, lty = 2, col = 'green')
mtext(paste0('R2 = ',round(summary(mod)$r.squared,2)),side = 3, line = -2, adj= 0.05)




#### testing site-level dirichlet w/ tedersoo prior data ####
#load data
d <- readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/ted_all_prior_data.rds')
d <- d[,.(Ectomycorrhizal,Saprotroph,Arbuscular,cn,pH,moisture,NPP,map,mat)]
d <- d[complete.cases(d),] #optional. This works with missing data.

#organize y data
y <- d[,.(Ectomycorrhizal,Saprotroph,Arbuscular)]
#make other column
y <- data.frame(lapply(y,crib_fun))
y$other <- 1 - rowSums(y)
y <- as.data.frame(y)
#reorder columns. other needs to be first.
y <- y[c('other','Ectomycorrhizal','Saprotroph','Arbuscular')]

#Drop in intercept, setup predictor matrix.
d$intercept <- rep(1,nrow(d))
x <- d[,.(intercept,cn,pH,moisture,NPP,mat,map)]
x$map <- log(x$map)

#source and run function
source('forecasting/development-testing_scripts/implement_missing_data_ddirch/ddirch_missing.data_site.level.only_function.r')
test <- site.level_dirlichet_jags(y=y,x_mu=x,adapt = 50, burnin = 50, sample = 100)

#plot predicted vs. observed, 1:1 line, vest fit line and R2 values.
#4 "species" modeled, setting up 2x2 grid.
par(mfrow = c(2,2))
y.names <- names(test$species_parameter_output)
for(i in 1:ncol(test$predicted)){
  plot(test$observed[,i] ~ test$predicted[,i], main = paste(y.names[i]))
  mod <- lm(test$observed[,i] ~ test$predicted[,i])
  beta.mod <- betareg::betareg(test$observed[,i] ~ test$predicted[,i])
  abline(0,1,lwd=2);abline(mod, lty = 2, col = 'green')
  mtext(paste0('R2 = ',round(summary(beta.mod)$pseudo.r.squared,3)),side = 3, line = -2, adj= 0.05)
}

#Does this beat simple beta-regression? Looks about the same. Thank god we did all this work to get the same answer...
#Note- both fit the Arbuscular data very poorly, regardless of what the R2 values say.
for(i in 1:ncol(y)){
  mod <- betareg::betareg(y[,i] ~ x$cn + x$pH + x$moisture + x$NPP + x$mat + x$map)
  cat(paste0(y.names[i],': R2 = ',round(summary(mod)$pseudo.r.squared,3)),'\n')
}



#### testing hierarchical dirichlet w/ subset of NEON data. ####
##Current hangup: some sites have observations but no core or plot level data, which throws off indexing.
rm(list=ls())
source('/home/caverill/NEFI_microbe/data_formatting/formatting_NEON_microbial/core_site_plot_aggregation_May.2018/0. aggregation paths.r')
#source functions
source('forecasting/development-testing_scripts/implement_missing_data_ddirch/ddirch_missing.data_function.r')
source('NEFI_functions/crib_fun.r')

#load data.
 obs.table <- readRDS( obs.table.path)
core.table <- readRDS(core.table.path)
plot.table <- readRDS(plot.table.path)
site.table <- readRDS(site.table.path)

##TESTING
#lets only look at complete cases for 2 sites right now.
 obs.table <-  obs.table[obs.table$siteID %in% c('BART','HARV','TALL'),]
core.table <- core.table[core.table$geneticSampleID %in% obs.table$geneticSampleID,]
plot.table <- plot.table[plot.table$plotID %in% obs.table$plotID,]
site.table <- site.table[site.table$siteID %in% obs.table$siteID,]

#put everything in correct order. 
#Convert all factor vectors to character vectors - this line sucks and is a good argument for data.table/tidyverse
 obs.table[sapply( obs.table, is.factor)] <-  lapply( obs.table[sapply( obs.table, is.factor)], as.character)
core.table[sapply(core.table, is.factor)] <-  lapply(core.table[sapply(core.table, is.factor)], as.character)
plot.table[sapply(plot.table, is.factor)] <-  lapply(plot.table[sapply(plot.table, is.factor)], as.character)
site.table[sapply(site.table, is.factor)] <-  lapply(site.table[sapply(site.table, is.factor)], as.character)

#hierarchically order the rows as site/plot/core.
 obs.table <-  obs.table[order( obs.table$siteID, obs.table$plotID, obs.table$geneticSampleID),]
core.table <- core.table[order(core.table$siteID,core.table$plotID,core.table$geneticSampleID),]
plot.table <- plot.table[order(plot.table$siteID,plot.table$plotID),]
site.table <- site.table[order(site.table$siteID),]
core_mu <- core.table
plot_mu <- plot.table
site_mu <- site.table

#setup indexing
core_plot <- as.factor(core.table$plotID)
core_site <- as.factor(core.table$siteID)
plot_plot <- as.factor(plot.table$plotID)
plot_site <- as.factor(substring(plot_plot,1,4))
site_site <- as.factor(site.table$siteID)

#setup y, core, plot and site data frames, and associated sd tables where applicable.
y <- obs.table[,2:5]
#dropping soil moisture and soil pH for now because all NA values for two site-dates of interest.
#core <- core.table[,c('siteID','plotID','soilMoisture','soilInWaterpH','soilTemp')]
#colnames(core)[3:5] <- c('moisture','pH','temperature')
core <- core.table[,c('siteID','plotID','soilTemp')]
colnames(core)[3] <- c('temperature')
core_sd <- core.table[,c('siteID','plotID')]
plot <- plot.table[,c('siteID','plotID','organicCPercent','CNratio')]
plot_sd <- plot.table[,c('siteID','plotID','organicCPercent_sd','CNratio_sd')]
site <- site.table[,c('siteID','elevation','map','mat','n.dep')]
site_sd <- site.table[,c('siteID','map_sd','mat_sd')]

#drop plot-site indices from data frames. Those saved in index vectors.
core$intercept <- rep(1,nrow(core))
core <- core[,-c(1:2)]; core <- core[,c('intercept','temperature')]
plot <- plot[,-c(1:2)]
site <- site[,-c(1:2)]

core_sd <- NA
plot_sd <- NA
site_sd <- NA

#reorder the y columns
#drop AM fungi- no reads.
y$other <- NULL
y$Arbuscular <- NULL
y <- data.frame(lapply(y,crib_fun))
y$other <- 1 - rowSums(y)
y <- y[c('other','Ectomycorrhizal','Saprotroph')]

#fit the model.
test <- hierarchical_dirlichet_jags(y = y,
                                    core_mu = core, plot_mu = plot, site_mu = site,
                                    core_plot = core_plot, core_site = core_site, plot_site = plot_site, 
                                    plot_plot = plot_plot, site_site = site_site,
                                    adapt = 100, burnin = 200, sample = 200)

#plot results
par(mfrow = c(2,2))
y.names <- names(test$species_parameter_output)
for(i in 1:ncol(test$predicted)){
  plot(test$observed[,i] ~ test$predicted[,i], main = paste(y.names[i]))
  mod <- lm(test$observed[,i] ~ test$predicted[,i])
  beta.mod <- betareg::betareg(test$observed[,i] ~ test$predicted[,i])
  abline(0,1,lwd=2);abline(mod, lty = 2, col = 'green')
  mtext(paste0('R2 = ',round(summary(beta.mod)$pseudo.r.squared,3)),side = 3, line = -2, adj= 0.05)
}

