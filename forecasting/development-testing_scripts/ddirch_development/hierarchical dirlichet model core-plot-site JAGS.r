#simulate some data with just core, plot and site intercepts.
#2 sites, 2 plots per site, 10 cores per plot = 40 obs
#This works with actual function. Currently in script hashed out.
#I also have the function broken out below which was important for trouble shooting.
rm(list = ls())

#source important functions.
source('/home/caverill/NEFI_microbe/NEFI_functions/crib_fun.r')
source('/home/caverill/NEFI_microbe/NEFI_functions/z_transform.r')
source('/home/caverill/NEFI_microbe/NEFI_functions/sd_to_precision.r')
source('/home/caverill/NEFI_microbe/NEFI_functions/precision_matrix_match.r')
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

#WORKS WITH FUNCTION.
#source('NEFI_functions/hierarchical_linear_dirlichet_jags.r')
#test <- hierarchical_dirlichet_jags.r(y = spp.y, core_mu = core_mu, plot_mu = plot_mu, site_mu = site_mu, core_plot = core_plot, core_site = core_site, plot_plot = plot_plot, site_site = site_site, adapt = 100, burnin = 200, sample = 200)
#par(mfrow=c(1,2)) #plot predicted vs. observed.
#line is shifted relative to true values. This implies the intercept is wrong.
#plot(test$observed[,1] ~ test$predicted[,1]);abline(0,1, lwd = 2)
#plot(test$observed[,2] ~ test$predicted[,2]);abline(0,1, lwd = 2)


y = spp.y; core_mu = core_mu; plot_mu = plot_mu; site_mu = site_mu; core_plot = core_plot; core_site = core_site; plot_plot = plot_plot; site_site = site_site
adapt = 100; burnin = 200; sample = 200; n.chains = 3; parallel = F
core_sd = NA; plot_sd = NA; site_sd = NA

#####BELOW IS INTERIOR OF ACTUAL DIRICHLET FUNCTION.

#grab names
y.names <- colnames(y)
core.names <- colnames(core_mu)
plot.names <- colnames(plot_mu)
site.names <- colnames(site_mu)

###massage your data together.
#deal with zero relative abundances.
y <- data.frame(y)
y <- data.frame(lapply(y, crib_fun))

#make sd objects if they were not supplied.
if(is.na(core_sd)){core_sd = data.frame(rep(1,nrow(core_mu)))}
if(is.na(plot_sd)){plot_sd = data.frame(rep(1,nrow(plot_mu)))}
if(is.na(site_sd)){site_sd = data.frame(rep(1,nrow(site_mu)))}

#Match up predictors and their SD. if no SD supplied we assign ~perfect precision.
core_sd <- precision_matrix_match(core_mu,core_sd)
plot_sd <- precision_matrix_match(plot_mu,plot_sd)
site_sd <- precision_matrix_match(site_mu,site_sd)

#covert sd to precision. output is matrix.
core_precision <- sd_to_precision(core_sd)
plot_precision <- sd_to_precision(plot_sd)
site_precision <- sd_to_precision(site_sd)

#make sure every else is a matrix.
y <- as.matrix(y)
core_mu <- as.matrix(core_mu)
plot_mu <- as.matrix(plot_mu)
site_mu <- as.matrix(site_mu)

#make sure index vectors are factors, drop levels.
core_plot <- droplevels(as.factor(core_plot))
core_site <- droplevels(as.factor(core_site))
plot_plot <- droplevels(as.factor(plot_plot))
site_site <- droplevels(as.factor(site_site))

###setup jags data object.
jags.data <- list(N.core = nrow(core_mu), N.plot = nrow(plot_mu), N.site = nrow(site_mu), N.spp = ncol(y),            #number of c/p/s observations and number of species
                  N.core.preds = ncol(core_mu), N.plot.preds = ncol(plot_mu), N.site.preds = ncol(site_mu),           #number of predictors at c/p/s level
                  core_mu = core_mu, plot_mu = plot_mu, site_mu = site_mu,                                            #x-value mean      matrices c/p/s
                  core_precision = core_precision, plot_precision = plot_precision, site_precision = site_precision,  #x-value precision matrices c/p/s
                  core_plot = core_plot, core_site = core_site, plot_plot = plot_plot, site_site = site_site,         #vectors that index across  c/p/s
                  y = y)                                                                                              #species matrix, y

###specify JAGS model.
###specify JAGS model.
jags.model = "
model {
#priors for each species*predictor*hierarchical-level
alpha ~ dnorm(0, 1.0E-3) 
#core level priors
for(i in 1:N.core.preds){
core.mm[i,1] <- 0
for (j in 2:N.spp) {core.mm[i,j] ~ dnorm(0, 1.0E-3)}
}
#plot level priors
for(i in 1:N.plot.preds){
plot.m[i,1] <- 0
for (j in 2:N.spp) {plot.m[i,j] ~ dnorm(0, 1.0E-3)}
}
#site level priors
for(i in 1:N.site.preds){
site.m[i,1] <- 0
for (j in 2:N.spp) {site.m[i,j] ~ dnorm(0, 1.0E-3)}
}

#predictor (x) values drawn from distributions.
for(j in 1:N.core.preds){for(i in 1:N.core){core[i,j] ~ dnorm(core_mu[i,j], core_precision[i,j])}} #core x values
for(j in 1:N.plot.preds){for(i in 1:N.plot){plot[i,j] ~ dnorm(plot_mu[i,j], plot_precision[i,j])}} #plot x values
for(j in 1:N.site.preds){for(i in 1:N.site){site[i,j] ~ dnorm(site_mu[i,j], site_precision[i,j])}} #site x values

#mean center all predictors (except core-level intercept).
for(i in 1:N.core){
core.center[i,1] <- 1
for(j in 2:N.core.preds){
core.center[i,j] <- core[i,j] - mean(core[,j])
}
}
for(i in 1:N.plot){for(j in 1:N.plot.preds){plot.center[i,j] <- plot[i,j] - mean(plot[,j])}}
for(i in 1:N.site){for(j in 1:N.site.preds){site.center[i,j] <- site[i,j] - mean(site[,j])}}

#save mean values for back transforming intercept values.
for(j in 1:N.core.preds){core.center.save[j] <- mean(core[,j])}
for(j in 1:N.plot.preds){plot.center.save[j] <- mean(plot[,j])}
for(j in 1:N.site.preds){site.center.save[j] <- mean(site[,j])}

#combine core/plot/site (x) values and parameters.
for(j in 1:N.spp){
for(i in 1:N.core){core.level[i,j] <- inprod(core.mm[,j], core.center[i,])           } #core level
for(p in 1:N.plot){plot.level[p,j] <- inprod(plot.m [,j], plot.center[plot_plot[p],])} #plot level
for(s in 1:N.site){site.level[s,j] <- inprod(site.m [,j], site.center[site_site[s],])} #site level
}

#fit species abundances as a linear combination of core, plot and site predictors.
for(i in 1:N.core){
for(j in 1:N.spp){
log(a0[i,j]) <- alpha + core.level[i,j] + plot.level[core_plot[i],j] + site.level[core_site[i],j]
#log(a0[i,j]) <- alpha + core.level[i,j]
}
y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp]) 
}

#map to original parameterization, assuming first column of core level predictors is intercept.
for(j in 1:N.spp){
core.norm[j] <- inprod(core.mm[2:N.core.preds,j], core.center.save[2:N.core.preds])
plot.norm[j] <- inprod(plot.m[,j], plot.center.save)
site.norm[j] <- inprod(site.m[,j], site.center.save)
}

for (j in 1:N.spp) {
#core.m[1,j] <- alpha + core.mm[1,j]
core.m[1,j] <- alpha + core.mm[1,j] - (core.norm[j] + plot.norm[j] + site.norm[j])
for (i in 2:N.core.preds){
core.m[i,j] <- core.mm[i,j]
}
}

} #close model loop.
"
###Fit JAGS model.
#parallel or not parallel.
run.method <- ifelse(parallel == F,'rjags','rjparallel')
#run jags model.
jags.out <- runjags::run.jags(   model = jags.model,
                                 data = jags.data,
                                 adapt = adapt,
                                 burnin = burnin,
                                 sample = sample,
                                 n.chains = n.chains,
                                 method = run.method,
                                 #monitor = c('alpha','core.mm','core.m','plot.m','site.m'))
                                 monitor = c('core.m','plot.m','site.m'))
#summarize output
out <- summary(jags.out)

#grab parmeters by species, make a list of species-parameter dataframes
output.list <- list()
for(i in 1:ncol(y)){
  z <- out[grep(paste0(',',i,']'),rownames(out)),]
  names <- c(core.names,plot.names,site.names)
  z <- cbind(names,z)
  colnames(z)[1] <- 'predictor'
  output.list[[i]] <- data.frame(z)
}
names(output.list) <- y.names

#get the matrix of predicted y values.
core.core <- data.frame(cbind(core_plot,core_site,core_mu))
plot.plot <- data.frame(cbind(plot_plot,plot_mu))
site.site <- data.frame(cbind(site_site,site_mu))
super.x <- merge(core.core,plot.plot, by.x = 'core_plot', by.y = 'plot_plot', all.x = T)
super.x <- merge(  super.x,site.site, by.x = 'core_site', by.y = 'site_site', all.x = T)
super.x <- super.x[,3:ncol(super.x)]

pred.list <- list()
for(i in 1:ncol(y)){
  pred <- exp(as.matrix(super.x) %*% as.numeric(as.character(output.list[[i]][,5])))
  pred.list[[i]] <- pred
}
predicted <- do.call('cbind',pred.list)
predicted <- predicted / rowSums(predicted)
colnames(predicted) <- colnames(y)
#get matrix of residuals
resid <- y - predicted

#Plot predicted vs. observed and the 1:1 line.
#regression (green) and 1:1 line (black) plot on top of eachother. YES.
par(mfrow = c(1,ncol(y)))
for(i in 1:ncol(y)){
  lab1 <- paste0('predicted y',i)
  lab2 <- paste0('observed y',i)
  plot(y[,i] ~ predicted[,i], xlab = lab1, ylab = lab2)
  mod <- lm(y[,i] ~ predicted[,i])
  r.sq <- summary(mod)$r.squared
  abline(0,1, lwd = 2) 
  abline(mod, lwd = 2, lty = 2, col = 'green')
  txt <- paste0('R2 = ',round(r.sq,2))
  mtext(txt, side = 3, line = -2, adj = 0.05)
}
