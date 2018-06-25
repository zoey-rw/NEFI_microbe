#clear environment, load runjags.
rm(list = ls())
library(runjags)

#### Generate pseudo data ####
#10 sites, 9 plots per site, 10 observations per plot. 
n.site <- 10
n.plot <- 9
n.core <- 10
siteID <- sort(rep(LETTERS[1:n.site],n.plot*n.core))
plotID <- rep(1:n.plot,n.site*n.core)
plotID <-paste(siteID,plotID,sep = '_')
site.means <- runif(n.site,10,25)
plot.means <- list()
for(i in 1:length(site.means)){
  plot.means[[i]] <- rnorm(n.plot,site.means[i],2)
}
plot.means <- unlist(plot.means)

#get observations for each plot, based on a mean, and a sd of 2
obs <- list()
for(i in 1:length(plot.means)){
  obs[[i]] <- rnorm(n.core, mean = plot.means[i], sd = 2)
}
x <- unlist(obs)

#create true data frame.
dat <- data.frame(siteID,plotID)
dat <- dat[order(plotID),]
dat$x <- x

#indexing.
plotindex <- as.factor(dat$plotID)
siteindex <- list()
for(i in 1:n.site){
  siteindex[[i]] <- rep(LETTERS[i],n.plot)
}
siteindex <- as.factor(unlist(siteindex))


#### Specify hierarchical JAGS model ####
jags.model = "
model {
#priors
global ~ dnorm(0,1.0E-4) #global level parameter prior.
site.tau ~ dgamma(0.1,0.1)
plot.tau ~ dgamma(0.1,0.1)
core.tau ~ dgamma(0.1,0.1)

for(i in 1:N.site){ site[i] ~ dnorm(global, site.tau)}            ## estimate global mean
for(j in 1:N.plot){ plot[j] ~ dnorm(site[siteindex[j]],plot.tau)} ## estimate site means
for(k in 1:N.core){ core[k] ~ dnorm(plot[plotindex[k]],core.tau)} ## estimate plot means

} #end model
"
##### specify JAGS data object. ####
jags.data <- list(core = dat$x, 
                  N.site = length(unique(dat$siteID)),
                  N.plot = length(unique(dat$plotID)),
                  N.core = nrow(dat),
                  plotindex = plotindex, siteindex=siteindex)

#### Fit JAGS model ####
jags.out <- runjags::run.jags(jags.model,
                              data = jags.data,
                              n.chains = 3,
                              monitor = c('plot','site'))
out <- summary(jags.out)

##### compare estimate plot/site means to true values. All w/in 95% credible interval. ####
par(mfrow = c(1,2))
fit.plot.means <- out[grep('plot',rownames(out)),][,4]
fit.site.means <- out[grep('site',rownames(out)),][,4]
plot(plot.means ~ fit.plot.means, main = 'plot-level fit'); abline(0,1,lwd=2, lty=2)
plot(site.means ~ fit.site.means, main = 'site-level fit'); abline(0,1,lwd=2, lty=2)

##### check uncertainty values #####
jags.out <- runjags::run.jags(jags.model,
                              data = jags.data,
                              n.chains = 3,
                              monitor = c('core.tau','plot.tau','site.tau','global'))
summary(jags.out)


##### create data with missing observations. #####
missing <- dat
#remove 3 observations within a plot in site A.
missing[2:4,3] <- NA
#remove an entire plot in site B.
missing[missing$plotID == 'B_1',3] <- NA
#remove all observations in site C.
missing[missing$siteID == 'C'  ,3] <- NA

##### create new jags data object for fitting missing data model. ####
jags.data <- list(core = missing$x, 
                  N.site = length(unique(dat$siteID)),
                  N.plot = length(unique(dat$plotID)),
                  N.core = nrow(dat),
                  plotindex = plotindex, siteindex=siteindex)

##### Fit missing data model. ####
jags.out <- runjags::run.jags(jags.model,
                              data = jags.data,
                              n.chains = 3,
                              monitor = c('core','plot','site'))

out <- summary(jags.out)
comp <- dat
comp$pred <- out[grep('core',rownames(out)),4]
comp$removed <- missing$x
z <- comp[is.na(comp$removed),]

##### how are the core level predcitions compared to true values? in range. #####
cbind(out[2:4,],dat$x[2:4])

##### how about the missing plot? ####
zz <- cbind(out[grep('core',rownames(out)),],missing,dat$x)
zz[zz$plotID == 'B_1',]

##### how about the missing site? ####
zz[zz$siteID == 'C',]
