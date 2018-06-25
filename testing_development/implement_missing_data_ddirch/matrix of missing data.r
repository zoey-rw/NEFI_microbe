#This takes a matrix of x values, fits c/p/s hierarchically.
#clear environment, load runjags.
rm(list = ls())
library(runjags)

#10 sites, 10 plots per site, 10 observations per plot. 
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
x1 <- unlist(obs)
obs <- list()
for(i in 1:length(plot.means)){
  obs[[i]] <- rnorm(n.core, mean = plot.means[i], sd = 2)
}
x2 <- unlist(obs)

#create true data frame.
dat <- data.frame(siteID,plotID)
dat <- dat[order(plotID),]
dat$x1 <- x1
dat$x2 <- x2
X <- as.matrix(dat[,c('x1','x2')])

#indexing.
plotindex <- as.factor(dat$plotID)
siteindex <- list()
for(i in 1:n.site){
  siteindex[[i]] <- rep(LETTERS[i],n.plot)
}
siteindex <- as.factor(unlist(siteindex))

jags.model = "
model {
  #priors
  for(m in 1:n.X){
      global[m] ~ dnorm(0,1.0E-4) #global level parameter prior.
    site.tau[m] ~ dgamma(0.1,0.1)
    plot.tau[m] ~ dgamma(0.1,0.1)
    core.tau[m] ~ dgamma(0.1,0.1)
    }
  
  #estimate means
  for(m in 1:n.X){
    for(i in 1:N.site){ site[i,m] ~ dnorm(global[m], site.tau[m])}           ## estimate global mean
    for(j in 1:N.plot){ plot[j,m] ~ dnorm(site[siteindex[j],m],plot.tau[m])} ## estimate site means
    for(k in 1:N.core){ core[k,m] ~ dnorm(plot[plotindex[k],m],core.tau[m])} ## estimate plot means
  }

} #end model
"
#specify data object.
jags.data <- list(core = X, n.X = ncol(X), 
                  N.site = length(unique(dat$siteID)),
                  N.plot = length(unique(dat$plotID)),
                  N.core = nrow(dat),
                  plotindex = plotindex, siteindex=siteindex)
#fit model.
jags.out <- runjags::run.jags(jags.model,
                              data = jags.data,
                              n.chains = 3,
                               monitor = c('plot','site'))
out <- summary(jags.out)
#compare estimate plot/site means to true values. All w/in 95% credible interval.
par(mfrow = c(2,2))
fit.plot.means <- out[grep('plot',rownames(out)),][,4]
fit.site.means <- out[grep('site',rownames(out)),][,4]
plot(plot.means ~ fit.plot.means[1:90], main = 'plot-level fit x1'); abline(0,1,lwd=2, lty=2)
plot(plot.means ~ fit.plot.means[91:180], main = 'plot-level fit x2'); abline(0,1,lwd=2, lty=2)
plot(site.means ~ fit.site.means[1:10], main = 'site-level fit x1'); abline(0,1,lwd=2, lty=2)
plot(site.means ~ fit.site.means[11:20], main = 'site-level fit x2'); abline(0,1,lwd=2, lty=2)

#check uncertainty values
jags.out <- runjags::run.jags(jags.model,
                              data = jags.data,
                              n.chains = 3,
                              monitor = c('core.tau','plot.tau','site.tau','global'))
summary(jags.out)


#create data with missing observations.
missing <- dat
#remove 3 observations within a plot in site A.
missing[2:4,3] <- NA
missing[4:7,4] <- NA
#remove an entire plot in site B.
missing[missing$plotID == 'B_1',3] <- NA
missing[missing$plotID == 'B_2',4] <- NA
#remove all observations in site C.
missing[missing$siteID == 'C'  ,3] <- NA
missing[missing$siteID == 'B'  ,4] <- NA

X <- as.matrix(missing[,3:4])

#creat new jags data object.
jags.data <- list(core = X, n.X = ncol(X), 
                  N.site = length(unique(dat$siteID)),
                  N.plot = length(unique(dat$plotID)),
                  N.core = nrow(dat),
                  plotindex = plotindex, siteindex=siteindex)

#fit, see what the model estimates for core level values.
jags.out <- runjags::run.jags(jags.model,
                              data = jags.data,
                              n.chains = 3,
                              monitor = c('core','plot','site'))

out <- summary(jags.out)
comp <- dat


###FROM here on has not been modded for multivariate case.
comp$pred <- out[grep('core',rownames(out)),4]
comp$removed <- missing$x
z <- comp[is.na(comp$removed),]

#how are the core level predcitions compared to true values? in range.
cbind(out[2:4,],dat$x[2:4])

#how about the missing plot?
zz <- cbind(out[grep('core',rownames(out)),],missing,dat$x)
zz[zz$plotID == 'B_1',]

#how about the missing site?
zz[zz$siteID == 'C',]
