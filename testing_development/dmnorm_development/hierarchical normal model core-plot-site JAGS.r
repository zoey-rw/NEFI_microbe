#simulate some data with just core, plot and site intercepts.
#2 sites, 2 plots per site, 10 cores per plot = 40 obs
rm(list = ls())

mat <- c(rep(2,10),rep(4,10), rep(6,10), rep(8,10))
map <- c(rep(400,20), rep(200,20))
core <- (data.frame(mat,map))
core$y <- mat*0.5 + map*-0.01 + rnorm(40,2)
core$plot <- c(rep('a1',10),rep('a2',10), rep('b1',10), rep('b2',10))
core$site <- c(rep('a',20), rep('b',20))
core$intercept <- rep(1, nrow(core))
core_plot <- as.factor(core$plot)
core_site <- as.factor(core$site)
plot_plot <- as.factor(unique(core$plot))
site_site <- as.factor(unique(core$site))
y <- as.matrix(core$y)

site <- as.matrix(c(400,200))
plot <- as.matrix(c(2,4,6,8))
#mat is at the plot level
#map is at the site level
summary(lm(y ~ mat + map, data = core))

#specify JAGS model. This is general for any linear combination of x's.
jags.model = "
model {
###priors for each species*predictor*hierarchical-level
#y uncertainty prior
tau <- pow(sigma, -2)
sigma ~ dunif(0, 100)

for(j in 1:N.spp){
for(a in 1:N.core.preds){core.m[a,j] ~ dnorm(0, 1.0E-4)} #core level priors
for(b in 1:N.plot.preds){plot.m[b,j] ~ dnorm(0, 1.0E-4)} #plot level priors
for(c in 1:N.site.preds){site.m[c,j] ~ dnorm(0, 1.0E-4)} #site level priors
}

###predictor (x) values drawn from distributions.
#core values
#for(j in 1:N.core.preds){
#for(i in 1:N.core){
#core[i,j] ~ dnorm(core_mu[i,j], core_precision[i,j])
#}
#}
#plot values
#for(j in 1:N.plot.preds){
#for(i in 1:N.plot){
#plot[i,j] ~ dnorm(plot_mu[i,j], plot_precision[i,j])
#}
#}
#site values
#for(j in 1:N.site.preds){
#for(i in 1:N.site){
#site[i,j] ~ dnorm(site_mu[i,j], site_precision[i,j])
#}
#}

#combine predictors and parameters.
for(j in 1:N.spp){
for(i in 1:N.core){core.level[i,j] <- inprod(core.m[,j], core[i,])           } #core level
for(p in 1:N.plot){plot.level[p,j] <- inprod(plot.m[,j], plot[plot_plot[p],])} #plot level
for(s in 1:N.site){site.level[s,j] <- inprod(site.m[,j], site[site_site[s],])} #site level
}


#go ahead and fit means of species abundances as a linear combination of predictors and parameters.
for(i in 1:N.core){
    #log link a0 for dirlichet
    for(j in 1:N.spp){
        a0[i,j] <- core.level[i,j] + plot.level[core_plot[i],j] + site.level[core_site[i],j]
    }
  #y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp]) 
   y[i,1:N.spp] ~dnorm(a0[i,1:N.spp], tau)
}

} #close model loop.
"

jags.data <- list(N.core = nrow(core), N.plot = length(plot_plot), N.site = length(site_site), N.spp = 1,
                  N.core.preds = 1, N.plot.preds = 1, N.site.preds = 1,
                  core = as.matrix(core$intercept), plot = plot, site = site,
                  core_plot = core_plot, core_site = core_site, plot_plot = plot_plot, site_site = site_site,
                  y = y)
                  #core_precision = core_precision, plot_precision = plot_precision, site_precision = site_precision,


jags.out <- run.jags(jags.model,
                     data=jags.data,
                     adapt = 200,
                     burnin = 1000,
                     sample = 2000,
                     n.chains=3,
                     method = 'rjags',
                     monitor = c('core.m','plot.m','site.m'))
summary(jags.out)
