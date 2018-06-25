#missing data model.
#fits a global level predictor for each independent variable.
#global level predictor is site level prior.
#site level estimate is plot level prior.
#plot level estimate is core level prior.
#if core level estimate is absent it is drawn from site level estimate.

#load some core level data missing at core and plot scale.
rm(list=ls())

#load data 
d <- readRDS('/fs/data3/caverill/NEFI_microbial/NEON_data_aggregation/core.table.rds')
d <- d[d$siteID %in% c('ORNL','DSNY'),]           #subset to two sites for development.
y <- as.matrix(d[,c('soilTemp','soilInWaterpH')]) #matrix of x values.

#drop factor levels. Important when subsetting to a few sites.
d$siteID <-droplevels(d$siteID)
d$plotID <-droplevels(d$plotID)

#setup categorical predictor matrices.
plot.preds <- model.matrix( ~ plotID - 1, data = d)
site.preds <- model.matrix( ~ siteID - 1, data = d)

#indexing.
plot_plot <- as.factor(unique(d$plotID))
plot_site <- as.factor(substr(plot_plot,1,4))

jags.model = "
model {
  #priors
  for(j in 1:N.x){
  glob[j] ~ dnorm(0,1.0E-4) #global level parameter prior.
    for(k in 1:n.site.preds){site[k,j] <- glob[j]}              #prior on site level parameters is glob level parameter.
    for(k in 1:n.plot.preds){plot[k,j] <- site[plot_site[k],j]} #prior on plot level parameters is site level parameter.
     sigma[j] ~ dunif(0, 100)
  glob.tau[j] <- pow(sigma[j], -2)
  site.tau[j] <- pow(sigma[j], -2)
  plot.tau[j] <- pow(sigma[j], -2)
  core.tau[j] <- pow(sigma[j], -2)
  }
  
  #combine priors and parameters
  for(j in 1:N.x){
    for(i in 1:N){
      glob.x[i,j] ~ dnorm(glob[j],glob.tau[j])
    site.hat[i,j] <- inprod(site[,j], site.preds[i,])
    plot.hat[i,j] <- inprod(plot[,j], plot.preds[i,])
      site.x[i,j] ~ dnorm(site.hat[i,j], site.tau[j])
      plot.x[i,j] ~ dnorm(plot.hat[i,j], plot.tau[j]) 
      core.x[i,j] ~ dnorm(  plot.x[i,j], core.tau[j]) #predict missing core values based on plot estimate.
    }
  }

} #end model
"
jags.data <- list(core.x = y, 
                  plot.preds=plot.preds, site.preds=site.preds, 
                  plot_site = plot_site,
                  n.plot.preds = ncol(plot.preds), n.site.preds = ncol(site.preds), 
                  N.x = ncol(y), N = nrow(d))

jags.out <- runjags::run.jags(jags.model,
                              data = jags.data,
                              n.chains = 3,
                              monitor = c('plot','site','core.x'))

#check estiamted values of temp and pH.
out <- summary(jags.out)
core.vals <- out[grep('core.x',rownames(out)),]
temp.new <- core.vals[grep(',1]',rownames(core.vals)),]
  pH.new <- core.vals[grep(',2]',rownames(core.vals)),]
temp <- data.frame(d$soilTemp,d$plotID,temp.new[,1:4])
  pH <- data.frame(d$soilInWaterpH,d$plotID,  pH.new[,1:4])
