#' hierarch_plot.means_JAGS
#' Get site level means and sd from NEON data observed at the plot scale.
#' Also return global mean and sd for a given predictor.
#' depends on the runjags package.
#'
#' @param x_mu         #vector of observations to modeled.
#' @param plot_site    #NEON plot IDs associated with x_mu observations.
#' @param n.adapt      #number of adaptive iterations. default 200.
#' @param n.burnin     #number of burnin iterations. default 1000.
#' @param n.sample     #number of sample iterations. default 1000.
#' @param n.chains     #number of MCMC chains. default 3.
#' @param parallel     #run chains in parallel or not? default F, notsetup to work yet.
#'
#' @return
#' @export
#'
#' @examples
hierarch_plot.means_JAGS <- function(x_mu, plot_site,
                                     n.adapt = 200, n.burnin = 1000, n.sample = 1000, n.chains = 3, parallel = F){
  
  #### specify JAGS model ####
  jags.model = "
  model {
  #get site level means.
  for(i in 1:N.plot) {
  plot_mu[i] ~ dnorm(site_mu[plot_site[i]], tau.site[plot_site[i]])
  }
  
  #get global mean and uncertainty.
  for(i in 1:N.site){
  site_mu[i] ~ dnorm(global_mu, tau.glob)
  }
  
  #Site level priors
  for(i in 1:N.site) {
  sigma.site[i] ~ dunif(0, 100)
  tau.site[i] <- pow(sigma.site[i], -2)
  }

  #Global level priors
  global_mu ~ dnorm(0,1E-3)  #I(0, ) #interval censored to be greater than zero.
  sigma.glob ~ dunif(0,100)
  tau.glob <- pow(sigma.glob, -2)
  
  }"
  
  #### setup JAGS data object. ####
  dat <- data.frame(x_mu,plot_site)
  dat <- dat[complete.cases(dat),]
  site.names <- unique((plot_site))
  
  jags.data <- list(N.plot = length(plot_site), N.site = length(unique(plot_site)),
                    plot_mu = dat$x_mu, 
                    plot_site = droplevels(as.factor(dat$plot_site))
                    )
  
  #### fit JAGS model ####
  #runmode <- ifelse(parallel == T,T,F)
  mod  <- runjags::run.jags(model = jags.model,
                            data = jags.data,
                            monitor = c('site_mu','global_mu'),
                            adapt = n.adapt,
                            burnin = n.burnin,
                            sample = n.sample,
                            n.chains = n.chains)
  out <- summary(mod)
  site.table <- data.frame(out[grep('site_mu', rownames(out)),])
  site.table$siteID <- site.names
  glob.table <- data.frame(t(out[grep('global_mu', rownames(out)),]))
  
  
  #### return output ####
  output.list <- list(out,site.table,glob.table)
  names(output.list) <- c('jags.summary','site.table','glob.table')
  return(output.list)
  
}