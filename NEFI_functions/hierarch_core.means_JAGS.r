#' hierarch_core.means_JAGS
#' Get site level means and sd from NEON data observed at the core scale.
#' Also return global mean and sd for a given predictor.
#' depends on the runjags package.
#'
#' @param x_mu         #vector of observations to modeled.
#' @param core_plot    #NEON plot IDs associated with x_mu observations.
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
hierarch_core.means_JAGS <- function(x_mu, core_plot,
                                n.adapt = 200, n.burnin = 1000, n.sample = 1000, n.chains = 3, parallel = F){
  
  #### specify JAGS model ####
  jags.model = "
  model {
  #get plot level means
  for(i in 1:N.core){
  core_mu[i] ~ dnorm(plot_mu[core_plot[i]], tau.plot[core_plot[i]])
  }
  
  #get site level means.
  for(i in 1:N.plot) {
  plot_mu[i] ~ dnorm(site_mu[plot_site[i]], tau.site[plot_site[i]])
  }
  
  #get global mean and uncertainty.
  for(i in 1:N.site){
  site_mu[i] ~ dnorm(global_mu, tau.glob)
  }
  
  #Plot level priors
  for(i in 1:N.plot){
  sigma.plot[i] ~ dunif(0, 100)
  tau.plot[i] <- pow(sigma.plot[i], -2)
  }
  #Site level priors
  for(i in 1:N.site) {
  sigma.site[i] ~ dunif(0, 100)
  tau.site[i] <- pow(sigma.site[i], -2)
  }
  #Global level priors
  global_mu ~ dnorm(0,1E-3)
  sigma.glob ~ dunif(0,100)
  tau.glob <- pow(sigma.glob, -2)
  
  }"
  
  #### setup JAGS data object. ####
  dat <- data.frame(x_mu,core_plot)
  dat <- dat[complete.cases(dat),]
  plot_site <- substring(unique(dat$core_plot),1,4)
  site.names <- unique((plot_site))
  
  jags.data <- list(N.core = nrow(dat), N.plot = length(plot_site), N.site = length(unique(plot_site)),
                    core_mu = dat$x_mu, core_plot = droplevels(as.factor(dat$core_plot)), plot_site = droplevels(as.factor(plot_site)))

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