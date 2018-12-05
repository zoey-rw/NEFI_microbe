#' hierarch_ddirch_means.r
#' Gets plot and site level means for a species matrix.
#' Fits a plot level model, then uses the output of the plot model for the site level model.
#' ideally this would be one model, but that thing won't converge.
#'
#' @param y           #species matrix.
#' @param core_plot   #vector linking core-level observations to plots.
#' @param plot_site   #vector linking plots to sites.
#'
#' @return            #plot and site mean matrices, as well as upper and lower 95% matrices.
#' @export
#'
#' @examples
hierarch_ddirch_means <- function(y, core_plot, plot_site, jags.method = 'rjags'){
  #specify plot data object.----
  jd.plot <- list(y=as.matrix(y), N.core = nrow(y), N.plot=length(unique(core_plot)), 
                  N.spp = ncol(y), core_plot = as.factor(core_plot))
  
  #jags model 1: plot level.----
  jags.model1 = "
  model {
    #plot means.
    for(i in 1:N.core){
      for(j in 1:N.spp){
        log(core.hat[i,j]) <- plot_mu[core_plot[i],j] + plot.intercept
      }
      y[i,1:N.spp] ~ ddirch(core.hat[i,1:N.spp]) 
    }
    
    #priors.
    plot.intercept ~ dnorm(0,1E-3)
    for(i in 1:N.plot){
      plot_mu[i,1] <- 0
      for(j in 2:N.spp){
        plot_mu[i,j] ~ dnorm(0,1E-3)
      }
    }
    #Map back for easy handling downstream.
    for(i in 1:N.plot){
      for(j in 1:N.spp){
        log(plot_pre[i,j]) <- plot_mu[i,j] + plot.intercept  
      }
      plot_out[i,1:N.spp] <- plot_pre[i,1:N.spp]/sum(plot_pre[i,1:N.spp])
    }
  }" #end jags model.

  #jags model 2: site level.----
  jags.model2 = "
  model {
    #site means.
    for(i in 1:N.plot){
    for(j in 1:N.spp){
    log(plot.hat[i,j]) <- site_mu[plot_site[i],j] + site.intercept
    }
    y[i,1:N.spp] ~ ddirch(plot.hat[i,1:N.spp]) 
    }
    
    #priors.
    site.intercept ~ dnorm(0,1E-3)
    for(i in 1:N.site){
      site_mu[i,1] <- 0
      for(j in 2:N.spp){
        site_mu[i,j] ~ dnorm(0,1E-3)
      }
    }
    #Map back for easy handling downstream.
    for(i in 1:N.site){
      for(j in 1:N.spp){
        log(site_pre[i,j]) <- site_mu[i,j] + site.intercept  
      }
      site_out[i,1:N.spp] <- site_pre[i,1:N.spp]/sum(site_pre[i,1:N.spp])
    }
  }" #end jags model.
  
  #fit plot level.----
  cat('Fitting plot level model...\n')
  pmod <- run.jags(jags.model1,
                   data = jd.plot,
                   adapt = 1500,
                   burnin = 5000,
                   sample = 3000,
                   n.chains = 3,
                   method = jags.method,
                   monitor = c('plot_out'))
  cat('plot level model fit.\n')
  #summarize plot fit.----
  plot.out <- summary(pmod, vars = 'plot_out')
  plot.mean.matrix <- matrix(plot.out[,4], nrow = length(plot_site), ncol = ncol(y))
  plot.lo95.matrix <- matrix(plot.out[,1], nrow = length(plot_site), ncol = ncol(y))
  plot.hi95.matrix <- matrix(plot.out[,3], nrow = length(plot_site), ncol = ncol(y))
  plot.fit <- list(plot.mean.matrix,plot.lo95.matrix,plot.hi95.matrix)
  names(plot.fit) <- c('mean','lo95','hi95')
  
  #specify site data object.----
  jd.site <- list(y = plot.mean.matrix, N.plot=length(unique(core_plot)),
                  N.spp = ncol(y), plot_site = as.factor(plot_site), 
                  N.site = length(unique(plot_site)))
  
  #fit site level.----
  cat('Fitting site level model...\n')
  smod <- run.jags(jags.model2,
                   data = jd.site,
                   adapt = 1500,
                   burnin = 5000,
                   sample = 3000,
                   n.chains = 3,
                   method = jags.method,
                   monitor = c('site_out'))
  cat('site level model fit.\n')
  
  #summarize site fit.----
  site.out <- summary(smod, var = 'site_out')
  site.mean.matrix <- matrix(site.out[,4], nrow = length(unique(plot_site)), ncol = ncol(y))
  site.lo95.matrix <- matrix(site.out[,1], nrow = length(unique(plot_site)), ncol = ncol(y))
  site.hi95.matrix <- matrix(site.out[,3], nrow = length(unique(plot_site)), ncol = ncol(y))
  site.fit <- list(site.mean.matrix,site.lo95.matrix,site.hi95.matrix)
  names(site.fit) <- c('mean','lo95','hi95')
  
  #return plot and site fits.----
  out <- list(plot.fit, site.fit)
  names(out) <- c('plot.fit','site.fit')
  return(out)
  cat('fitting complete!\n')
}