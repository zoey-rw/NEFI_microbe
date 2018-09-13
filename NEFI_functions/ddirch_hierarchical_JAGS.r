#' hierarchical_dirlichet_jags()
#' Fits a hierarchical dirlichet model with x-values observed at the core, plot and site level.
#'
#' @param y species matrix as a dataframe. Column names are species names.
#' @param core_mu core level observation means as a dataframe. column names are predictor names.
#' @param plot_mu plot level observation means as a dataframe. column names are predictor names.
#' @param site_mu site level observation means as a dataframe. column names are predictor names.
#' @param core_plot factor vector that links core level observations to plot IDs. Must be a factor.
#' @param core_site factor vector that links core level observations to site IDs. Must be a factor.
#' @param plot_site factor vector that links plot level observations to site IDs. Must be a factor.
#' @param plot_plot factor vector of unique plot IDs, same order as plot data. Must be a factor.
#' @param site_site factor vector od unique site IDs, same order as site data. Must be a factor.
#' @param core_sd core level observation sd as a dataframe. column names should match core_mu. optional.
#' @param plot_sd plot level observation sd as a dataframe. column names should match plot_mu. optional.
#' @param site_sd site level observation sd as a dataframe. column names should match site_mu. optional.
#' @param adapt   number of adaptive iterations for JAGS simulation.
#' @param burnin  number of burnin iterations for JAGS simulation.
#' @param sample  number of sample iterations for JAGS simulation.
#' @param n.chains number of chains for JAGS simulation.
#' @param parallel whether or not to run JAGS chains in parallel.
#'
#' @return returns a list with the fitted model, a list of species by parameter tables, and matrices of predicted, observed and residual values on the observation scale (0,1)
#' @export
#'
#' @examples
#' adapt = 50; burnin = 100; sample = 200; n.chains = 3; parallel = F;
hierarchical_dirlichet_jags   <- function(y,
                                          core_mu, plot_mu, site_mu,
                                          core_sd = NA, plot_sd = NA, site_sd = NA,
                                          core_plot, core_site, plot_site, plot_plot, site_site,
                                          adapt = 500, burnin = 1000, sample = 2000, n.chains = 3, parallel = F){
  #Load some important dependencies.
  source('/home/caverill/NEFI_microbe/NEFI_functions/crib_fun.r')
  source('/home/caverill/NEFI_microbe/NEFI_functions/sd_to_precision.r')
  source('/home/caverill/NEFI_microbe/NEFI_functions/precision_matrix_match.r')
  
  #Some checks before we get started.
  y    <- as.data.frame(y)
  core_mu <- as.data.frame(core_mu)
  plot_mu <- as.data.frame(plot_mu)
  site_mu <- as.data.frame(site_mu)
  #core level. First column needs to be intercept, a vector of 1s.
  if(mean(core_mu[,1]) != 1){stop('First column in core_mu is not a vector of 1s for the intercept. This needs to be or this function doesnt work.')}
  if(ncol(core_mu) < 2){stop('We need an intercept and at least one predictor at the core level for this to work. core_mu has less than 2 columns. Try again buddy.')}
  
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
  
  #Match up predictors and their SD. if no SD supplied we assign perfect (really really high) precision.
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
  plot_site <- droplevels(as.factor(plot_site))
  site_site <- droplevels(as.factor(site_site))
  
  ###setup jags data object.
  jags.data <- list(N.core = nrow(core_mu), N.plot = nrow(plot_mu), N.site = nrow(site_mu), N.spp = ncol(y),            #number of c/p/s observations and number of species
                    N.core.preds = ncol(core_mu), N.plot.preds = ncol(plot_mu), N.site.preds = ncol(site_mu),           #number of predictors at c/p/s level
                    core_mu = core_mu, plot_mu = plot_mu, site_mu = site_mu,                                            #x-value mean      matrices c/p/s
                    core_precision = core_precision, plot_precision = plot_precision, site_precision = site_precision,  #x-value precision matrices c/p/s
                    core_plot = core_plot, core_site = core_site,                                                       #vectors that index across  c/p/s
                    plot_plot = plot_plot, plot_site = plot_site, site_site = site_site,                                #vectors that index across  c/p/s
                    y = y)                                                                                              #species matrix, y
  
  ###specify JAGS model.
  jags.model = "
  model {
  #priors for each species*predictor*hierarchical-level.----
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

  ### Begin missing data model ###
  #missing X data priors, core-level.
  for(m in 1:N.core.preds){
      core.global[m] ~ dnorm(0,1.0E-4) #global level parameter prior.
    core.site.tau[m] ~ dgamma(0.1,0.1)
    core.plot.tau[m] ~ dgamma(0.1,0.1)
    core.core.tau[m] ~ dgamma(0.1,0.1)
  }
  #missing X data priors, plot-level.
  for(m in 1:N.plot.preds){
      plot.global[m] ~ dnorm(0,1.0E-4) #global level parameter prior.
    plot.site.tau[m] ~ dgamma(0.1,0.1)
    plot.plot.tau[m] ~ dgamma(0.1,0.1)
  }
  #missing X data priors, site-level.
  for(m in 1:N.site.preds){
      site.global[m] ~ dnorm(0,1.0E-4) #global level parameter prior.
    site.site.tau[m] ~ dgamma(0.1,0.1)
  }
  
  #fill in any missing X values at core level.
  for(m in 1:N.core.preds){
    for(i in 1:N.site){ site.core[i,m] ~ dnorm(core.global[m], core.site.tau[m])}           ## estimate global mean
    for(j in 1:N.plot){ plot.core[j,m] ~ dnorm(site.core[plot_site[j],m],core.plot.tau[m])} ## estimate site means
    for(k in 1:N.core){   core_mu[k,m] ~ dnorm(plot.core[core_plot[k],m],core.core.tau[m])} ## estimate plot means
  }
  #fill in any missing X values at plot level.
    for(m in 1:N.plot.preds){
    for(i in 1:N.site){ site.plot[i,m] ~ dnorm(plot.global[m], plot.site.tau[m])}           ## estimate global mean
    for(j in 1:N.plot){   plot_mu[j,m] ~ dnorm(site.plot[plot_site[j],m],plot.plot.tau[m])} ## estimate site means
  }
  #fill in any missing X values at site level.
  for(m in 1:N.plot.preds){
    for(i in 1:N.site){   site_mu[i,m] ~ dnorm(site.global[m], site.site.tau[m])}           ## estimate global mean
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
  
  #make a super output that also returns model
  super.list <- list(jags.out, output.list,predicted,y,resid)
  names(super.list) <- c('jags_model','species_parameter_output','predicted','observed','residual')
  
  #return model and output
  return(super.list)
  
}