#core-level forecast to NEON.
#once data have been aggregated at core-plot-site level, we need to fill in missing values and uncertainties hierarchically.
#this product is then shipped to the forecast, which draws missing data.
#1. hierarchically execute your site-level forecast, it should basically be the same, but at the core level.
#2. No need to update yet, Just get hierarchical means.
#3. You in principle can use this as your plot and site level forecast? Once you figure out hierarchical ddirch means.
#4. The relEM problem- you solve this in forecast because missing data happens outside of JAGS.
#5. Once you get this done, hope ddirch aggregation is solved. If not, ask Mike. THen we can validate.
#6. Replication *should* be quick on the bacteria side. Clustered, binned taxonomically.
#7. Fit prior, then send prior through core-level forecast, aggregate plot and site level using something.
#clearn environment, loads paths.
rm(list=ls())
source('paths.r')
source('NEFI_functions/precision_matrix_match.r')

#1. data prep.----
#load prior model fit- model fit at site level.
mod <- readRDS(ted_ITS.prior_fg_JAGSfit)
mod <- mod$all.preds
#how many species (taxonomic or functional groups)
n.spp <- length(mod$species_parameter_output)
#merge mcmc output
mcmc <- mod$jags_model$mcmc
mcmc <- do.call(rbind, mcmc)
mcmc <- mcmc[,grep("^x\\.m\\[",colnames(mcmc))]

#predictors.
preds <- mod$species_parameter_output$other$predictor
indices <- c('sampleID','plotID','siteID')

#load NEON data, missing values have been hierarchically estimated and assigned uncertainties.
dat <- readRDS(hierarch_filled.path)

#choose your data products. Since this is core level forecast we use lowest level of everything.
core_mu <- dat$core.core.mu
core_sd <- dat$core.core.sd
plot_mu <- dat$plot.plot.mu
plot_sd <- dat$plot.plot.sd
site_mu <- dat$site.site.mu
site_sd <- dat$site.site.sd

#getr rid of relEM from plot_mu and plot_sd if present
plot_mu <- plot_mu[,!(colnames(plot_mu) %in% 'relEM')]
plot_sd <- plot_sd[,!(colnames(plot_sd) %in% 'relEM')]

#2. Draw covariates.----
#generic function.
sim_covs <- function(x,y){
  to.return <- list()
  mu <- as.matrix(x[,!(colnames(x) %in% indices)])
  sd <- as.matrix(y[,!(colnames(y) %in% indices)])
  if(ncol(mu) == 1){
    colnames(mu) <- colnames(x)[ncol(x)]
    colnames(sd) <- colnames(x)[ncol(x)]
  }
  for(i in 1:ncol(mu)){
    to.return[[i]] <- rnorm(nrow(x), mu[,i], sd[,i])
  }
  to.return <- do.call(cbind,to.return)
  colnames(to.return) <- colnames(mu)
  index <- data.frame(x[,(colnames(x) %in% indices)])
  to.return <- data.frame(cbind((index), to.return))
  colnames(to.return)[1:ncol(index)] <- colnames(x)[1:ncol(index)]
  return(to.return)
}
core_now <- sim_covs(core_mu,core_sd)
plot_now <- sim_covs(plot_mu,plot_sd)
site_now <- sim_covs(site_mu,site_sd)

#some one off transformations.
#add intercept at core level.
core_now <- cbind(rep(1,nrow(core_now)), core_now)
colnames(core_now)[1] <- 'intercept'

#map needs to be log-transformed at the site level.
if('map' %in% colnames(site_now)){site_now$map <- log(site_now$map)}

#b.relEM needs to be inverse logit transformed then multiplied by 100 (% EM trees, 0-100).
#then update name
if('b.relEM' %in% colnames(plot_now)){
  plot_now$b.relEM <- (boot::inv.logit(plot_now$b.relEM))*100
  colnames(plot_now)[colnames(plot_now)=="b.relEM"] <- "relEM"
  }

#3. Draw parameters, organize by core-plot-site.----
pars <- mcmc[sample(nrow(mcmc), 1, replace = T),]
#oraganize into a matrix.
n.spp = 4
pars <- matrix(pars, ncol = n.spp)
rownames(pars) <- mod$species_parameter_output$other$predictor



