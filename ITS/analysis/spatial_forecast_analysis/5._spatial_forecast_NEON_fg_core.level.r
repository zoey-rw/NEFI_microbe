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
source('NEFI_functions/sd_to_precision.r')
source('NEFI_functions/precision_matrix_match.r')
library(runjags)

#set number of simulations and output path.
n.sim <- 1000

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

#get rid of relEM from plot_mu and plot_sd if present
plot_mu <- plot_mu[,!(colnames(plot_mu) %in% 'relEM')]
plot_sd <- plot_sd[,!(colnames(plot_sd) %in% 'relEM')]

#get indexing for downstream plot and site aggregation of predictions.
core_plot <- droplevels(core_mu$plotID)
core_site <- droplevels(core_mu$siteID)
plot_site <- droplevels(unique(core_mu$plotID))
plot_site <- substring(plot_site, 1,4)

#load observed fungal data for validation.
val <- readRDS(NEON_taxa_fg.path)
val <- val$rel.counts
val$sampleID <- gsub('-GEN','',val$geneticSampleID)
nrow(core_mu[core_mu$sampleID %in% val$sampleID,])


#2. Begin simulation loop!----
#objects to store predictions.
pred.x.m.out <- list()
cred.out <- list()
pred.out <- list()

for(i in 1:n.sim){
  #2a. Draw covariates.----
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
  
  #2b. Draw parameters, organize by core-plot-site.----
  pars <- mcmc[sample(nrow(mcmc), 1, replace = T),]
  #oraganize into a matrix.
  pars <- matrix(pars, ncol = n.spp)
  rownames(pars) <- mod$species_parameter_output$other$predictor
  
  #merge together core-plot-site predictors into a single data frame.
  preds <- merge(core_now, plot_now)
  preds <- merge(preds, site_now)
  #put preds in same order as core_mu for downstream matching.
  preds <- preds[order(match(preds$sampleID, core_mu$sampleID)),]
  
  #2c. combine parameters and covariates to get predicted values for all species.----
  pred.x.m <- list()
  for(k in 1:ncol(pars)){
    pred.spp <- preds[colnames(preds) %in% rownames(pars)]
    pred.spp <- pred.spp[,order(match(colnames(pred.spp), rownames(pars)))]
    pred.x.m[[k]] <- exp(as.matrix(pred.spp) %*% pars[,k])
  }
  pred.x.m <- do.call(cbind,pred.x.m)
  
  #2d. save credible and predictive interval on original scale.----
  pred.x.m.out[[i]] <- pred.x.m
  cred.out[[i]] <- pred.x.m / rowSums(pred.x.m)
  pred.out[[i]] <- DirichletReg::rdirichlet(nrow(pred.x.m),pred.x.m)
  
  #2e. close loop.----
}

#3. aggregate credible and predictive interval at core, plot and site scale.----
#3a. calculate credible and predictive intervals core scale.----
pred.mean       <- apply(simplify2array(cred.out), 1:2, mean)
pred.ci.0.025   <- apply(simplify2array(cred.out), 1:2, quantile, probs = c(0.025))
pred.ci.0.975   <- apply(simplify2array(cred.out), 1:2, quantile, probs = c(0.975))
pred.pi.0.025   <- apply(simplify2array(pred.out), 1:2, quantile, probs = c(0.025))
pred.pi.0.975   <- apply(simplify2array(pred.out), 1:2, quantile, probs = c(0.975))
pred.sd <- apply(array(unlist(cred.out), c(nrow(cred.out[[1]]), ncol(cred.out[[1]]), length(cred.out))), c(1,2), sd)
pred.precision <- sd_to_precision(pred.sd)
#bind together in list.
core.pred <- list(pred.mean,pred.ci.0.025,pred.ci.0.975,pred.pi.0.025,pred.pi.0.975)
for(i in 1:length(core.pred)){
  rownames(core.pred[[i]]) <- core_mu$sampleID
  colnames(core.pred[[i]]) <- names(mod$species_parameter_output)
}

#3b. validate core-cast.----
#validation <- core.pred[[1]]
#validation <- validation[rownames(validation) %in% val$sampleID,]
#validation <- validation[order(rownames(validation)),]
#validation <- data.frame(validation)
#val <- val[order(val$sampleID),]
#plot it.
#plot(val$Ectomycorrhizal ~ validation$Ectomycorrhizal)


#3b. Specify JAGS aggregation model.----
jags.model = "
model {
  #draw y observations from their distributions.
  #for(i in 1:ncol(y_mu)){
  #  y[,i] ~ dnorm(y_mu[,i], y_precision[,i])
  #}

  # Observations (single set per core):
  for(i in 1:N.core){
    y[i,1:N.spp] ~ ddirch(plot_mu[core_plot[i],1:N.spp] * core_alpha) 
  }
  # Plot means:
  for(i in 1:N.plot){
    plot_mu[i,1:N.spp] ~ ddirch(site_mu[plot_site[i],1:N.spp] * plot_alpha)
  }
  
  # Site means:
  for(i in 1:N.site){
    site_mu[i,1:N.spp] ~ ddirch(priors[1:N.spp] * site_alpha)
  }
  
  # Priors:
  for(s in 1:N.spp){
    priors[s] <- 1
  }
  site_alpha ~ dgamma(0.01, 0.01)
  plot_alpha ~ dgamma(0.01, 0.01)
  core_alpha ~ dgamma(0.01, 0.01)
  
}" #end jags model.

#3c. setup JAGS data objects.----
jd <- list(y=pred.mean, core_plot = as.factor(core_plot), plot_site = as.factor(plot_site),
           N.spp = ncol(pred.mean), N.core = nrow(pred.mean), N.plot = length(plot_site), N.site = length(unique(plot_site)))

#3d. run JAGS model in runjags.----
jags.out <- run.jags(model = jags.model,
                    data = jd,
                    n.chains = 3,
                    adapt = 500,
                    burnin = 1000,
                    sample = 1000,
                    monitor = c('plot_mu','site_mu'))

#3e. check site-aggreagted core-cast vs. old site-cast.----
#WHAT ORDER ARE THESE SITES IN???
site.out <- summary(jags.out, var = 'site_mu')
site.mean <- matrix(site.out[,4], ncol = ncol(pred.mean))
colnames(site.mean) <- names(mod$species_parameter_output)
rownames(site.mean) <- unique(plot_site)

#check against our original site-cast.
site_cast <- readRDS(NEON_site_fcast_fg.path)
site_cast <- site_cast$all.preds$mean

check <- site.mean[rownames(site.mean) %in% rownames(site_cast),]

