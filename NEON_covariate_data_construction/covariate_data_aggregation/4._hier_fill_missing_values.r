# Filling in missing values at NEON cores, plots, and sites, hierarchically.
# once data have been aggregated at core-plot-site level, we need to fill in missing values and uncertainties hierarchically.
# this product is then shipped to the forecast, which draws missing data.
# NOTE: CPER, STER and WOOD are plains. You need to tell the filling that there are no ectomycorrhizal trees here.

rm(list=ls())
source('paths.r')
source('NEFI_functions/precision_matrix_match.r')

#set output path.
output.path <- hierarch_filled_data.path

#specify sites with zero ectomycorrhizal trees (plains sites)
no.ecm <- c('CPER','STER','WOOD')

#load prior model fit- model fit at site level.
# mod <- readRDS(bahram_16S_prior_dmulti.ddirch_all.group_JAGSfits) # can't read newer objects??
# mod <- mod$phylum
# preds <- mod$species_parameter_output$other$predictor

# adding predictors manually since they differ between fungi and bacteria.
preds <- c("NPP", "cn", "conifer", "forest", "intercept", "map", "mat", "pC", "pH", "pH_water", "relEM", "ndep.glob")
if('relEM' %in% preds){preds <- c(as.character(preds),'b.relEM')}

#add plot, site and sampleID to preds.
keys <- c('sampleID','plotID','siteID')
check <- c(as.character(preds),keys)
check_sd <-c(paste0(as.character(preds), '_sd'), keys)

#load NEON data.
core_obs  <- readRDS( core_obs_data.path)
core_core <- readRDS(core_core_data.path)
core_plot <- readRDS(core_plot_data.path)
core_site <- readRDS(core_site_data.path)
core_glob <- readRDS(core_glob_data.path)
plot_plot <- readRDS(plot_plot_data.path)
plot_site <- readRDS(plot_site_data.path)
plot_glob <- readRDS(plot_glob_data.path)
site_site <- readRDS(site_site_data.path)
site_glob <- readRDS(site_glob_data.path)

#cores, plots and sites that need to be present.
needed <- list()
needed[[1]] <- as.character(unique(core_obs$geneticSampleID))
needed[[2]] <- as.character(unique(core_obs$plotID))
needed[[3]] <- as.character(unique(core_obs$siteID))
names(needed) <- c('geneticSampleID','plotID','siteID')

#all y-cores present in core_core geneticSampleIDs?
sum(unique(core_obs$geneticSampleID) %in% core_core$geneticSampleID) == length(unique(core_obs$geneticSampleID))
nrow(core_obs) == nrow(core_core)
#all plots in core_core present in core_plot?
sum(unique(core_core$plotID) %in% core_plot$plotID) == length(unique(core_core$plotID))
#all sites in core_core present in core_site?
sum(unique(core_core$siteID) %in% core_site$siteID) == length(unique(core_core$siteID))
#all y-plots present in plot_plot plotIDs? 
sum(unique(core_obs$plotID) %in% plot_plot$plotID) == length(unique(core_obs$plotID))
#all sites in plot_plot present in plot_site?
sum(unique(plot_plot$siteID) %in% plot_site$siteID) == length(unique(plot_plot$siteID))
#all y-sites present in site_site siteIDs?
sum(unique(core_obs$siteID) %in% site_site$siteID) == length(unique(core_obs$siteID))



#Get core-plot-site indexing variables.
sites <- site_site$siteID
plots <- plot_plot$plotID
cores <- core_core$sampleID
plot.site <- substring(plots,1,4)
core.plot <- substring(cores,1,8)
core.site <- substring(cores,1,4)

#subset based on predictors actually in the model.
core_site_mu <- core_site[,colnames(core_site) %in% (check   )]
core_site_sd <- core_site[,colnames(core_site) %in% (check_sd)]
core_site_sd <- precision_matrix_match(core_site_mu,core_site_sd)
core_plot_mu <- core_plot[,colnames(core_plot) %in% (check   )]
core_plot_sd <- core_plot[,colnames(core_plot) %in% (check_sd)]
core_plot_sd <- precision_matrix_match(core_plot_mu,core_plot_sd)
core_core_mu <- core_core[,colnames(core_core) %in% (check   )]
core_core_sd <- core_core[,colnames(core_core) %in% (check_sd)]
core_core_sd <- precision_matrix_match(core_core_mu,core_core_sd)
plot_site_mu <- plot_site[,colnames(plot_site) %in% (check   )]
plot_site_sd <- plot_site[,colnames(plot_site) %in% (check_sd)]
plot_site_sd <- precision_matrix_match(plot_site_mu,plot_site_sd)
plot_plot_mu <- plot_plot[,colnames(plot_plot) %in% (check   )]
plot_plot_sd <- plot_plot[,colnames(plot_plot) %in% (check_sd)]
plot_plot_sd <- precision_matrix_match(plot_plot_mu,plot_plot_sd)
site_site_mu <- site_site[,colnames(site_site) %in% (check   )]
site_site_sd <- site_site[,colnames(site_site) %in% (check_sd)]
site_site_sd <- precision_matrix_match(site_site_mu,site_site_sd)

#Step 1. Missing data model. Assign missing means and sd hierarchically as necessary
#core_site
for(j in 1:ncol(core_site_mu)){
  name <- colnames(core_site_mu)[j]
  if(name %in% preds){
    core_site_mu[is.na(core_site_mu[,j]),j] <- core_glob[core_glob$pred == name,'Mean']
    core_site_sd[is.na(core_site_sd[,j]),j] <- core_glob[core_glob$pred == name,'SD'  ]
  }
}
#core_plot
for(j in 1:ncol(core_plot_mu)){
  name <- colnames(core_plot_mu)[j]
  if(name %in% preds){
    for(i in 1:nrow(core_plot_mu)){
      c.site <- core_plot_mu$siteID[i]
      if(is.na(core_plot_mu[i,j])){
        core_plot_mu[i,j] <- core_site_mu[core_site_mu$siteID == c.site,name]
        core_plot_sd[i,j] <- core_site_sd[core_site_sd$siteID == c.site,name]
      }
    }
  }
}
#core_core
for(j in 1:ncol(core_core_mu)){
  name <- colnames(core_core_mu)[j]
  if(name %in% preds){
    for(i in 1:nrow(core_core_mu)){
      c.plot <- core_core_mu$plotID[i]
      if(is.na(core_core_mu[i,j])){
        core_core_mu[i,j] <- core_plot_mu[core_plot_mu$plotID == c.plot,name]
        core_core_sd[i,j] <- core_plot_sd[core_plot_sd$plotID == c.plot,name]
      }
    }
  }
}
#plot_site
for(j in 1:ncol(plot_site_mu)){
  name <- colnames(plot_site_mu)[j]
  if(name %in% preds){
    plot_site_mu[is.na(plot_site_mu[,j]),j] <- plot_glob[plot_glob$pred == name,'Mean']
    plot_site_sd[is.na(plot_site_sd[,j]),j] <- plot_glob[plot_glob$pred == name,'SD'  ]
  }
  if(name == 'b.relEM'){
    plot_site_mu[plot_site_mu$siteID %in% no.ecm,j] <- -10
    plot_site_sd[plot_site_sd$siteID %in% no.ecm,j] <- 0.01
  }
}
#plot_plot
for(j in 1:ncol(plot_plot_mu)){
  name <- colnames(plot_plot_mu)[j]
  if(name %in% preds){
    for(i in 1:nrow(plot_plot_mu)){
      c.site <- plot_plot_mu$siteID[i]
      if(is.na(plot_plot_mu[i,j]) & name %in% colnames(plot_site_mu)){
        plot_plot_mu[i,j] <- plot_site_mu[plot_site_mu$siteID == c.site,name]
        plot_plot_sd[i,j] <- plot_site_sd[plot_site_sd$siteID == c.site,name]
      }
    }
  }
  if(name == 'b.relEM'){
    plot_plot_mu[plot_plot_mu$siteID %in% no.ecm, j] <- -10
    plot_plot_sd[plot_plot_sd$siteID %in% no.ecm, j] <- 0.01
  }
}


#site_site
site_site_mu <- site_site_mu[site_site_mu$siteID %in% core_obs$site,]
site_site_sd <- site_site_sd[site_site_sd$siteID %in% core_obs$site,]

# reorder core_obs and add deprecatedVialID to core_core_mu and core_core_sd
core_obs$geneticSampleID <- as.factor(core_obs$geneticSampleID)
core_obs <- core_obs[order(match(core_obs$geneticSampleID, core_core_mu$sampleID)),]
core_core_mu$deprecatedVialID <- core_obs$deprecatedVialID
core_core_sd$deprecatedVialID <- core_obs$deprecatedVialID

output.list <- list(core_obs,
                    core_core_mu,core_core_sd,
                    core_plot_mu,core_plot_sd,
                    core_site_mu,core_site_sd,
                    plot_plot_mu,plot_plot_sd,
                    plot_site_mu,plot_site_sd,
                    site_site_mu,site_site_sd)
names(output.list) <- c('core.obs',
                        'core.core.mu','core.core.sd',
                        'core.plot.mu','core.plot.sd',
                        'core.site.mu','core.site.sd',
                        'plot.plot.mu','plot.plot.sd',
                        'plot.site.mu','plot.site.sd',
                        'site.site.mu','site.site.sd')
#save output.
saveRDS(output.list, output.path)
