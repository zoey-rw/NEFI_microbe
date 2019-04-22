# remove missing data from NEON
# alternative to filling missing data hierarchically

source('paths.r')
source('NEFI_functions/precision_matrix_match.r')

#set output path.
output.path <- missing_data_removed_16S.path

#specify sites with zero ectomycorrhizal trees (plains sites)
no.ecm <- c('CPER','STER','WOOD')

#load prior model fit- model fit at site level.
mod <- readRDS(paste0(scc_gen_16S_dir, "/JAGS_output/prior_phylo_JAGSfit_phylum_fewer_taxa_more_burnin.rds"))
mod <- mod$phylum
preds <- mod$all.preds$species_parameter_output$other$predictor
if('relEM' %in% preds){preds <- c(as.character(preds),'b.relEM')}

#add plot, site and sampleID to preds.
keys <- c('sampleID','plotID','siteID')
check <- c(as.character(preds),keys)
check_sd <-c(paste0(as.character(preds), '_sd'), keys)

#load NEON data.
core_obs  <- readRDS( core_obs_16S.path)
core_core <- readRDS(core_core_16S.path)
core_plot <- readRDS(core_plot_16S.path)
core_site <- readRDS(core_site_16S.path)
core_glob <- readRDS(core_glob_16S.path)
plot_plot <- readRDS(plot_plot_16S.path)
plot_site <- readRDS(plot_site_16S.path)
plot_glob <- readRDS(plot_glob_16S.path)
site_site <- readRDS(site_site_16S.path)
site_glob <- readRDS(site_glob_16S.path)

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
#all y-plots present in plot_plot plotIDs? FALSE
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

# fill in no-ecto data for 3 sites.

# plot_site
for(j in 1:ncol(plot_site_mu)){
  name <- colnames(plot_site_mu)[j]
  if(name == 'b.relEM'){
    plot_site_mu[plot_site_mu$siteID %in% no.ecm,j] <- -10
    plot_site_sd[plot_site_sd$siteID %in% no.ecm,j] <- 0.01
  }
}
#plot_plot
for(j in 1:ncol(plot_plot_mu)){
  name <- colnames(plot_plot_mu)[j]
  if(name == 'b.relEM'){
    plot_plot_mu[plot_plot_mu$siteID %in% no.ecm, j] <- -10
    plot_plot_sd[plot_plot_sd$siteID %in% no.ecm, j] <- 0.01
  }
}

# remove any rows with missing data
df.list <- list(core_core_mu,core_core_sd,
                core_plot_mu,core_plot_sd,
                core_site_mu,core_site_sd,
                plot_plot_mu,plot_plot_sd,
                plot_site_mu,plot_site_sd,
                site_site_mu,site_site_sd)
#output.list <- lapply(df.list, function(x) x[complete.cases(x),])

core_obs$sampleID <- substr(core_obs$geneticSampleID, 1, nchar(core_obs$geneticSampleID)-4) 
core_obs <- core_obs[core_obs$sampleID %in% core_core_mu$sampleID,]

output.list <- append(list(core_obs),df.list)
names(output.list) <- c('core.obs',
                        'core.core.mu','core.core.sd',
                        'core.plot.mu','core.plot.sd',
                        'core.site.mu','core.site.sd',
                        'plot.plot.mu','plot.plot.sd',
                        'plot.site.mu','plot.site.sd',
                        'site.site.mu','site.site.sd')
saveRDS(output.list, output.path)

