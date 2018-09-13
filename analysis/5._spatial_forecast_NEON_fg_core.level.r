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

#load prior model fit- model fit at site level.
mod <- readRDS(ted_ITS.prior_fg_JAGSfit)

#load NEON data.
core_core <- readRDS(core_core.path)
core_plot <- readRDS(core_plot.path)
core_site <- readRDS(core_site.path)
core_glob <- readRDS(core_glob.path)
plot_plot <- readRDS(plot_plot.path)
plot_site <- readRDS(plot_site.path)
plot_glob <- readRDS(plot_glob.path)
site_site <- readRDS(site_site.path)
site_glob <- readRDS(site_glob.path)


#tell function with preds you are using and indexing
sites <- site_site$siteID
plots <- plot_plot$plotID
cores <- core_core$sampleID
plot.site <- substring(plots,1,4)
core.plot <- substring(cores,1,8)
core.site <- substring(cores,1,4)

#Some plots in core-level not in plot-level.
to_add <- as.character(unique(core_core[!(core_core$plotID %in% core_plot$plotID),]$plotID))
to_add <- data.frame(to_add)
colnames(to_add) <- 'plotID'
test <- plyr::rbind.fill(core_plot,to_add)


preds <- mod$all.preds$species_parameter_output$other$predictor
if('relEM' %in% preds){preds <- c(as.character(preds),'b.relEM')}
#add plot, site and sampleID to preds.
keys <- c('sampleID','plotID','siteID')
check <- c(as.character(preds),keys)
check_sd <-c(paste0(as.character(preds), '_sd'), keys)
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

#drop sites that aren't present at core-level
test <- site_site_mu[site_site_mu$siteID %in% core_core_mu$siteID,]


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
}
#plot_plot
for(j in 1:ncol(plot_plot_mu)){
  name <- colnames(plot_plot_mu)[j]
  if(name %in% preds){
    for(i in 1:nrow(plot_plot_mu)){
      c.site <- plot_plot_mu$siteID[i]
      if(is.na(plot_plot_mu[i,j])){
        plot_plot_mu[i,j] <- plot_site_mu[plot_site_mu$siteID == c.site,name]
        plot_plot_sd[i,j] <- plot_site_sd[plot_site_sd$siteID == c.site,name]
      }
    }
  }
}



