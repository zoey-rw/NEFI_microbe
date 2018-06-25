#Fit a hierarchical core/plot/site jags model.
#clearn envirionment, source paths.
rm(list=ls())
source('/home/caverill/NEFI_microbe/data_formatting/formatting_NEON_microbial/core_site_plot_aggregation_May.2018/0. aggregation paths.r')
#source some other functions.
source('/home/caverill/NEFI_microbe/NEFI_functions/crib_fun.r')
source('/home/caverill/NEFI_microbe/NEFI_functions/z_transform.r')
source('/home/caverill/NEFI_microbe/NEFI_functions/precision_matrix_match.r')
source('/home/caverill/NEFI_microbe/NEFI_functions/z_transform_table.r')
source('/home/caverill/NEFI_microbe/NEFI_functions/sd_to_precision.r')

#load data.
 obs.table <- readRDS( obs.table.path)
core.table <- readRDS(core.table.path)
plot.table <- readRDS(plot.table.path)
site.table <- readRDS(site.table.path)

##TESTING
#lets only look at complete cases for 2 sites right now.
 obs.table <-  obs.table[obs.table$siteID %in% c('BART','HARV'),]
core.table <- core.table[core.table$geneticSampleID %in% obs.table$geneticSampleID,]
plot.table <- plot.table[plot.table$plotID %in% obs.table$plotID,]
site.table <- site.table[site.table$siteID %in% obs.table$siteID,]

#put everything in correct order. 
#Convert all factor vectors to character vectors - this line sucks and is a good argument for data.table/tidyverse
 obs.table[sapply( obs.table, is.factor)] <-  lapply( obs.table[sapply( obs.table, is.factor)], as.character)
core.table[sapply(core.table, is.factor)] <-  lapply(core.table[sapply(core.table, is.factor)], as.character)
plot.table[sapply(plot.table, is.factor)] <-  lapply(plot.table[sapply(plot.table, is.factor)], as.character)
site.table[sapply(site.table, is.factor)] <-  lapply(site.table[sapply(site.table, is.factor)], as.character)

#hierarchically order the rows as site/plot/core.
 obs.table <-  obs.table[order( obs.table$siteID, obs.table$plotID, obs.table$geneticSampleID),]
core.table <- core.table[order(core.table$siteID,core.table$plotID,core.table$geneticSampleID),]
plot.table <- plot.table[order(plot.table$siteID,plot.table$plotID),]
site.table <- site.table[order(site.table$siteID),]

#setup indexing
core_plot <- core.table$plotID
core_site <- core.table$siteID
plot_plot <- plot.table$plotID
site_site <- site.table$siteID

#setup y, core, plot and site data frames, and associated sd tables where applicable.
y <- obs.table[,2:5]
#dropping soil moisture and soil pH for now because all NA values for two site-dates of interest.
#core <- core.table[,c('siteID','plotID','soilMoisture','soilInWaterpH','soilTemp')]
#colnames(core)[3:5] <- c('moisture','pH','temperature')
core <- core.table[,c('siteID','plotID','soilTemp')]
colnames(core)[3] <- c('temperature')
core_sd <- core.table[,c('siteID','plotID')]
plot <- plot.table[,c('siteID','plotID','organicCPercent','CNratio')]
plot_sd <- plot.table[,c('siteID','plotID','organicCPercent_sd','CNratio_sd')]
site <- site.table[,c('siteID','elevation','map','mat','n.dep')]
site_sd <- site.table[,c('siteID','map_sd','mat_sd')]


#transform y values.
y <- data.frame(lapply(y, crib_fun))

#fill out sd tables.
core_sd <- precision_matrix_match(core,core_sd)
plot_sd <- precision_matrix_match(plot,plot_sd)
site_sd <- precision_matrix_match(site,site_sd)

#If you are z-transforming then use these functions to transform mu and sd values.
z_trans = F
if(z_trans == T){
  z.core <- z_transform_table(core,core_sd)
  z.plot <- z_transform_table(plot,plot_sd)
  z.site <- z_transform_table(site,site_sd)
  core <- z.core[[1]]
  plot <- z.plot[[1]]
  site <- z.site[[1]]
  core_sd <- z.core[[2]]
  plot_sd <- z.plot[[2]]
  site_sd <- z.site[[2]]
}

#convert sd tables to precision tables.
core_precision <- sd_to_precision(core_sd)
plot_precision <- sd_to_precision(plot_sd)
site_precision <- sd_to_precision(site_sd)

#subset to numeric columns.
core <- Filter(is.numeric, core)
plot <- Filter(is.numeric, plot)
site <- Filter(is.numeric, site)
core_precision <- Filter(is.numeric, core_precision)
plot_precision <- Filter(is.numeric, plot_precision)
site_precision <- Filter(is.numeric, site_precision)

##TESTING
#just have one predictor at each level.
#this loses the column names. okay for now. 
core <- data.frame(core[,1])
plot <- data.frame(plot[,1])
site <- data.frame(site[,1])
core_precision <- data.frame(core_precision[,1])
plot_precision <- data.frame(plot_precision[,1])
site_precision <- data.frame(site_precision[,1])

#specify JAGS model. This is general for any linear combination of x's.
jags.model = "
model {
  ###priors for each species*predictor*hierarchical-level
  #core priors
  for(j in 1:N.spp){
    for(k in 1:N.core.preds){
      core.m[k,j] ~ dnorm(0, 1.0E-4)
    }
  }
  #plot priors
    for(j in 1:N.spp){
    for(k in 1:N.plot.preds){
      plot.m[k,j] ~ dnorm(0, 1.0E-4)
    }
  }
  #site priors
    for(j in 1:N.spp){
    for(k in 1:N.site.preds){
      site.m[k,j] ~ dnorm(0, 1.0E-4)
    }
  }

  ###predictor (x) values drawn from distributions.
  #core values
  for(j in 1:N.core.preds){
    for(i in 1:N.core){
      core[i,j] ~ dnorm(core_mu[i,j], core_precision[i,j])
    }
  }
  #plot values
  for(j in 1:N.plot.preds){
    for(i in 1:N.plot){
      plot[i,j] ~ dnorm(plot_mu[i,j], plot_precision[i,j])
    }
  }
  #site values
  for(j in 1:N.site.preds){
    for(i in 1:N.site){
      site[i,j] ~ dnorm(site_mu[i,j], site_precision[i,j])
    }
  }

  #combine predictors and parameters.
  #core level
  for(i in 1:N.core){
    for(j in 1:N.spp){
      core.level[i,j] <- inprod(m.core[,j], core[i,])
    }
  }
  #plot level
  for(i in 1:N.plot){
    for(j in 1:N.spp){
      plot.level[i,j] <- inprod(m.plot[,j], plot[plot_plot[i],])
    }
  }
  #site level
  for(i in 1:N.site){
    for(j in 1:N.spp){
      site.level[i,j] <- inprod(m.site[,j], site[site_site[i],])
    }
  }
  
  #go ahead and fit means of species abundances as a linear combination of predictors and parameters.
  ###I THINK THIS LOOP NEEDS MORE HIERARCHY.
  for(i in 1:N.core){
    for(j in 1:N.spp){
      log(a0[i,j]) <- core.level[i,j] + plot.level[core_plot[i],j] + site.level[core_site[i],j]
      }
    y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp]) 
  }

} #close model loop.
"

jags.data <- list(N.core = nrow(core), N.plot = nrow(plot), N.site = nrow(site), N.spp = ncol(y),
                  N.core.preds = ncol(core), N.plot.preds = ncol(plot), N.site.preds = ncol(site),
                  core_mu = core, plot_mu = plot, site_mu = site,
                  core_precision = core_precision, plot_precision = plot_precision, site_precision = site_precision,
                  y = y)