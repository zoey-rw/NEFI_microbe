# Stacked barplots for bacteria @ prior observations, NEON forecasts, and NEON observations
rm(list=ls())
source('paths.r')
source('NEFI_functions/crib_fun.r')
library(data.table)

# read in prior abundances
prior.all <- readRDS(bahram_16S_common_phylo_fg_abun.path)

# read in NEON forecasts and observations
pl.cast <- readRDS(NEON_cps_fcast_ddirch_16S.path)
pl.truth <- readRDS(NEON_phylo_fg_plot.site_obs_16S.path)

p <- 2
for (p in c(6)){
  preds <- pl.cast[[p]]  
  prior <- prior.all[[p]]$rel.abundances
  raw.truth <- pl.truth[[p]]
  
  # only site-level
  fcast <- preds$site.fit
  truth <- raw.truth$site.fit
  
  # convert to data frame and add siteIDs
  pred <- as.data.frame(fcast$mean)
  pred$siteID <- rownames(pred)
  obs <- as.data.frame(truth$mean)
  obs$siteID <- rownames(obs)
  prior <- as.data.frame(prior)
  prior$siteID <- rownames(prior)
  
  obs <- obs[,match(colnames(pred), colnames(obs))]
  
  # melt to long format
  df_obs <- melt(obs, id.vars = c("siteID"))
  df_pred <- melt(pred, id.vars = c("siteID"))
  df_prior <- melt(prior, id.vars = c("siteID"))
  
  o <- ggplot(df_obs, aes(x = siteID, y = value, fill = variable)) + 
    geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 320, hjust=-.1)) + 
    ggtitle("NEON observed site-level averages of bacterial groups") + 
    labs(aesthetic=names(pl.cast)[[p]])  
  
  f <- ggplot(df_pred, aes(x = siteID, y = value, fill = variable)) + 
    geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 320, hjust=-.1))+ #legend.position = "none") + 
    ggtitle("NEON predicted site-level averages of bacterial groups")# + facet_grid(siteID ~ .)
  
  pri <- ggplot(df_prior, aes(x = siteID, y = value, fill = variable)) + 
    geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 320, hjust=-.1))+
                                        #, legend.position = "none") + 
    ggtitle("Prior observed averages of bacterial groups")# + facet_grid(siteID ~ .)
  
  
  grid.arrange(pri, f, o,
               layout_matrix = rbind(c(1, 1),
                                     c(2, 3))
  )
}