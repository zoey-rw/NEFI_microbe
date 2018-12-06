rm(list=ls())
source('paths.r')
library(runjags)
library(betareg)
source('NEFI_functions/crib_fun.r')

#Check cosmo genera.
#fit <- readRDS(bahram_16S.prior_12gen_JAGSfit)
#fit <- fit$all.preds

# Check top 15 cosmo phyla.
#fit <- readRDS(bahram_16S.prior_15phyla_JAGSfit)
#fit <- fit$all.preds

# Check by phylogenetic group.
fit <- readRDS(bahram_16S_prior_phylo.group_JAGSfits)
fit <- fit$phylum

#check convergence plots and prsf scores.
fit$species_parameter_output
summary(fit$jags_model)
#plot(fit$jags_model)

# save plots.
#pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/prior_fit_cosmo_16S.pdf")
#pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/prior_fit_phyla_16S.pdf")
#pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/prior_fit_phylo.group_phyla_no_mois_16S.pdf")

#check the plots.
par(mfrow = c(3,3))
#par(mar=c(3,3,3,3))
for(i in 1:ncol(fit$predicted)){
  plot(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i], pch = 16)
  Axis(x="predicted", side=2)
  abline(0,1,lwd = 2)
  abline(lm(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i]), lty = 2, col = 'purple')
  mod <- betareg::betareg(crib_fun(fit$observed[,i]/rowSums(fit$observed)) ~ crib_fun(fit$predicted[,i]))
  rsq <-round(summary(mod)$pseudo.r.squared, 2)
  mtext(colnames(fit$predicted)[i], side = 3)
  mtext(paste0('R2=',rsq), side = 3, line = -1.5, adj = 0.05)
  
}
dev.off()
