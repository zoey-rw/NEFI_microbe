# check fits of 16S Bahram functional group  models.

rm(list=ls())
source('paths.r')
library(runjags)
library(betareg)
library(coda)
library(ddpcr) # for quiet()
source('NEFI_functions/crib_fun.r')

#Load JAGS model.
# allfits <- readRDS(prior_16S_all.fg.groups_JAGSfits.path)
# allfits <- readRDS(bahram_16S_prior_phylo.group_JAGSfits)
allfits <- readRDS(paste0(scc_gen_16S_dir,"/JAGS_output/prior_phylo_JAGSfit_fewer_taxa.rds"))
phyla <- readRDS(paste0(scc_gen_16S_dir,"/JAGS_output/prior_phylo_JAGSfit_phylum_fewer_taxa_more_burnin.rds"))

allfits$phylum <- phyla$phylum

# allfits <- readRDS(paste0(scc_gen_16S_dir,"/JAGS_output/prior_phylo_JAGSfit_class.rds")) # did not converge with no nutrients
# allfits <- readRDS(paste0(scc_gen_16S_dir,"/JAGS_output/prior_phylo_JAGSfit_order.rds")) # order doesn't converge when nutrients are included
# allfits <- readRDS(paste0(scc_gen_16S_dir,"/JAGS_output/prior_phylo_JAGSfit_family.rds")) # didn't converge at all for either
# allfits <- readRDS(paste0(scc_gen_16S_dir,"/JAGS_output/prior_phylo_JAGSfit_genus.rds")) # genus mostly DID converge for all predictors, but not at all when nutrients were excluded

#check convergence plots and prsf scores.
#fit$species_parameter_output
#summary(fit$jags_model)
#plot(fit$jags_model)
#fit.mcmc <- as.mcmc(fit$jags_model) 

# has everything converged? print any high prsf scores.
for (i in 1:length(allfits)) {
  print(names(allfits)[i])
  fit <- allfits[[i]]
  fit <- fit$no.nutr.preds
  #print(colnames(fit$observed)[2])
  quiet(s <- summary(fit$jags_model))
  print("no nutr")
  print(s[which(s[,11] > 1.1),])
  
  fit <- allfits[[i]]
  fit <- fit$all.preds
  #print(colnames(fit$observed)[2])
  quiet(a <- summary(fit$jags_model))
  print("all preds")
  print(a[which(a[,11] > 1.1),])
}
fit$species_parameter_output
# save plots.
#pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/figures/prior_fit_ddirch_fg_no_nutr_fig.pdf")
pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/figures/prior_fit_ddirch_fg_fig.pdf")

#check the plots.
par(mfrow = c(2,2))
for (p in 1:length(allfits)) {
fit <- allfits[[p]]
#fit <- fit$all.preds
fit <- fit$no.nutr.preds

for(i in 2:ncol(fit$predicted)) {
  plot(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i], 
       pch = 16, xlab="predicted abundance", ylab="observed abundance")
  Axis(x="predicted", side=2)
  abline(0,1,lwd = 2)
  abline(lm(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i]), lty = 2, col = 'purple')
  mod <- betareg::betareg(crib_fun(fit$observed[,i]/rowSums(fit$observed)) ~ crib_fun(fit$predicted[,i]))
  rsq <-round(summary(mod)$pseudo.r.squared, 3)
  mtext(colnames(fit$predicted)[i], side = 3)
  mtext(paste0('R2=',rsq), side = 3, line = -1.5, adj = 0.05)
}
}
dev.off()

