# check fits of 16S Bahram N-cycler  models.

rm(list=ls())
source('paths.r')
library(runjags)
library(betareg)
library(coda)
source('NEFI_functions/crib_fun.r')

#Load JAGS model.
allfits <- readRDS(bahram_16S_prior_N_cycle_JAGSfits)

#fit.mcmc <- as.mcmc(fit$jags_model) 

#check convergence plots and prsf scores.
#fit$species_parameter_output
#summary(fit$jags_model)
#plot(fit$jags_model)

#pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/prior_fit_cop_olig_16S.pdf")
pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/figures/prior_fit_ddirch_N_cycler.pdf")

#check the plots.
par(mfrow = c(1,3))
for (i in 1:length(allfits)){
fit <- allfits[[i]]
#fit <- fit$all.preds
#fit <- fit$climate.preds
fit <- fit$site.preds
#par(mar=c(3,3,3,3))
#for(i in 1:ncol(fit$predicted)){
i <- 2 #don't need the "other" model since it is the same data (abundance subtracted from 1)
  plot(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i], 
       pch = 16, xlab="predicted abundance", ylab="observed abundance")
  Axis(x="predicted", side=2)
  abline(0,1,lwd = 2)
  abline(lm(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i]), lty = 2, col = 'purple')
  mod <- betareg::betareg(crib_fun(fit$observed[,i]/rowSums(fit$observed)) ~ crib_fun(fit$predicted[,i]))
  rsq <-round(summary(mod)$pseudo.r.squared, 3)
  mtext(colnames(fit$predicted)[i], side = 3)
  mtext(paste0('R2=',rsq), side = 3, line = -1.5, adj = 0.05)
  
#}
}
dev.off()
