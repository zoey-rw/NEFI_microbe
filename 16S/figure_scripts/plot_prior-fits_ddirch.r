#check fits, convergence plots and prsf scores for calibration models

rm(list=ls())
source('paths.r')
library(runjags)
library(betareg)
library(coda)
library(ddpcr) # for quiet() function
source('NEFI_functions/crib_fun.r')

allfits <- readRDS(prior_delgado_ddirch_16S.path)

# has everything converged? print any high prsf scores.
for (i in 1:length(allfits)) {
  print(names(allfits)[[i]])
  fit <- allfits[[i]]
  summary(fit$jags_model)
  quiet(s <- summary(fit$jags_model))
  print(s[which(s[,11] > 1.1),])
}
fit$species_parameter_output

# save plots.
pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/figures/prior_fit_ddirch_all_groups.pdf")

#check the plots.
allfits <- output.list
r2_all <- list()
par(mfrow = c(3,3))

for (p in 1:length(allfits)) {
fit <- allfits[[p]]
r2_lev <- list()
for(i in 1:ncol(fit$predicted)) {
  if (colnames(fit$predicted)[i]=="other") next()
  plot(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i],
       pch = 16, xlab="predicted abundance", ylab="observed abundance")
  Axis(x="predicted", side=2)
  abline(0,1,lwd = 2)
  abline(lm(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i]), lty = 2, col = 'purple')
  mod <- betareg::betareg(crib_fun(fit$observed[,i]/rowSums(fit$observed)) ~ crib_fun(fit$predicted[,i]))
  rsq <-round(summary(mod)$pseudo.r.squared, 3)
  mtext(colnames(fit$predicted)[i], side = 3)
  mtext(paste0('R2=',rsq), side = 3, line = -1.5, adj = 0.05)
  r2_lev[[i]] <- rsq 
}
r2_all[[p]] <- r2_lev
}
dev.off()

# view summary by rank
lev.mu <- lapply(r2_all[1:5], unlist)
lev.mu <- lapply(lev.mu, mean)
tax.mu <- mean(unlist(r2_all[1:5]))
fg.mu <- mean(unlist(r2_all[6:18]))
n.mu <- mean(unlist(r2_all[7:13]))
c.mu <- mean(unlist(r2_all[c(6,14:16)]))
co.mu <- mean(unlist(r2_all[c(17:18)]))
