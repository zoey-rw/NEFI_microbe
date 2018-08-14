#summarizing model fits.
rm(list=ls())
library(boot)
source('paths.r')

fun <- readRDS(ted_ITS.prior_fg_JAGSfit)
gen <- readRDS(ted_ITS.prior_20gen_JAGSfit)

#get R2 values
 all <- fun$all.preds
site <- fun$site.preds
clim <- fun$climate.preds

r2.list <- list()
for(i in 1:ncol(all$observed)){
  name.spp <- names(all$species_parameter_output)[i]
  mod1 <- lm(logit( all$observed[,i]) ~ logit( all$predicted[,i]))
  mod2 <- lm(logit(site$observed[,i]) ~ logit(site$predicted[,i]))
  mod3 <- lm(logit(clim$observed[,i]) ~ logit(clim$predicted[,i]))
  rsq1 <- summary(mod1)$r.squared
  rsq2 <- summary(mod2)$r.squared
  rsq3 <- summary(mod3)$r.squared
  rsq <- c(rsq1,rsq2,rsq3)
  rsq <- round(rsq, 3)
  r2.list[[i]] <- c(name.spp,rsq)
}
out <- data.frame(do.call(rbind,r2.list))
colnames(out) <- c('group','all','site','climate')
as.numeric(out[,2:4])

#get genera level R2 values
all.gen <- gen$all.preds
r2.list <- list()
for(i in 1:ncol(all.gen$observed)){
  name.spp <- names(all.gen$species_parameter_output)[i]
  mod <- lm(logit(all.gen$observed[,i]) ~ logit(all.gen$predicted[,i]))
  rsq <- summary(mod)$r.squared
  r2.list[[i]] <- c(name.spp,rsq)
}
r2 <- data.frame(do.call(rbind,r2.list))
r2[,2] <- as.numeric(as.character(r2[,2]))
hist(r2[,2])

#fg fits.
par(mfrow = c(2,2))
for(i in 2:5){
  name.spp <- names(all$species_parameter_output)[i]
  mod1 <- lm(inv.logit( all$observed[,i]) ~ inv.logit( all$predicted[,i]))
  mod2 <- lm(inv.logit(site$observed[,i]) ~ inv.logit(site$predicted[,i]))
  mod3 <- lm(inv.logit(clim$observed[,i]) ~ inv.logit(clim$predicted[,i]))
  rsq1 <- summary(mod)$r.squared
  plot(inv.logit(all$observed[,i]) ~ inv.logit(all$predicted[,i]), cex = 0.7, main = name.spp)
  abline(mod, lty =2, lwd = 2, col = 'purple')
  mtext(paste0('R2 = ',round(rsq, 2)), side = 3, line = -2, adj = 0.05)
}


#20 genera fits.
#best fit is inocybe, -4.
plot(logit(all.gen$observed[,4]) ~ logit(all.gen$predicted[,4]))
#next best fit is tomentella - 5.
plot(logit(all.gen$observed[,4]) ~ logit(all.gen$predicted[,4]))
hist(r2[,2])
