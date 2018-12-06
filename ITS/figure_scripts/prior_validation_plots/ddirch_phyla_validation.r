#checking phyla level prior fits.
rm(list=ls())
source('paths.r')

#load phyla fit.
d <- readRDS(ted_ITS.prior_phyla_JAGSfit)

fitted <- d$all.preds$predicted
obs <- d$all.preds$observed

#ascomycetes
par(mfrow=c(1,2))
plot(obs[,2] ~ fitted[,2])
mod <- lm(obs[,2] ~ fitted[,2])
rsq <- summary(mod)$r.squared
rsq <- round(rsq,2)
abline(0,1,lwd = 2)
abline(mod, lty =2 , col = 'purple')
mtext(paste0('R2=',rsq),side = 3, line = -2, adj=0.05)
mtext(colnames(obs)[2], side = 3)
#basidiomycetes
plot(obs[,3] ~ fitted[,3])
mod <- lm(obs[,3] ~ fitted[,3])
rsq <- summary(mod)$r.squared
rsq <- round(rsq,2)
abline(0,1,lwd = 2)
abline(mod, lty =2 , col = 'purple')
mtext(paste0('R2=',rsq),side = 3, line = -2, adj=0.05)
mtext(colnames(obs)[3], side = 3)
