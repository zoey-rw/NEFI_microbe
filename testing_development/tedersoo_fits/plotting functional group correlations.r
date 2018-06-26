#plot function group correlations
rm(list=ls())
library(data.table)
library(betareg)
library(PerformanceAnalytics)
source('NEFI_functions/crib_fun.r')
d <- data.table(readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/ted_all_prior_data.rds'))
of_interest <- d[,.(Ectomycorrhizal,Arbuscular,Saprotroph,hydrophillic,hydrophobic)]
#setup data
of_interest <- lapply(of_interest,crib_fun)
of_interest <- lapply(of_interest,logit)
of_interest <- as.data.frame(of_interest)

png(filename='figures/ted_fungal_groups_correlations.png',width=10,height=10,units='in',res=300)
chart.Correlation(of_interest,  pch = 16, cex.labels = 2)
dev.off()

#Just plot of ECM vs. SAP.
png(filename='figures/ted_ECM.SAP_corr.png',width=5,height=5,units='in',res=300)
of_interest <- as.data.frame(of_interest)
par(mfrow= c(1,1),
    mai = c(0.1,0.1,0.1,0.1),
    oma = c(4,4,0,0))
plot(Saprotroph ~ Ectomycorrhizal, data = of_interest, cex = 0.7, pch = 16, ylab = NA, xlab = NA)
abline(lm(Saprotroph ~ Ectomycorrhizal, data = of_interest), lwd = 2)
r.sq <- summary(lm(Saprotroph ~ Ectomycorrhizal, data = of_interest))$r.sq
mtext(paste('r2=',round(r.sq,2)), side=3, line=-1.5, cex = 1.2, adj = 0.9)
mtext('Saprotroph Abundance'     , side = 2, line = 2.5, cex = 1.2)
mtext('Ectomycorrhizal Abundance', side = 1, line = 2.5, cex = 1.2)
dev.off()