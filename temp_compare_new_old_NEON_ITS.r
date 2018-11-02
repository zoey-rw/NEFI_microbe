#Comparing ITS functional groups and genera with new and old aggregation strategies.
rm(list=ls())
source('paths.r')
source('NEFI_functions/crib_fun.r')
library(boot)

#Functional groups.----
#check different ways of data processing. old vs. new.
new <- readRDS(NEON_ITS_fastq_taxa_fg.path)
old <- readRDS(NEON_taxa_fg.path)
#for fucntional group abundances.
new <- new$rel.abundances
old <- old$rel.abundances

#for seq totals.
#new <- new$seq_total
#names(old$seq_total) <- old$abundances$geneticSampleID
#old <- old$seq_total

#match them up- functional groups.
old <- old[old$geneticSampleID %in% new$geneticSampleID,]
new <- new[new$geneticSampleID %in% old$geneticSampleID,]
old <- old[order(old$geneticSampleID),]
new <- new[order(new$geneticSampleID),]

#matching seq totals
old <- old[names(old) %in% names(new)]
new <- new[names(new) %in% names(old)]
old <- old[order(names(old))]
new <- new[order(names(new))]

a <- new$Saprotroph
b <- old$Saprotroph
c <- new$Ectomycorrhizal
d <- old$Ectomycorrhizal
a <- crib_fun(a)
b <- crib_fun(b)
c <- crib_fun(c)
d <- crib_fun(d)

mod.sap <- lm(logit(a) ~ logit(b))
mod.ecm <- lm(logit(c) ~ logit(d))

par(mfrow=c(1,2))
plot(logit(a) ~ logit(b))
abline(mod.sap, lwd = 2)
mtext(paste0('r2= ',round(summary(mod.sap)$r.squared,2)), side = 3, adj = 0.05, line = -2)
plot(logit(c) ~ logit(d))
abline(mod.ecm, lwd = 2)
mtext(paste0('r2= ',round(summary(mod.ecm)$r.squared,2)), side = 3, adj = 0.05, line = -2)


#Check taxa correlations.----