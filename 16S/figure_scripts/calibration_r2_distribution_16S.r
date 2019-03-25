#Plot in-sample calibration R2 values as a density plot.
rm(list=ls())
source('paths.r')

#Set output path.----
output.path <- prior_16S_r2_distribution.density_figure.path

#load data.----
#fg_all <- readRDS(prior_16S_all.fg.groups_JAGSfits.path)
#pl <- readRDS("/fs/data3/caverill/NEFI_data/16S/scc_gen/JAGS_output/bahram_16S.prior_phylo_new_test.rds")
pl <- readRDS(paste0(scc_gen_16S_dir,"/JAGS_output/prior_phylo_JAGSfit_fewer_taxa.rds"))
phyla <- readRDS(paste0(scc_gen_16S_dir,"/JAGS_output/prior_phylo_JAGSfit_phylum_fewer_taxa_more_burnin.rds"))

pl$phylum <- phyla$phylum

#loop over lists and get R2 values.----
#functional groups.
fg.r2.all <- list()
for (p in 1:length(fg_all)){
  fg <- fg_all[[p]]
# obs <- fg$all.preds$observed
# pred <- fg$all.preds$predicted
  obs <- fg$no.nutr.preds$observed
  pred <- fg$no.nutr.preds$predicted
fg.r2 <- list()
for(i in 1:ncol(obs)){
  fg.r2[[i]] <- summary(lm(obs[,i] ~ pred[,i]))$r.squared
  }
fg.r2 <- unlist(fg.r2)
names(fg.r2) <- colnames(obs)
#drop 'other'.
fg.r2 <- fg.r2[names(fg.r2) != 'other']
fg.r2.all[[p]] <- fg.r2
}
fg.r2.all <- unlist(fg.r2.all)
#phylogenetic groups.
pl.r2 <- list()
for(i in 1:length(pl)){
  lev <- pl[[i]]
  lev <- lev$all.preds
  obs <- lev$observed
  pred <- lev$predicted
  lev.r2 <- list()
  for(j in 1:ncol(obs)){lev.r2[[j]] <- summary(lm(obs[,j] ~ pred[,j]))$r.squared}
  lev.r2 <- unlist(lev.r2)
  names(lev.r2) <- colnames(obs)
  lev.r2 <- lev.r2[names(lev.r2) != 'other']
  pl.r2[[i]] <- lev.r2
  names(pl.r2)[[i]] <- names(pl)[i]
}
pl.r2 <- unlist(pl.r2)

#Merge all r2 values.
all.r2 <- c(fg.r2.all, pl.r2)

#Get truncated density so that we don't have densities less than zero.----
h <- density(all.r2)$bw  #get dansity bandwith
# Compute edge weights.
w <- 1 / pnorm(0, mean=all.r2, sd=h, lower.tail=FALSE)
#Generate truncated distribution.
h <- density(all.r2, bw=h, kernel="gaussian", weights=w / length(all.r2))
h$y[h$all.r2 < 0] <- 0
#Check: the integral ought to be close to 1:
if(sum(h$y * diff(h$x)[1]) > 1.1 | sum(h$y * diff(h$x)[1]) < 0.9){
  warning('Zero-truncated density sum is not close to 1. Consider doing something else...')
}


#Make density plot.----
#png save settings.
#png(filename=output.path,width=6,height=5,units='in',res=300)

#global plot settings.
trans <- 0.3 #shading transparency.
o.cex <- 1.3 #outer label size.
par(mfrow = c(1,1), mar = c(4.2,4.2,1.5,1.5))

#plot.
plot(h,xlim = c(0, 1), ylim = c(0, round(max(h$y), 0)), bty = 'n', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1)
polygon(h, col = adjustcolor('purple',trans), border = NA)
mtext('Density', side = 2, line = 2.2, cex = o.cex)
mtext(expression(paste("Calibration R"^"2")), side = 1, line = 2.5, cex = o.cex)

#end plot.
dev.off()
