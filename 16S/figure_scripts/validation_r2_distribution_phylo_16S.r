# distribution of validation R2 by functional and taxonomic group - 
# right now, just phylogenetic

rm(list=ls())
source('paths.r')
#source('NEFI_functions/zero_truncated_density.r')

# source ddirch_forecast
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/zero_truncated_density.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

output.path <- paste0(pecan_gen_16S_dir,"figures/valR2_by.phylo_density_fig.png")

# load in phylo fits
pl <- readRDS(paste0(scc_gen_16S_dir, "JAGS_output/prior_phylo_JAGSfit_fewer_taxa.rds"))
phylum.mod <- readRDS(paste0(scc_gen_16S_dir,"JAGS_output/prior_phylo_JAGSfit_phylum_fewer_taxa_more_burnin.rds"))
pl$phylum <- phylum.mod$phylum

map <- readRDS(core_obs_16S.path) # linking core IDs

#Get calibration r2 values.----
pl.r2 <- list()
lev.r2.out <- list()
for(i in 1:length(pl)){
  lev <- pl[[i]]
  lev <- lev$no.nutr.preds
  obs <- lev$observed
  pred <- lev$predicted
  lev.r2 <- list()
  for(j in 1:ncol(obs)){lev.r2[[j]] <- summary(lm(obs[,j] ~ pred[,j]))$r.squared}
  lev.r2 <- unlist(lev.r2)
  names(lev.r2) <- colnames(obs)
  lev.r2 <- lev.r2[names(lev.r2) != 'other']
  pl.r2[[i]] <- lev.r2
  #names(pl.r2)[[i]] <- names(pl)[i]
}
all.r2 <- unlist(pl.r2)

#load forecasts predicted and observed.----
pl.cast <- readRDS(NEON_cps_fcast_all_phylo_16S.path)
#pl.cast <- readRDS(paste0(pecan_gen_16S_dir,"/NEON_forecast_data/NEON_fcast_comp_cases.rds"))

pl.truth <- readRDS(NEON_all.phylo.levels_plot.site_obs_16S.path)
pl.core <- readRDS(NEON_16S_phylo_groups_abundances.path)

#Re-order.
# pl.cast <- pl.cast[c('fg','phylum','class','order','family','genus')]
# pl.truth <- pl.truth[c('fg','phylum','class','order','family','genus')]
# pl.core <- pl.core[c('fg','phylum','class','order','family','genus')]

#names.
#names(pl.cast)[1] <- 'functional'

#get core, plot site R2 values out of sample.----
all.core.rsq <- list()
all.plot.rsq <- list()
all.site.rsq <- list()
for(i in 1:length(pl.cast)){
  fcast <- pl.cast[[i]]
  core.rsq <- list()
  plot.rsq <- list()
  site.rsq <- list()
  #core.level----
  y <- (pl.core[[i]]$abundances + 1)/pl.core[[i]]$seq_total
  x <- fcast$core.fit$mean
  #make sure row and column orders match.
  y <- as.data.frame(y)
  y$deprecatedVialID <- rownames(y)
  y <- merge(y, map[,c("deprecatedVialID", "geneticSampleID")], by = "deprecatedVialID")
  rownames(y) <- gsub('-GEN','',y$geneticSampleID)
  y <- y[rownames(y) %in% rownames(x),]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[,colnames(y) %in% colnames(x)]
  x <- x[,colnames(x) %in% colnames(y)]
  x <- x[order(match(rownames(x),rownames(y))),]
  x <- x[,order(match(colnames(x),colnames(y)))]
  
  
  #fit model, grab r2.
  for(k in 1:ncol(x)){
    fungi_name <- colnames(x)[k]
    rsq <- summary(lm(y[,k] ~ x[,k]))$r.squared
    names(rsq) <- fungi_name
    core.rsq[[k]] <- rsq
  }
  
  #plot.level----
  x <- fcast$plot.fit$mean
  y <- pl.truth[[i]]$plot.fit$mean
  #make sure row and column order match.
  rownames(y) <- gsub('\\.','_',rownames(y))
  #rownames(y) <- gsub('-GEN','',rownames(y))
  y <- y[rownames(y) %in% rownames(x),]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[,colnames(y) %in% colnames(x)]
  x <- x[,colnames(x) %in% colnames(y)]
  x <- x[order(match(rownames(x),rownames(y))),]
  x <- x[,order(match(colnames(x),colnames(y)))]
  #fit model, grab r2.
  for(k in 1:ncol(y)){
    fungi_name <- colnames(x)[k]
    rsq <- summary(lm(y[,k] ~ x[,k]))$r.squared
    names(rsq) <- fungi_name
    plot.rsq[[k]] <- rsq
  }
  #site.level----
  x <- fcast$site.fit$mean
  y <- pl.truth[[i]]$site.fit$mean
  #make sure row and column order match.
  rownames(y) <- gsub('-GEN','',rownames(y))
  y <- y[rownames(y) %in% rownames(x),]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[,colnames(y) %in% colnames(x)]
  x <- x[,colnames(x) %in% colnames(y)]
  x <- x[order(match(rownames(x),rownames(y))),]
  x <- x[,order(match(colnames(x),colnames(y)))]
  #fit model, grab r2.
  for(k in 1:ncol(y)){
    fungi_name <- colnames(x)[k]
    rsq <- summary(lm(y[,k] ~ x[,k]))$r.squared
    names(rsq) <- fungi_name
    site.rsq[[k]] <- rsq
  }
  #wrap up for return.----
  all.core.rsq[[i]] <- unlist(core.rsq)
  all.plot.rsq[[i]] <- unlist(plot.rsq)
  all.site.rsq[[i]] <- unlist(site.rsq)
  
}
lev.mu <- unlist(lapply(all.site.rsq, mean  ))
lev.sd <- unlist(lapply(all.site.rsq, sd    ))
lev.N  <- unlist(lapply(all.site.rsq, length))
lev.se <- lev.sd / sqrt(lev.N)
names(lev.mu) <- names(pl.cast)
core.rsq <- unlist(all.core.rsq)
plot.rsq <- unlist(all.plot.rsq)
site.rsq <- unlist(all.site.rsq)
core.rsq <- core.rsq[-grep('other',names(core.rsq))]
plot.rsq <- plot.rsq[-grep('other',names(plot.rsq))]
site.rsq <- site.rsq[-grep('other',names(site.rsq))]

#Subset to observations that have a minimum calibration R2 value.----
pass <- all.r2[all.r2 > .1]
core.rsq <- core.rsq[names(core.rsq) %in% names(pass)]
plot.rsq <- plot.rsq[names(plot.rsq) %in% names(pass)]
site.rsq <- site.rsq[names(site.rsq) %in% names(pass)]
core.d <- zero_truncated_density(core.rsq)
plot.d <- zero_truncated_density(plot.rsq)
site.d <- zero_truncated_density(site.rsq)

#png save line.----
png(filename=output.path,width=8,height=5,units='in',res=300)

#global plot settings.----
par(mfrow = c(1,2))
limx <- c(0,1)
limy <- c(0, 3)
trans <- 0.2 #shading transparency.
o.cex <- 1.3 #outer label size.
cols <- c('purple','cyan','yellow')
par(mfrow = c(1,2), mar = c(5,4.2,1.5,1.5))

#Density plot.----
plot(site.d,xlim = c(0, 1), ylim = limy, bty = 'n', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
polygon(site.d, col = adjustcolor(cols[1],trans))
mtext('Density', side = 2, line = 2.2, cex = o.cex)
mtext(expression(paste("Validation R"^"2")), side = 1, line = 2.5, cex = o.cex)

#Validation rsq ~ function/phylo scale.----
x <- 1:length(lev.mu)
limy <- c(0,max(lev.mu + lev.se))
plot(lev.mu ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='n', xaxt = 'n')
segments(x,lev.mu-lev.se,x,lev.mu+lev.se)
lines(x, lev.mu, lty = 2)
mtext(expression(paste("Validation R"^"2")), side = 2, line = 2.2, cex = o.cex)
axis(1, labels = F)
text(x=x, y = -.06, labels= names(lev.mu), srt=45, adj=1, xpd=TRUE)
dev.off()
