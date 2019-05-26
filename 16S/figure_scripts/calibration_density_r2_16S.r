rm(list=ls())
source('paths.r')

#Set output path.----
output.path <- paste0(pecan_gen_16S_dir, "figures/calibration_density_16S.png")

#load data.----
# load fits for phylogenetic groups
pl <- readRDS(paste0(scc_gen_16S_dir,"/JAGS_output/prior_phylo_JAGSfit_fewer_taxa.rds"))
phyla <- readRDS(paste0(scc_gen_16S_dir,"/JAGS_output/prior_phylo_JAGSfit_phylum_fewer_taxa_more_burnin.rds"))
pl$phylum <- phyla$phylum
# load fits for functional groups
fg_all <- readRDS(prior_16S_all.fg.groups_JAGSfits.path)

fg.r2.all <- list()
for (p in 1:length(fg_all)){
  fg <- fg_all[[p]]
  obs <- fg$no.nutr.preds$observed
  pred <- fg$no.nutr.preds$predicted
  fg.r2 <- list()
  for(i in 1:ncol(obs)){
    fg.r2[[i]] <- summary(lm(obs[,i] ~ pred[,i]))$r.squared
  }
  fg.r2 <- unlist(fg.r2)
  names(fg.r2) <- colnames(obs)
  fg.r2 <- fg.r2[names(fg.r2) != 'other'] #drop 'other'.
  fg.r2.all[[p]] <- fg.r2
}
fg.r2.all <- list(unlist(fg.r2.all))
names(fg.r2.all)<- "functional group"
#phylogenetic groups.
pl.r2 <- list()
lev.r2.out <- list()
for(i in 1:length(pl)){
  #if (p==1) lev <- pl[[i]]$
  lev <- pl[[i]]
  obs <- lev$no.nutr.preds$observed
  pred <- lev$no.nutr.preds$predicted
  lev.r2 <- list()
  for(j in 1:ncol(obs)){lev.r2[[j]] <- summary(lm(obs[,j] ~ pred[,j]))$r.squared}
  lev.r2 <- unlist(lev.r2)
  names(lev.r2) <- colnames(obs)
  lev.r2 <- lev.r2[names(lev.r2) != 'other']
  pl.r2[[i]] <- lev.r2
  names(pl.r2)[[i]] <- names(pl)[i]
}
all.r2 <- c(fg.r2.all, pl.r2)
lev.mu <- unlist(lapply(all.r2, mean))
lev.sd <- unlist(lapply(all.r2, sd))
lev.N  <- unlist(lapply(all.r2, length))
lev.se <- lev.sd / sqrt(lev.N)
all.r2.unlist <- unlist(all.r2)
#names(lev.mu)[1] <- 'functional'


#Get truncated density so that we don't have densities less than zero.----
h <- density(all.r2.unlist)$bw  #get dansity bandwith
# Compute edge weights.
w <- 1 / pnorm(0, mean=all.r2.unlist, sd=h, lower.tail=FALSE)
#Generate truncated distribution.
h <- density(all.r2.unlist, bw=h, kernel="gaussian", weights=w / length(all.r2.unlist))
h$y[h$all.r2.unlist < 0] <- 0
#Check: the integral ought to be close to 1:
if(sum(h$y * diff(h$x)[1]) > 1.1 | sum(h$y * diff(h$x)[1]) < 0.9){
  warning('Zero-truncated density sum is not close to 1. Consider doing something else...')
}


#Make density plot.----
#png save settings.
png(filename=output.path,width=8,height=5,units='in',res=300)

#global plot settings.
trans <- 0.3 #shading transparency.
o.cex <- 1.3 #outer label size.
par(mfrow = c(1,2), mar = c(5.5,4.2,1.5,1.5))

#R2 density plot.
plot(h,xlim = c(0, 1), ylim = c(0, round(max(h$y), 0)), bty = 'n', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1)
polygon(h, col = adjustcolor('purple',trans), border = NA)
mtext('Density', side = 2, line = 2.2, cex = o.cex)
mtext(expression(paste("Calibration R"^"2")), side = 1, line = 2.5, cex = o.cex)

# #calbiration R2 as a function of phylo/function scale.
x <- 1:length(lev.mu)
limy <- c(0,max(lev.mu + lev.se))
plot(lev.mu ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='n', xaxt = 'n')
#segments(x,lev.mu-lev.se,x,lev.mu+lev.se)
lines(x, lev.mu, lty = 2)
mtext(expression(paste("Calibration R"^"2")), side = 2, line = 2.2, cex = o.cex)
axis(1, labels = F)
text(x=x, y = -.04, labels= names(lev.mu), srt=45, adj=1, xpd=TRUE)
# 
# #calbiration R2 as a function of phylo/function scale, with all points shown.
# 
# df <- lapply(all.r2, data.frame)
# levs <- list()
# for (i in 1:6){
#   lev <- df[[i]]
#   lev$level <- names(df[i])
#   lev$taxon <- rownames(df[[i]])
#   lev$val <- lev$X..i..
#   lev$X..i.. <- NULL
#   lev$num <- i
#   levs[[i]] <- lev
# }
# r2.df <- dplyr::bind_rows(levs)
# r2.df$level <- as.factor(r2.df$level)
# #levels(r2.df$level) <- as.factor(unique(r2.df$level))
# 
# x <- 1:length(lev.mu)
# limy <- c(0,1)
# plot(lev.mu ~ x, cex = 2, ylim = limy, pch = 18, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.7),
#      ylab = NA, xlab = NA, bty='n', xaxt = 'n')
# 
# points(r2.df$val ~ jitter(r2.df$num,.2), cex=.6, pch = 16)
# 
# mtext(expression(paste("Calibration R"^"2")), side = 2, line = 2.2, cex = o.cex)
# axis(1, labels = F)
# text(x=x, y = -.04, labels= names(lev.mu), srt=45, adj=1.1, xpd=TRUE)
# 

#end plot.
dev.off()
