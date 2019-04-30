#calculating 16S spatial signal in taxa by phylo scale across NEON, using entire network.

rm(list=ls())
source('paths.r')
library(boot)
library(ape)

#logit transform observed values and model residuals?----
do_logit <- F

#Load data.----
# soil core data
loc_all <- readRDS(dp1.10086.00_output_16S.path)
# abundance data
fg <- readRDS(NEON_all.fg_plot.site_obs_16S.path)
d_all <- readRDS(NEON_all.phylo.levels_plot.site_obs_16S.path)
d_all <- c(fg, d_all)
#spatial forecast.
fcast_fg <- readRDS(NEON_cps_fcast_fg_16S.path)
fcast <- readRDS(NEON_cps_fcast_all_phylo_16S.path)
fcast <- c(fcast_fg, fcast)
# read in obs table that links deprecatedVialID and geneticSampleID
#map <- readRDS(obs.table_16S.path)
map <- readRDS(core_obs_16S.path)

#Calculate spatial statistics for all levels on forecast residuals.----
a_out <- list() #spatial signal of raw values.
b_out <- list() #spatial signal of model residuals.
for(k in 1:length(d_all)){
  #Grab tax/functional level.
  d <- d_all[[k]] 
  d <- d$core.fit
  loc <- loc_all
  pred <- fcast[[k]]$core.fit$mean
  
  #match up data, subset to mineral soil, drop stuff from identical locations.----
  map <- map[,c("deprecatedVialID", "geneticSampleID")]
  map$geneticSampleID <- gsub('-GEN','',map$geneticSampleID)
  d <- as.data.frame(d)
  d$deprecatedVialID <- rownames(d)
  d <- merge(d, map, by = "deprecatedVialID")
  rownames(d) <- d$geneticSampleID
  d$geneticSampleID <- NULL
  d$deprecatedVialID <- NULL
  #drop organic horizons. Generates zero distances.
  loc <- loc[loc$horizon == 'M',]
  loc$lat_lon <- paste0(loc$adjDecimalLongitude, loc$adjDecimalLatitude)
  drop <- loc[duplicated(loc$lat_lon),]$sampleID #drop some duplicated locations.
  loc <- loc[!(loc$sampleID %in% drop),]
  loc <- loc[loc$sampleID %in% rownames(d),]
  d <-   d[rownames(d) %in% loc$sampleID,]
  loc <- loc[order(match(loc$sampleID, rownames(d))),]

  #make sure pred/obs are in both dataframes
  pred <- pred[rownames(pred) %in% rownames(d),]
  d <- d[rownames(d) %in% rownames(pred),]
  pred <- pred[colnames(pred) %in% colnames(d),]
  d <- d[colnames(d) %in% colnames(pred),,drop=FALSE]
  #make sure predictions and observation columns are in the same order.
  pred <- pred[,order(match(colnames(pred), colnames(d)))]
  d <- d[,order(match(colnames(d), colnames(pred)))]
  loc <- loc[loc$sampleID %in% rownames(d),]
  
  #remove "other" column
  pred <- pred[,colnames(pred) != "other", drop=FALSE]
  d <- d[,colnames(d) != "other", drop=FALSE]
  
  #generate model residuals as a logit difference.
  resid <- (d) - (pred)
  if(do_logit == T){
    resid <- logit(d) - logit(pred)
    d <- logit(d)
  }
  
  #Generate spatial distance matrix.----
  tax.loc <- loc[,c('adjDecimalLongitude','adjDecimalLatitude')]
  tax.dist <- geosphere::distm(tax.loc)
  inv.tax.dist <- 1/tax.dist
  diag(inv.tax.dist) <- 0
  
  #calulate moran statistics.----
  moran_out <- list()
  resid_out <- list()
  for(i in 1:ncol(d)){
    #Get r.abundance distance matrix on logit scale.
    cat('Calculating ',colnames(d)[i],' moran statistics...\n')
    #make sure there is variance in taxa being calculated.
    if(length(unique(d[,i])) == 1){
      moran_out[[i]] <- NA
      next
    }
    #moran stat for raw residuals.
    moran.stat  <- ape::Moran.I((d[,i]), inv.tax.dist)
    moran.stat  <- c(moran.stat$observed, moran.stat$p.value)
    #moran stat for post-model residuals.
    resid.stat <- ape::Moran.I(resid[,i], inv.tax.dist)
    resid.stat <- c(resid.stat$observed, resid.stat$p.value)
    #moran.resid
    names(moran.stat) <- c('moran','p.val')
    names(resid.stat) <- c('moran','p.val')
    moran_out[[i]] <- moran.stat
    resid_out[[i]] <- resid.stat
  }
  
  #Drop taxa that couldn't fit, wrap up level output.----
  names(moran_out) <- colnames(d)
  names(resid_out) <- colnames(resid)
  moran_out <- moran_out[!(is.na(moran_out))]
  resid_out <- resid_out[!(is.na(resid_out))]
  moran_out <- data.frame(do.call(rbind, moran_out))
  resid_out <- data.frame(do.call(rbind, resid_out))
  
  #return output.----
  a_out[[k]] <- moran_out
  b_out[[k]] <- resid_out
  
}
names(a_out) <- names(d_all)
names(b_out) <- names(d_all)

# flatten 12 functional group models into one list item
fg_a <- a_out[1:12]
fg_a_group <- dplyr::bind_rows(fg_a)
rownames(fg_a_group) <- c(names(fg_a)[1:11], "copiotroph", "oligotroph") #so cop/olig are separate
fg_a_group <- list(fg_a_group)
names(fg_a_group) <- "functional"
a_out <- c(fg_a_group, a_out[13:17]) # recombine fg with phylo groups

fg_b <- b_out[1:12]
fg_b_group <- dplyr::bind_rows(fg_b)
rownames(fg_b_group) <- c(names(fg_b)[1:11], "copiotroph", "oligotroph") #so cop/olig are separate
fg_b_group <- list(fg_b_group)
names(fg_b_group) <- "functional"
b_out <- c(fg_b_group, b_out[13:17]) # recombine fg with phylo groups

#calculate average morans I for all groups.----
a_avg <- list()
b_avg <- list()
for(i in 1:length(a_out)){
  #no model.
  a.mu <- mean(a_out[[i]]$moran)
  a.sd <-   sd(a_out[[i]]$moran)
  a.se <- a.sd/sqrt(nrow(a_out[[i]]))
  to_return.a <- c(a.mu, a.sd, a.se)
  names(to_return.a) <- c('mu','sd','se')
  a_avg[[i]] <- to_return.a
  #residuals of model.
  b.mu <- mean(b_out[[i]]$moran)
  b.sd <-   sd(b_out[[i]]$moran)
  b.se <- b.sd/sqrt(nrow(b_out[[i]]))
  to_return.b <- c(b.mu, b.sd, b.se)
  names(to_return.b) <- c('mu','sd','se')
  b_avg[[i]] <- to_return.b
  
}
a_avg <- data.frame(do.call(rbind, a_avg))
b_avg <- data.frame(do.call(rbind, b_avg))
rownames(a_avg) <- names(a_out)
rownames(b_avg) <- names(b_out)

#assign x positions, functional groups first.
a_avg$x <- c(1:(nrow(a_avg)))
b_avg$x <- c(1:(nrow(b_avg)))

#plot.-----
par(mfrow=c(1,2))

#Raw spatial signal.----
limy <- c(0, max(b_avg$mu + b_avg$se))
plot(mu ~ x, data = a_avg, cex = 2, pch = 16, ylim = limy,
     ylab = "Moran's I", xlab = NA, xaxt = 'n')
#error bars.
mu <- a_avg$mu
x <- a_avg$x
upr <- mu + a_avg$se
lwr <- mu - a_avg$se
arrows(c(x), lwr, c(x), upr, length=0.00, angle=90, code=3, col = 'black', lwd = 2)
#x-axis.
x.lab <- rownames(a_avg)
x.lab[x.lab == 'fg'] <- 'function_group'
axis(1, at=a_avg$x, labels=x.lab, cex = 0.8)
mtext('Spatial Signal',side = 3, line = 0.5, cex = 1.4)

#Post model spatial signal.----
plot(mu ~ x, data = b_avg, cex = 2, pch = 16, ylim = limy,
     ylab = "Moran's I", xlab = NA, xaxt = 'n')
#error bars.
mu <- b_avg$mu
x <- b_avg$x
upr <- mu + b_avg$se
lwr <- mu - b_avg$se
arrows(c(x), lwr, c(x), upr, length=0.00, angle=90, code=3, col = 'black', lwd = 2)
#x-axis.
x.lab <- rownames(b_avg)
x.lab[x.lab == 'fg'] <- 'function_group'
axis(1, at=b_avg$x, labels=x.lab, cex = 0.8)
mtext('Spatial Signal Model Residuals',side = 3, line = 0.5, cex = 1.4)
