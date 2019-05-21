#compare out of sample forecast to number of SVs within a group.
rm(list=ls())
source('paths.r')
library(data.table)

#load PHYLO data.----
#load in-/out-of-sample diversity stats.
cal.div <- readRDS(bahram_16S_common_phylo_groups_list.path)
val.div <- readRDS(NEON_16S_phylo_groups_abundances.path)

#Load in calibration/validation models.
cal <- readRDS(bahram_16S_prior_phylo.group_JAGSfits)
val <- readRDS(NEON_cps_fcast_all_phylo_16S.path)

#Load in validation data.
val.obs <- readRDS(NEON_all.phylo.levels_plot.site_obs_16S.path)


# 
# #load FG data.----
# #load in-/out-of-sample diversity stats.
# cal.div <- readRDS(prior_fg_abundances_16S.path)
# val.div <- readRDS(NEON_fg_abundances_16S.path)
# 
# #Load in calibration/validation models.
# cal <- readRDS(prior_16S_all.fg.groups_JAGSfits.path)
# val <- readRDS(NEON_cps_fcast_fg_16S.path)
# 
# #Load in validation data.
# val.obs <- readRDS(NEON_all.fg_plot.site_obs_16S.path)
# 






#Get calibration R2 values, merge in diversity stats.----
cal.r2_table <- list()
for(i in 1:length(cal)){
  cal_level <- cal[[i]]
  div <- cal.div[[i]]$group_frequencies
  pred <- cal_level$predicted
  obs <- cal_level$observed
  cal.lev_r2   <- list()
  cal.lev_r2.1 <- list()
  for(j in 1:ncol(pred)){
    mod <- lm(( obs[,j]) ~ (pred[,j]))
    cal.lev_r2[[j]] <- summary(mod)$r.squared
    #relative to 1:1 line.
    rss <- sum((pred[,j] - obs[,j]) ^ 2)       ## residual sum of squares
    tss <- sum((obs[,j] - mean(obs[,j])) ^ 2)  ## total sum of squares
    rsq.1 <- 1 - rss/tss
    rsq.1 <- ifelse(rsq.1 < 0, 0, rsq.1)
    cal.lev_r2.1[[j]] <- rsq.1
  }
  cal.lev_r2 <- data.frame(unlist(cal.lev_r2))
  cal.lev_r2 <- cbind(colnames(obs), cal.lev_r2)
  colnames(cal.lev_r2) <- c('group','cal_rsq')
  cal.lev_r2$cal_rsq.1 <- unlist(cal.lev_r2.1)
  #merge in diversity stats.
  colnames(div) <- c('group', 'cal_samp_freq','cal_N.SV','cal_diversity','cal_evenness')
  cal.lev_r2 <- merge(cal.lev_r2, div, all.x = T)
  cal.lev_r2 <- cal.lev_r2[-grep('other',cal.lev_r2$group),]
  cal.r2_table[[i]] <- cal.lev_r2
}
names(cal.r2_table) <- names(cal)

#get validation R2 values (site-level), merge in diversity stats.----
val.r2_table <- list() #full table
for(i in 1:length(val)){
  pred <-     val[[i]]$site.fit$mean
  obs <- val.obs[[i]]$site.fit$mean
  pred <- pred[,order(match(colnames(pred), colnames(obs)))]
  
  #make sure same order and identity of rownames
  obs <-  obs[rownames(obs ) %in% rownames(pred),]
  pred <- pred[rownames(pred) %in% rownames( obs),]
  rsq   <- list()
  rsq.1 <- list()
  abundance <- list()
  variance  <- list()
  range     <- list()
  for(j in 1:ncol(obs)){
    #get best-fit rsq.
    mod <- lm(obs[,j] ~ pred[,j])
    rsq[[j]] <- summary(mod)$r.squared
    #get rsq relative to 1:1 line.
    rss <- sum((pred[,j] - obs[,j]) ^ 2)       ## residual sum of squares
    tss <- sum((obs[,j] - mean(obs[,j])) ^ 2)  ## total sum of squares
    rsq_1 <- 1 - (rss/tss)
    rsq_1 <- ifelse(rsq_1 < 0, 0, rsq_1)
    rsq.1[[j]] <- rsq_1
    abundance[[j]] <- boot::inv.logit(mean(boot::logit(obs[,j])))
    variance[[j]] <- boot::inv.logit(  sd(boot::logit(obs[,j])))
    range[[j]] <- max(obs[,j]) - min(obs[,j])
  }
  rsq <- data.frame(unlist(rsq))
  rsq.1 <- data.frame(unlist(rsq.1))
  abundance <- data.frame(unlist(abundance))
  variance <- data.frame(unlist(variance))
  range <- data.frame(unlist(range))
  rsq <- cbind(colnames(obs), rsq, rsq.1,abundance, variance,range)
  colnames(rsq) <- c('group','val_rsq','val_rsq.1','val_abundance','val_variance','val_range')
  
  rsq$group <- as.character(rsq$group)
  div <- val.div[[i]]$group_frequencies
  colnames(div) <- c('group','val_samp_freq','val_N.SV','val_diversity','val_evenness')
  rsq <- rsq[-grep('other',rsq$group),]
  rsq <- merge(rsq, div, all.x = T)
  rsq$phylo_level <- names(val)[i]
  val.r2_table[[i]] <- rsq
}
names(val.r2_table) <- names(val)

#merge the tables.----
cal.stats <- do.call(rbind, cal.r2_table)
val.stats <- do.call(rbind, val.r2_table)
d <- merge(cal.stats, val.stats)

#Model the R2 values.----
#vlaidation rsq of taxa with cal_rsq > 0.10.
#calibration and validation rsq are not related.
#validation rsq is related to validation variance, N.SVs (which are correlated). Best model just includes val_variance.
#diversity, N.SVs, variance all correlated, variance best predictor.
#also related to phylo level, but again, this mkaes val_variance model worse.
#true whether we model rsq or rsq.1 
#--------
#### all of the above was for fungi, not bacteria ####
mod <- lm(val_rsq ~ cal_rsq + phylo_level,data = d[d$cal_rsq > 0.33,])
summary(mod)

par(mfrow = c(2,2)) # no clear patterns
plot(val_rsq ~ cal_rsq        , data = d[d$cal_rsq > 0.1,])
plot(val_rsq ~ log10(val_N.SV), data = d[d$cal_rsq > 0.1,])
plot(val_rsq ~ val_diversity  ,data = d[d$cal_rsq > 0.1,])
plot(val_rsq ~ val_variance   , data = d[d$cal_rsq > 0.1,])

summary(lm(val_rsq ~ cal_diversity, data = d))
summary(lm(val_diversity ~ val_variance, data = d[d$cal_rsq > 0.1,]))
summary(lm(val_rsq ~ log10(cal_N.SV), data = d[d$cal_rsq > 0.1,]))
summary(lm(val_rsq ~ log10(cal_N.SV)+val_variance, data = d[d$cal_rsq > 0.1,]))
summary(lm(val_rsq ~ val_variance, data = d[d$cal_rsq > 0.1,])) #p=.0003, adjR2=.1118
summary(lm(val_rsq ~ phylo_level, data = d[d$cal_rsq > 0.1,])) 
summary(lm(val_rsq ~ val_samp_freq, data = d)) 
summary(lm(val_rsq ~ phylo_level, data = d[d$cal_rsq > 0.1,])) #p=

summary(lm(val_rsq ~ phylo_level, data = d))

plot(val_rsq ~ log10(cal_N.SV), data = d[d$cal_rsq > 0.1,])
plot(val_rsq ~ cal_diversity, data = d[d$cal_rsq > 0.1,])
h <- d[d$cal_rsq > 0.1,]
x <- h[,!colnames(h) %in% c("val_rsq", "group")]

library(MASS)
fit <- lm(h$val_rsq~.,data=x)
step <- stepAIC(fit, direction="both")
summary(step)
step$anova # display results


summary(lm(d$val_rsq ~ cal_rsq.1 + val_rsq.1 + val_abundance + val_variance + 
             val_range + val_samp_freq + val_N.SV + val_diversity + val_evenness, data=d)) # best model according to stepAIC, but doesn't account for all the collinearity

library(corrplot) # correlation plot
library(caret)
dm = data.matrix(d[,!colnames(d) %in% c("val_rsq", "val_rsq1","phylo_level", "group")]) 
corrplot(cor(dm), order = 'hclust')
findCorrelations(cor(dm), 0.75)


# can we see why phylum forecasts split neatly into good and bad???
# no. n=11 is too small?
p <- d[d$phylo_level=="phylum",]
t <- p[,!colnames(p) %in% c("val_rsq", "phylo_level", "group")]
summary(lm(val_rsq ~ val_samp_freq, data = p)) #p=.46
fit <- lm(p$val_rsq~.,data=t)
step <- stepAIC(fit, direction="both")
summary(step)
vif(fit)
sqrt(vif(fit)) > 2

# no phyla
np <- d[d$phylo_level != "phylum",]
head(np)
summary(lm(val_rsq ~ cal_rsq, data = np[np$cal_rsq > 0.33,])) #p=.46


#THIS IS THE ANSWER. 
#When you subset to calibration R2 > 0.33, calibration does predict validation, R2 = 0.49, unadjusted. <- not true for bacteria at all :(
par(mfrow=(c(1,2)))
mod1 <- lm(val_rsq ~ cal_rsq, data = d[d$cal_rsq > 0.33,])
plot(val_rsq ~ cal_rsq, data = d[d$cal_rsq > 0.33,], ylab = 'validation R2', xlab = 'calibration R2')
abline(mod1, lwd = 2)
mtext(paste0('R2 = ',round(summary(mod1)$r.squared,2)), side = 3, line = -2, adj = 0.05)

#You can explain 88% (R2 = 0.88, unadjusted) if you add phylo/functional group to the model. <- also not true for bacteria
mod2 <- lm(val_rsq ~ cal_rsq + phylo_level, data = d[d$cal_rsq > 0.33,])
y <- d[d$cal_rsq > 0.33,]$val_rsq
x <- fitted(mod2)
plot(y ~ x, ylab = 'observed validation R2', xlab = 'predicted validation R2')
abline(lm(y~x), lwd = 2)
mtext(paste0('R2 = ',round(summary(mod2)$r.squared,2)), side = 3, line = -2, adj = 0.05)
