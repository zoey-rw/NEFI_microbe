#NEON plot-scale cross-validation.
#This has a massive amount of burnin, so takes a while to run.
#Fit MULTINOMIAL dirlichet models to all groups of bacteria from 95% of NEON core-scale observations.
#Not going to apply hierarchy, because it would not be a fair comparison to the Tedersoo model.
#Missing data are allowed.
#clear environment
rm(list = ls())
library(data.table)
library(doParallel)
source('paths.r')
#source('NEFI_functions/dmulti-ddirch_site.level_JAGS.r')
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/tic_toc.r')

# source function from colin's github
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions//ddirch_site.level_JAGS.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

# source Colin's paths.r
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/paths.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

#detect and register cores.----
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#set output path.----
output.path <- plot.CV_NEON_ddirch_16S_JAGSfit
calval_data.path <- plot.CV_NEON_cal.val_data_16S.path

#load NEON plot-scale data.----
dat <- readRDS(hierarch_filled_data.path)
y <- readRDS(NEON_all.phylo.levels_plot.site_obs_16S.path)

#get core-level covariate means and sd.----
core_mu <- dat$core.plot.mu
plot_mu <- dat$plot.plot.mu
site_mu <- dat$site.site.mu

#merge together.
plot_mu$siteID <- NULL
core.preds <- merge(core_mu   , plot_mu)
core.preds <- merge(core.preds, site_mu)
core.preds$relEM <- NULL
names(core.preds)[names(core.preds)=="b.relEM"] <- "relEM"

#get core-level SD.
core_sd <- dat$core.plot.sd
plot_sd <- dat$plot.plot.sd
site_sd <- dat$site.site.sd
#merge together.
plot_sd$siteID <- NULL
core.sd <- merge(core_sd   , plot_sd)
core.sd <- merge(core.sd, site_sd)
core.sd$relEM <- NULL
names(core.sd)[names(core.sd)=="b.relEM"] <- "relEM"

# change plotIDs to match those from Y obs
core.preds$plotID <- gsub('_','\\.',core.preds$plotID)
core.sd$plotID <- gsub('_','\\.',core.sd$plotID)

#Split into calibration / validation data sets.----
set.seed(420)
#ID <- rownames(yphylumplot.fit$mean)
#cal.ID <- sample(ID, round(length(ID)/ 2))
#val.ID <- ID[!(ID %in% cal.ID)]
#Subset by plot and site.
cal.p <- 0.7 #how much data in calibration vs. validation.
plotID <- rownames(y$phylum$plot.fit$mean)
siteID <- substr(plotID,1, 4)
plots <- data.frame(plotID, siteID) 
plots <- plots[plotID %in% core.preds$plotID,]
sites <- unique(plots$siteID)
cal <- list()
val <- list()
for(i in 1:length(sites)){
  sub <- plots[plots$siteID == sites[i],]
  cal_sub <- sub[sub$plotID %in% sample(sub$plotID, round(nrow(sub) * cal.p)),]
  val_sub <- sub[!(sub$plotID %in% cal_sub$plotID),]
  cal[[i]] <- cal_sub
  val[[i]] <- val_sub
}
cal <- do.call(rbind, cal)
val <- do.call(rbind, val)
cal.ID <- as.character(cal$plotID)
val.ID <- as.character(val$plotID)
#cal.ID <- cal.ID[cal.ID %in% core.preds$plotID]
#val.ID <- val.ID[val.ID %in% core.preds$plotID]


#loop through y values. Wasn't an easy way to loop through levels of list.
y.cal <- list()
y.val <- list()
for(i in 1:length(y)){
  lev <- y[[i]]$plot.fit
  lev.cal <- list()
  lev.val <- list()
  lev.cal$mean <- lev$mean[rownames(lev$mean) %in% cal.ID,]
  lev.val$mean <- lev$mean[rownames(lev$mean) %in% val.ID,]
  lev.cal$lo95 <- lev$lo95[rownames(lev$lo95) %in% cal.ID,]
  lev.val$lo95 <- lev$lo95[rownames(lev$lo95) %in% val.ID,]
  lev.cal$hi95 <- lev$hi95[rownames(lev$hi95) %in% cal.ID,]
  lev.val$hi95 <- lev$hi95[rownames(lev$hi95) %in% val.ID,]
  #return to larger list.
  y.cal[[i]] <- lev.cal
  y.val[[i]] <- lev.val
}
names(y.cal) <- names(y)
names(y.val) <- names(y)

#split x means and sd's.
x_mu.cal <- core.preds[core.preds$plotID %in% cal.ID,]
x_mu.val <- core.preds[core.preds$plotID %in% val.ID,]
x_sd.cal <- core.sd   [core.sd   $plotID %in% cal.ID,]
x_sd.val <- core.sd   [core.sd   $plotID %in% val.ID,]

#match the order.
x_mu.cal <- x_mu.cal[order(match(x_mu.cal$plotID, rownames(y.cal$phylum$mean))),]
x_sd.cal <- x_sd.cal[order(match(x_sd.cal$plotID, rownames(y.cal$phylum$mean))),]
x_mu.val <- x_mu.val[order(match(x_mu.val$plotID, rownames(y.val$phylum$mean))),]
x_sd.val <- x_sd.val[order(match(x_sd.val$plotID, rownames(y.val$phylum$mean))),]

#subset to predictors of interest, drop in intercept.----
rownames(x_mu.cal) <- rownames(y.cal$phylum$abundances)
intercept <- rep(1, nrow(x_mu.cal))
x_mu.cal <- cbind(intercept, x_mu.cal)
x_mu.cal <- x_mu.cal[,c('intercept','pH','pC','cn','relEM','map','mat','NPP','forest','conifer')]


#save calibration/valiation data sets.----
dat.cal <- list(y.cal, x_mu.cal, x_sd.cal)
dat.val <- list(y.val, x_mu.val, x_sd.val)
names(dat.cal) <- c('y.cal','x_mu.cal','x_sd.cal')
names(dat.val) <- c('y.val','x_mu.val','x_sd.val')
dat.out <- list(dat.cal, dat.val)
names(dat.out) <- c('cal','val')
saveRDS(dat.out, calval_data.path)

#fit model using function in parallel loop.-----
#for running production fit on remote.
cat('Begin model fitting loop...\n')
tic()
output.list<-
  foreach(i = 1:length(y)) %dopar% {
    #y.group <- round(y.cal[[i]]$mean * 3000) #Should perhaps draw from uncertainties, but these supplied hi/lo95 values dont account for covariance among taxa. # <- what
    y.group <- y.cal[[i]]$mean
  
    fit <- site.level_dirlichet_jags(y=y.group,x_mu=x_mu.cal, x_sd=x_sd.cal, #seq.depth = rowSums(y.group),
                                        #adapt = 2000, burnin = 20000, sample = 6000, 
                                        adapt = 200, burnin = 200, sample = 200,   #testing
                                        parallel = T, parallel_method = 'parallel') #setting parallel rather than rjparallel. 
    return(fit)                                                                     #allows nested loop to work.
  }
cat('Model fitting loop complete! ')
toc()


#name the items in the list
names(output.list) <- names(y.cal)

#save output.----
cat('Saving fit...\n')
saveRDS(output.list, output.path)
cat('Script complete. \n')