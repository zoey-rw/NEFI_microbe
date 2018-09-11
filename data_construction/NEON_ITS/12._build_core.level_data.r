#build core-level observation data set.
#clear environment, source paths.
rm(list=ls())
library(runjags)
source('paths.r')

#1. observations made at the core level.----
#samples with ITS data, along with DNA identifiers. Only take samples measured during peak greeness.
dp1.10801 <- readRDS(dp1.10108.00_output.path)
dp1.10801$geneticSampleID <- as.character(dp1.10801$geneticSampleID)
dp1.10801 <- dp1.10801[!(dp1.10801$geneticSampleID == ''),] #many samples missing geneticSampleID. throw them out.
dp1.10801 <- dp1.10801[order(dp1.10801$geneticSampleID),]

#moisture and pH data.
dp1.10086 <- readRDS(dp1.10086.00_output.path)
dp1.10086$site_date <- paste0(dp1.10086$siteID,'-',dp1.10086$dateID)
dp1.10086$geneticSampleID <- as.character(dp1.10086$geneticSampleID)

#soil C and N data.
dp1.10078 <- readRDS(dp1.10078.00_output.path)
dp1.10078$site_date_plot <- paste0(dp1.10078$siteID,'-',dp1.10078$dateID,'-',dp1.10078$plotID)

#merge these 3 sets of observations together.
to_merge <- dp1.10086[,!c(colnames(dp1.10086) %in% colnames(dp1.10801))]
to_merge$geneticSampleID <- dp1.10086$geneticSampleID
merged <- merge(dp1.10801,to_merge)
to_merge <- dp1.10078[,!c(colnames(dp1.10078) %in% colnames(merged))]
to_merge$sampleID <- dp1.10078$sampleID
to_merge <- to_merge[!(duplicated(to_merge$sampleID)),] #get rid of analytical replicates.
merged <- merge(merged,to_merge, by = 'sampleID', all.x=T)
merged$year <- substring(merged$dateID,1,4)
merged <- merged[!is.na(merged$siteID),]

#Peak Greenness in 2014. 551 observations. 13 sites.
merged <- merged[merged$sampleTiming == 'peakGreenness' & merged$year == '2014',]

#2. Get observation uncertainty for soil C and N, and their ratio.----
duped <- dp1.10078[duplicated(dp1.10078$sampleID),]
duped <- dp1.10078[dp1.10078$sampleID %in% duped$sampleID,]

#fit an observation uncertainty.
pC <- duped[,c('sampleID','organicCPercent')]
to_ag <- aggregate(organicCPercent ~ sampleID, data = pC, FUN = mean)
colnames(to_ag)[2] <- 'mean_C'
to_analyze <- merge(to_ag, pC[,c('sampleID','organicCPercent')], all.y=T)
plot(organicCPercent ~ mean_C, data = to_analyze)

#jags model
jags.model = "
model {
  #prior on relationship between mean and observed pC values.
  #preliminary analysis says this is slope of 1, intercept of zero, as it should be. Not fitting an intercept.
    m ~ dnorm(0, .0001) #flat uninformative priors.
  
  #magnitude of %C observation affects uncertainty.
  #tau is still a function of sigma.
  #sigma is now a funtion of an intercept (as before) and magnitude of C observation, with a parameter that gets fitted.
  k0 ~ dunif(0,100)
  k1 ~ dunif(0,100)
  
  #fit a linear model. x is mean C for an observation, y is observed.
  for(i in 1:N){
    y.hat[i] <-  m*x[i]
    sigma[i] <- k0 + k1*x[i]
      tau[i] <- pow(sigma[i], -2)
        y[i] ~ dnorm(y.hat[i], tau[i])
  }
} #close model loop.
"
jd <- list(y=to_analyze$organicCPercent, x=to_analyze$mean_C, N = nrow(to_analyze))
fit   <- run.jags(model = jags.model,
                 data = jd,
                 adapt = 200,
                 burnin = 2000,
                 sample = 3000,
                 monitor = c('m','k0','k1')
                  )
out <- summary(fit)
