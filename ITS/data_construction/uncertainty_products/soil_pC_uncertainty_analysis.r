#building uncertainty products for NEON soils data where available.
#clear environment, load paths.
rm(list=ls())
source('paths.r')
library(runjags)
#load data.
dp1.10078 <- readRDS(dp1.10078.00_output.path)

#variable generating ucnertainty for.
of_interest <- 'organicCPercent'
#set output path.
output.path <- pC_NEON_uncertainty_model.path

#organize analysis product.----
#grab within core replicated observations.
duped <- dp1.10078[duplicated(dp1.10078$sampleID),]
duped <- dp1.10078[dp1.10078$sampleID %in% duped$sampleID,]
duped <- duped[,c('sampleID',of_interest)]
#take core-level means.
to_ag <- aggregate(duped[,c(of_interest)] ~ sampleID, data = duped, FUN = mean)
colnames(to_ag)[2] <- 'mean'
to_analyze <- merge(to_ag, duped[,c('sampleID',of_interest)], all.y=T)
colnames(to_analyze)[3] <- 'obs'
to_analyze <- to_analyze[complete.cases(to_analyze),]
plot(obs ~ mean, data = to_analyze)

#jags model.----
jags.model = "
model {
#prior on relationship between mean and observed pC values.
#preliminary analysis says this is slope of 1, intercept of zero, as it should be. Not fitting an intercept.
m ~ dnorm(0, .0001) #flat uninformative priors.

#magnitude of %C observation affects uncertainty.
#tau is still a function of sigma.
#sigma is now a funtion of an intercept (as before) and magnitude of C observation, with a parameter that gets fitted.
#interval censored. variances cannot be zero.
k0 ~ dnorm(0,1E-3) I(0,)
k1 ~ dnorm(0,1E-3) I(0,)

#fit a linear model. x is mean C for an observation, y is observed.
for(i in 1:N){
y.hat[i] <-  m*x[i]
sigma[i] <- k0 + k1*x[i]
tau[i] <- pow(sigma[i], -2)
y[i] ~ dnorm(y.hat[i], tau[i])
}
} #close model loop.
"
#specify jags data object and fit.----
jd <- list(y=to_analyze$obs, x=to_analyze$mean, N = nrow(to_analyze))
library(runjags)
fit  <- run.jags(model = jags.model,
                 data = jd,
                 adapt = 200,
                 burnin = 2000,
                 sample = 3000,
                 monitor = c('m','k0','k1'),
                 n.chains = 3)
out <- summary(fit)

#Save output.
saveRDS(out,output.path)
