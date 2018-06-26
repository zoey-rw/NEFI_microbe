#fitting semivariograms to data.
#clear environment, load pacakges
rm(list=ls())
library (sp)
library(gstat)
library(boot)
library(data.table)

#Load Talbot 2014 PNAS data.
d <- data.table(readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/tal_core_data.rds'))
d$relEM <- (d$relEM * (length(d$relEM) - 1) + 0.5) / length(d$relEM)
d$relEM <- logit(d$relEM)
d <- d[,.(relEM,pH, moisture, cn, mat30, map30, doy, horizon, epoch_date, site, latitude, longitude)]
d <- d[complete.cases(d),]

#set coordinates
coordinates(d) <- ~longitude + latitude
bubble(d, zcol = 'relEM', fill=T, do.sqrt = F, maxsize=3)

#fit variogram
vario.1 <- variogram(relEM~1, data = d) #get semivariance values
vario.model.1 <- fit.variogram(vario.1, vgm(model='Sph', nugget=T)) #fit a variogram model with a Spherical error structure and a "nugget" (nugget = intercept)
#plot your semivariance and your model.
plot(vario.1, model=vario.model.1)

#We can repeat this with the residuals of a fitted model.
m <-  lm(relEM ~ pH + cn + moisture + mat30 + map30 + horizon + doy, data = d)
summary(m)
d$resid <- residuals(m)

#fit variogram to residuals
vario.2 <- variogram(resid~1, data = d) #get semivariance values
vario.model.2 <- fit.variogram(vario.2, vgm(model='Sph', nugget=T)) #fit a variogram model with a Spherical error structure and a "nugget" (nugget = intercept)
#plot your semivariance and your model.
<<<<<<< HEAD
plot(vario.2, model=vario.model.1)
=======
plot(vario.2, model=vario.model.1)
>>>>>>> 86f6aad5e093b910fe7cc8ff293b9089204f114c
