#load data, fit a bunch of hierarchical models
rm(list=ls())
library(nlme)
library(MuMIn)
library(data.table)
library(boot)

#load data
d <- data.table(readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/tal_core_data.rds'))
#d <- data.table(readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/tal_all_prior_data.rds'))

d$relEM  <- (d$relEM  * (nrow(d) - 1) + 0.5) / nrow(d)
d$relAM  <- (d$relAM  * (nrow(d) - 1) + 0.5) / nrow(d)
d$relWR  <- (d$relWR  * (nrow(d) - 1) + 0.5) / nrow(d)
d$relSAP <- (d$relSAP * (nrow(d) - 1) + 0.5) / nrow(d)


#Fit ecto models
m <-  lm(logit(relEM) ~ doy + mat30 + map30 + pH + cn + horizon + wet.dep + latitude,  data = d)
m <- lme(logit(relEM) ~ doy + mat30 + map30 + pH + cn + horizon + wet.dep + latitude, random = ~1|site/plotID, data = d, na.action = na.omit)
m <- lme(logit(relEM) ~ n.dep + mat30 + moisture + horizon + latitude, random = ~1|site/plotID, data = d, na.action = na.omit)
summary(m)
r.squaredGLMM(m)
#plotting
test <- d[complete.cases(d[,.(relEM,n.dep,mat30,moisture,horizon,latitude)])]
plot(logit(test$relEM) ~ (fitted(m)), cex= 0.5, pch = 16)

#Fit an AM model
m <- lme(logit(relAM) ~ doy + mat30 + map30 + wet.dep + latitude + pH + cn + moisture + horizon, random = ~1|site/plotID, data = d, na.action = na.omit)
m <- lme(logit(relAM) ~ doy + mat30 + map30 + wet.dep + latitude + moisture , random = ~1|site/plotID, data = d, na.action = na.omit)
summary(m)
r.squaredGLMM(m)
test <- d[complete.cases(d[,.(doy,mat30,map30,wet.dep,latitude,moisture)])]
plot(logit(test$relAM) ~ fitted(m), cex=0.5, pch = 16)

#Fit a SAP model
m <- lme(logit(relSAP) ~ wet.dep + pH + moisture + horizon, random = ~1|site/plotID, data = d, na.action = na.omit)
summary(m)
r.squaredGLMM(m)
test <- d[complete.cases(d[,.(wet.dep,pH,moisture,horizon)])]
plot(logit(test$relSAP) ~ fitted(m), cex=0.5, pch = 16)


#Fit WR model
m <- lme(logit(relWR) ~ doy + mat30 + map30 + pH + cn + horizon + wet.dep, random = ~1|site/plotID, data = d, na.action = na.omit)
summary(m)
r.squaredGLMM(m)
test <- d[complete.cases(d[,.(doy,mat30,map30,pH,cn,horizon,wet.dep)])]
plot(logit(test$relWR) ~ fitted(m), cex=0.5, pch = 16)
