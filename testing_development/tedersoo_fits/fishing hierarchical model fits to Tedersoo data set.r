#load data, fit a bunch of hierarchical models
rm(list=ls())
library(nlme)
library(MuMIn)
library(data.table)
library(boot)
library(betareg)

#output path.
output.path <- '/fs/data3/caverill/NEFI_microbial/model_fits/ted_frequentist.rds'

#load data
d <- data.table(readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/ted_all_prior_data.rds'))

start.position <- grep("Ectomycorrhizal", colnames(d))

#function to transform [0,1] to (0,1) a la Cribari-Neto & Zeileis 2010
crib.fun <- function(x){(x * (length(x) - 1) + 0.5) / length(x)}
d[,start.position:ncol(d)] <- lapply(d[,start.position:ncol(d)],crib.fun)
#then logit transform.
d[,start.position:ncol(d)] <- lapply(d[,start.position:ncol(d)],logit)

#fit models of relative abundance
model.out.list <- list()
#specify columns to model- all relative abundances
to_model <- (d[,start.position:ncol(d)])
#specify fixed and random effects data.
fixed.effects <- d[,.(doy,MAT,MAP,MAT_CV,MAP_CV,pH,cn,latitude)]
random.effects <- d[,.(Site)]

#Run a model loop.
for(i in 1:ncol(to_model)){
  y <- to_model[,to_model[,..i]] #grab dependent variable
  data <- cbind(y,fixed.effects,random.effects) #build dataframe for analysis
  #specify model formula
  formula <- as.formula(paste(paste(colnames(y)),'~', paste(colnames(fixed.effects), collapse="+")))
  #run model
  mod <-  lm(formula, data = data)
  sum <- data.frame(summary(mod)[[4]])
  sum <- sum[-1,]
  #perform stepAIC on model, subset fixed effects to only include these predictors.
  s.AIC <- data.frame(summary(stepAIC(mod))[[4]])[-1,]
  fixed.effects.AIC <- fixed.effects[,rownames(s.AIC), with = F]
  #specify a new formula of fixed effects, if any predictors are significant
  if(length(fixed.effects.AIC) > 0){
    formula <- as.formula(paste(paste(colnames(y)),'~', paste(colnames(fixed.effects.AIC), collapse="+")))
    mod <- lm(formula, data = data)        
  }
  #get model summary, R2 values, and a complete dataset with fitted values.
  out.data <- data.frame(y,fixed.effects.AIC,random.effects)
  out.data <- out.data[complete.cases(out.data),]
  out.data$fitted <- fitted(mod)
  #r.sq <- r.squaredGLMM(mod)
  r.sq <- summary(mod)$r.squared
  report.list <- list(mod,r.sq,out.data)
  names(report.list) <- c('model','r.sq','data')
  model.out.list[[i]] <- report.list
  names(model.out.list)[[i]] <- colnames(y)
}

#sort by r2 values because you are awful.
z <- data.frame(names(model.out.list))
out.r2m <- list()
out.r2c <- list()
out.preds <- list()
for(k in 1:length(model.out.list)){
  r2m <- model.out.list[[k]]$r.sq[1]
  r2c <- model.out.list[[k]]$r.sq[2]
  preds <- rownames(summary(model.out.list[[k]]$model)$tTable)
  out.r2c[[k]] <- r2c
  out.r2m[[k]] <- r2m
  out.preds[[k]] <- preds
}
z$r2m <- do.call('rbind',out.r2m)
z$r2c <- do.call('rbind',out.r2c)
z[order(-z$r2m),]

d <- model.out.list$Oidiodendron$data
plot(Oidiodendron ~ fitted, data = d, cex = 0.5, pch = 16)

saveRDS(model.out.list,output.path)
