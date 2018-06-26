#load data, fit a bunch of hierarchical models
rm(list=ls())
library(nlme)
library(MuMIn)
library(data.table)
library(boot)

#load data
d <- data.table(readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/tal_all_prior_data.rds'))

start.position <- 20 #column number where relative abundances start

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
fixed.effects <- d[,.(doy,MAT,MAP,MAT_CV,MAP_CV,pH,cn,Horizon,latitude)]
random.effects <- d[,.(Site,Plot)]

#Run a model loop.
for(i in 1:ncol(to_model)){
  y <- to_model[,to_model[,..i]] #grab dependent variable
  data <- cbind(y,fixed.effects,random.effects) #build dataframe for analysis
  #specify model formula
  formula <- as.formula(paste(paste(colnames(y)),'~', paste(colnames(fixed.effects), collapse="+")))
  #run model
  mod <- lme(formula, random = ~1|Site/Plot, data = data, na.action = na.omit)
  sum <- data.frame(summary(mod)$tTable)
  sum <- sum[-1,]
  #remove non-significant predictors.
  new.fixed <- rownames(sum[sum$p.value < 0.05,])
  #fix categorical variable names
  new.fixed[new.fixed == 'HorizonOH'] <- 'Horizon'
  #specify a new formula of fixed effects, if any predictors are significant
  if(length(new.fixed) > 0){
    formula <- as.formula(paste(paste(colnames(y)),'~', paste(new.fixed, collapse="+")))
    mod <- lme(formula, random = ~1|Site/Plot, data = data, na.action = na.omit)        
  }
  #get model summary, R2 values, and a complete dataset with fitted values.
  out.data <- data[,c(y,data[,..new.fixed],random.effects)]
  out.data <- out.data[complete.cases(out.data),]
  out.data$fitted <- fitted(mod)
  r.sq <- r.squaredGLMM(mod)
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

d <- model.out.list$hydrophobic$data
plot(hydrophobic ~ fitted, data = d, cex = 0.5, pch = 16)
