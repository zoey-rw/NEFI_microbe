#fitting as much as I can to tedersoo
#clear environment, source functions.
rm(list=ls())
source('NEFI_functions/crib_fun.r')
source('paths.r')

#load data
d <- data.table(readRDS(ted.ITSprior_data))


#transform [0,1] to (0,1) a la Cribari-Neto & Zeileis 2010
start <- grep("Ectomycorrhizal", colnames(d))
  end <- grep('hydrophobic', colnames(d))
d[,start:end] <- lapply(d[,start:end],crib_fun)

#functional group y data.
y_fg <- d[,.(Ectomycorrhizal,Saprotroph,Arbuscular,Pathogen,hydrophillic,hydrophobic)]
y_fg <- as.data.frame(y_fg)
#genus level y data.
start <- grep('Russula',colnames(d))
  end <- grep('Tricholoma',colnames(d))
y_ge <- d[,start:end]
y_ge <- as.data.frame(y_ge)

#Get matrix of predictors.
x <- d[,.(pH,moisture,pC,cn,doy,NPP,mat,map,mat_CV,map_CV,mdr,aridity,latitude,forest,conifer,relEM)]
x <- d[,.(pH,moisture,pC,cn    ,NPP,mat,map,mat_CV,map_CV,mdr,aridity,latitude,forest,conifer,relEM)]
#x$mat_relEM <- x$mat * x$relEM
x <- as.data.frame(x)

fg_models <- list()
for(i in 1:ncol(y_fg)){
  name <- colnames(y_fg)[i]
  data <- cbind(y_fg[,i],x)
  colnames(data)[1] <- 'y'
  data <- data[complete.cases(data),]
  preds <- paste(colnames(x), collapse = ' + ')
  form <- as.formula(paste0('y ~ ',preds))
  mod <- betareg::betareg(form, data = data)
   pred <- fitted(mod)
    obs <- data$y
  resid <- obs - pred
  super.list <- list(mod,pred,obs,resid)
  names(super.list) <- c('model','predicted','observed','residual')
  fg_models[[i]] <- super.list
}
names(fg_models) <- colnames(y_fg)

ge_models <- list()
for(i in 1:ncol(y_ge)){
  name <- colnames(y_ge)[i]
  data <- cbind(y_ge[,i],x)
  colnames(data)[1] <- 'y'
  data <- data[complete.cases(data),]
  preds <- paste(colnames(x), collapse = ' + ')
  form <- as.formula(paste0('y ~ ',preds))
  mod <- betareg::betareg(form, data = data)
  pred <- fitted(mod)
  obs <- data$y
  resid <- obs - pred
  super.list <- list(mod,pred,obs,resid)
  names(super.list) <- c('model','predicted','observed','residual')
  ge_models[[i]] <- super.list
}
names(ge_models) <- colnames(y_ge)

#png(filename='figures/ted_ITSprior_fg_fits.png',width=12,height=10,units='in',res=300)
par(mfrow = c(2,3))
for(i in 1:length(fg_models)){
  name <- names(fg_models)[i]
  rsq <- round(summary(fg_models[[i]]$model)$pseudo.r.squared,3)
  plot(fg_models[[i]]$observed ~ fg_models[[i]]$predicted, cex = 2) #on relative scale
  abline(lm(fg_models[[i]]$observed ~ fg_models[[i]]$predicted), lty = 2, col = 'green')
  #plot(boot::logit(fg_models[[i]]$observed) ~ boot::logit(fg_models[[i]]$predicted), cex = 0.8) #logit transform to linear.
  #abline(lm(boot::logit(fg_models[[i]]$observed) ~ boot::logit(fg_models[[i]]$predicted)), lty = 2, col = 'green')
  abline(0,1,lwd=2)
  mtext(name, line = -2, adj = 0.05, col = 'purple')
  mtext(paste0('R2 = ',rsq), line = -3.2, adj = 0.05, col = 'purple')
}
#dev.off()

#png(filename='figures/ted_ITSprior_20.top.genera_fits.png',width=12,height=10,units='in',res=300)
par(mfrow = c(5,4),
    mai = c(0.3,0.3,0.3,0.3))
for(i in 1:length(ge_models)){
  name <- names(ge_models)[i]
  rsq <- round(summary(ge_models[[i]]$model)$pseudo.r.squared,3)
  plot(ge_models[[i]]$observed ~ ge_models[[i]]$predicted, cex = 0.8)
  abline(lm(ge_models[[i]]$observed ~ ge_models[[i]]$predicted), lty = 2, col = 'green')
  #plot(boot::logit(ge_models[[i]]$observed) ~ boot::logit(ge_models[[i]]$predicted), cex = 0.8) #logit transform to linear.
  #abline(lm(boot::logit(ge_models[[i]]$observed) ~ boot::logit(ge_models[[i]]$predicted)), lty = 2, col = 'green')
  abline(0,1,lwd=2)
  mtext(name, line = -2, adj = 0.05, col = 'purple')
  mtext(paste0('R2 = ',rsq), line = -3.2, adj = 0.05, col = 'purple')
}
#dev.off()
