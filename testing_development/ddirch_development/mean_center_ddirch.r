#May 25, 2018
#Fitting a dirichlet to 3 species abundances, predicted by an intercept and two covariates.
#This solution presented by Martyn Plummer.
#Precdictors are mean centered for fitting, and this is backed out whne reporting parameter estimates.
#clear environment, load packages.
rm(list = ls())
library(runjags)

#simulate data.
n <- 100
y1 <- round(rnorm(n, 100, 5)) 
mat <- rnorm(length(y1), 10, 3)
map <- rnorm(length(y1),100,10)
y1 <- y1 + round(mat*-3) + round(map*1)
y2 <- round(rnorm(n, 100, 5))
y3 <- round(rnorm(n, 100, 5))
spp.y <- as.matrix((data.frame(y1,y2,y3)))
y <- spp.y / rowSums(spp.y)
x <- as.matrix(cbind(rep(1,nrow(y)),mat,map)) #x predictors, including intercept

test <- matrix(ncol = ncol(x), nrow = nrow(x))
for(i in 1:nrow(x)){
  for(k in 1:ncol(x)){
    test[i,k] <- x[i,k] - mean(x[,k])
  }
}

#Specify JAGS model.
dirichlet.model = "
model {
#setup priors for each species
alpha ~ dnorm(0, 1.0E-3) 
for(i in 1:N.x){
beta[i,1] <- 0
for (j in 2:N.spp) {
beta[i,j] ~ dnorm(0, 1.0E-3)
}
}

#mean center predictors.
for(i in 1:N){
  x_center[i,1] <- 1
  for(j in 2:N.x){
    x_center[i,j] <- x[i,j] - mean(x[,j])
    }
}

#save centers for downstream transform.
for(j in 1:N.x){
  center_save[j] <- mean(x[,j])
}

#implement dirlichet
for (i in 1:N){
  for (j in 1:N.spp){
    #log(a0[i,j]) <- alpha + inprod(beta[,j], x[i,])
    log(a0[i,j]) <- alpha + inprod(beta[,j], x_center[i,])
    }
  y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp]) 
}

#map to original parameterization
for (j in 1:N.spp) {
#m[1,j] <- alpha + beta[1,j]
m[1,j] <- alpha + beta[1,j] - inprod(beta[2:N.x,j], center_save[2:N.x])
for (i in 2:N.x){
m[i,j] <- beta[i,j]
}
}

} #close model loop.
"

#Setup JAGS data object, fit JAGS model.
jags.data <- list(y = y, x = x,
                  N = nrow(y), N.spp = ncol(y), N.x = ncol(x))
jags.out <- runjags::run.jags(dirichlet.model,
                              data=jags.data,
                              adapt = 200,
                              burnin = 200,
                              sample = 200,
                              n.chains=3,
                              monitor=c('m'))
#summarize output.
out <- summary(jags.out)


###visualize.
#extract mean parameter matrix.
parms <- (matrix(out[grep('m', rownames(out)),4],nrow = (ncol(x)),ncol = (ncol(y))))
#multiply parameters by predictors to get predicted values.
pred.list <- list()
for(i in 1:ncol(y)){
  pred.list[[i]] <- exp(x %*% parms[,i])
}
pred.list <- (as.matrix(do.call('cbind', pred.list)))
pred.list <- pred.list / rowSums(pred.list)

#Plot predicted vs. observed and the 1:1 line.
par(mfrow = c(1,3))
for(i in 1:ncol(y)){
  plot(y[,i] ~ pred.list[,i])
  mod <- lm(y[,i] ~ pred.list[,i])
  r.sq <- summary(mod)$r.squared
  abline(0,1, lwd = 2) 
  abline(mod, lwd = 2, lty = 2, col = 'green')
  txt <- paste0('R2 = ',round(r.sq,2))
  mtext(txt, side = 3, line = -2, adj = 0.05)
}