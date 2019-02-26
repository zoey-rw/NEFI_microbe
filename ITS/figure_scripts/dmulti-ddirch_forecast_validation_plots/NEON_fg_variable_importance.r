#plotting variable importance: functional groups.
#Clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
library(ggplot2)

#load data.
d <- readRDS(NEON_dmulti.ddirch_var.importance_fg.path)
mu <- d$mean
lo_95 <- d$ci_0.025
hi_95 <- d$ci_0.975

#drop "other" group.----
mu <- mu[,!(colnames(mu) %in% c('other'))]
lo_95 <- lo_95[,!(colnames(lo_95) %in% c('other'))] 
hi_95 <- hi_95[,!(colnames(hi_95) %in% c('other'))]

#Global plot settings.----
#i = 2 #2=Ectomycorrhizal fungi.
par(mfrow = c(2,2))
outer.lab.cex <- 1.5

#loop over different fungal groups.----
for(i in 2:ncol(d$mean)){
  #plot importance in decreasing order.
  mu   <- d$mean[,i][order(d$mean[,i], decreasing = T)]
  hi95 <- d$ci_0.975[,i][order(match(names(d$ci_0.975[,i]), names(mu)))]
  lo95 <- d$ci_0.025[,i][order(match(names(d$ci_0.025[,i]), names(mu)))]
  limy <- c(0,max(mu) * 1.01)
  lab <- names(mu)
  fungi.name <- colnames(d$mean)[i]
  J<-barplot(height = mu,
             beside = true, las = 2,
             ylim = limy,
             cex.names = 0.75, xaxt = "n",
             main = paste0(fungi.name,' fungi'),
             border = "black", axes = TRUE)
  #add error bars.
  #segments(J, lo95, J, hi95, lwd = 1.5)
  #add variable labels.
  ypos <- -limy[2] / 10
  text(x = J, y = ypos, srt = 45,
       adj = 1, labels = lab, xpd = TRUE)
}
#Outer labels.----
mtext('Variable Importance',side = 2, line = -1.5, cex = outer.lab.cex, outer = T)




