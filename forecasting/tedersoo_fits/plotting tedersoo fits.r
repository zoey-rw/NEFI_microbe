#plot tedersoo fits
rm(list=ls())
library(data.table)
source('NEFI_functions/crib_fun.r')
d <- readRDS('/fs/data3/caverill/NEFI_microbial/model_fits/ted_frequentist.rds')


png(filename='figures/tedersoo_fits.png',width=10,height=10,units='in',res=300)

par(mfrow = c(5,5),
    mai = c(0,0,0,0),
    oma = c(5,5,1,1))

for(i in 1:length(d)){
  z <- d[[i]]
  plot(z$data[,1] ~ z$data$fitted, ylab = NA,xlab=NA, cex = 0.5, pch = 16)
  abline(lm(z$data[,1] ~ z$data$fitted), lwd =2)
  mtext(paste('r2=',round(z$r.sq,2)), side = 3, line = -5, adj = 0.05, cex = 0.8, col = 'purple')
  mtext(paste(names(d)[i])          , side = 3, line = -3, adj = 0.05, cex = 0.8, col = 'purple')
}

mtext('observed', side = 2, line = 3, outer = T)
mtext('fitted'  , side = 1, line = 3, outer = T)

dev.off()