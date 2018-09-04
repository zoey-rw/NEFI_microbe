#plottingSITE level forecasts to NEON.
#clear environment, source paths, functions and packages.
rm(list=ls())
source('paths.r')
source('NEFI_functions/forecast_plot.r')

#output path.
figure.path <- paste0('figures/NEON_site.level_20gen_forecast.png')

#load and format forecast.----
d <- readRDS(NEON_site_fcast_20gen.path)

of_interest <- d$all.preds
of_interest <- lapply(of_interest,as.data.frame)
#get percentages by multiplying by 100.
for(i in 1:length(of_interest)){
  of_interest[[i]] <- of_interest[[i]] * 100
}


#set some plot features.----
png(filename=figure.path,width=14,height=13,units='in',res=300)
limy <- c(0,100)
trans <- 0.3
par(mfrow=c(4,5),
    mai = c(.8,.8,.1,.1))

#Actual plot code.----
#loop over groups.
for(i in 2:ncol(of_interest$mean)){
  grp.name <- colnames(of_interest$mean)[i]
  lab <- rownames(of_interest$mean)
  plot(of_interest$mean[,i], ylim = limy, xlab = NA, xaxt='n', bty = 'n', ylab = NA, cex=0)
  axis(1, at=1:nrow(of_interest$mean), labels=lab, las = 2)
  lines(smooth.spline(of_interest$mean[,i], spar = 0), lwd = 2)
  range <- c(1:nrow(of_interest$mean))
  polygon(c(range, rev(range)),c(of_interest$pi_0.975[,i], rev(of_interest$pi_0.025[,i])), col=adjustcolor('green', trans), lty=0)
  polygon(c(range, rev(range)),c(of_interest$ci_0.975[,i], rev(of_interest$ci_0.025[,i])), col=adjustcolor('blue' , trans), lty=0)
  label <- paste0('Predicted % ',grp.name,' Fungi')
  mtext(label, side = 2, line = 2.4, cex = 1.1)
}

#End plot.----
dev.off()
