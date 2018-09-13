#' forecast_plot.r
#' Takes a dataframe with observed data as `obs`, then predicted, 95%CI and 95%PI as columns, returns a plot that looks nice.
#'
#' @param to.plot   #dataframe with observations, mean prediction, uper and lower bounds of PI and CI for each prediction.
#' @param limy      #y limit values. Default 0,100
#' @param trans     #transparanecy of shaded region, 0-1. Default 0.3.
#' @param ylab      #y-axis label.
#'
#' @return
#' @export
#'
#' @examples
forecast_plot <- function(to.plot, limy=c(0,100),trans=0.3, ylab=NA){
  par(mai = c(0.1,0.8,0.1,0.1))
  plot(to.plot$mean, ylim = limy, cex = 0, ylab = NA, xlab = NA, xaxt='n',bty='n')
  lines(smooth.spline(to.plot$mean), lwd = 2) #plot mean line
  #plot and shade 95% credible interval.
  range <- c(1:nrow(to.plot))
  polygon(c(range, rev(range)),c(to.plot$pi_0.975, rev(to.plot$pi_0.025)), col=adjustcolor('green', trans), lty=0)
  polygon(c(range, rev(range)),c(to.plot$ci_0.975, rev(to.plot$ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
  #drop observed points.
  par(new = T)
  plot(to.plot$obs, cex = 0.6, pch =16, ylim=limy,xlab=NA,ylab=NA,xaxt='n',yaxt='n', bty='n')
  
  #get some summary stats on here.
  lpos <- 0.06
  pi_match <- round((nrow(to.plot[to.plot$obs < to.plot$pi_0.975 & to.plot$obs > to.plot$pi_0.025,]) / nrow(to.plot))*100,1)
  rsq <- round(summary(lm(to.plot$obs ~ to.plot$mean))$r.squared,2)
  txt <- paste0(pi_match,'% of obs. in 95% predictive interval.')
  mtext(txt,adj = lpos, line = -1.8)
  txt <- bquote(R^2 == .(format(rsq, digits = 2)))
  mtext(txt,adj = lpos, line = -3.1)
  mtext(ylab, side = 2, line = 2)
}