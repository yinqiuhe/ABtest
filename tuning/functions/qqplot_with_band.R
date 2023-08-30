library(qqplotr)

qqplot_with_band<-function(plot_boot_sample, omega){
  library(BoutrosLab.plotting.general)
  ks_alpha=ks.test.critical.value( length(plot_boot_sample) , 1-omega, alternative = "two.sided");
  all_prob = seq(0, 1, length.out=N_boot_out+1)
  plot(qunif(all_prob), 
       quantile(plot_boot_sample, probs = all_prob), 
       pch = 1, cex = 0.8, cex.main=1, 
       col="red", ylab = "Sample Quantiles", xlab = "Theoretical Quantiles", 
       main = "qq-plot-boot-p-values", xlim=c(0,1), ylim=c(0,1))
  abline(0,1)
  e_cdf_sample = ecdf(plot_boot_sample)
  y1=e_cdf_sample(all_prob)+ks_alpha
  y2=e_cdf_sample(all_prob)-ks_alpha
  
  y1poly=y1 * (y1<=1) + 1*(y1>1)
  y2poly=y2 * (y2>=0)
  library(scales)
  mycolor=alpha("blue", alpha = 0.1)
  polygon(c(y2poly, rev(y1poly)), c(all_prob, rev(all_prob)), 
          col = mycolor, lty = 0)
  n_d_1=which(y1<=1)
  n_d_2=which(y2>=0)
  lines(y1[n_d_1], all_prob[n_d_1], col="blue", lty=2)      
  lines(y2[n_d_2], all_prob[n_d_2], col="blue", lty=2) 
}







qqplot_with_band_v2<-function(plot_boot_sample, omega, figtext){
  library(BoutrosLab.plotting.general)
  pdf(file = figtext, width=5.5, height=5.5) #for adding box in legend
  par(mar=c(4.5, 5, 1.4, 0.5))
  ks_alpha=ks.test.critical.value( length(plot_boot_sample) , 1-omega, alternative = "two.sided");
  N_boot_out = length(plot_boot_sample)
  all_prob = seq(0, 1, length.out=N_boot_out+1)
  plot(qunif(all_prob), 
       quantile(plot_boot_sample, probs = all_prob), 
       pch = 1, cex = 1, cex.main=1, 
       col="red", ylab = "Sample Quantiles", xlab = "Theoretical Quantiles", 
       cex.lab = 2, cex.axis = 2 , 
       # main = "qq-plot-boot-p-values", 
       xlim=c(0,1), ylim=c(0,1))
  abline(0,1)
  e_cdf_sample = ecdf(plot_boot_sample)
  y1=e_cdf_sample(all_prob)+ks_alpha
  y2=e_cdf_sample(all_prob)-ks_alpha
  
  y1poly=y1 * (y1<=1) + 1*(y1>1)
  y2poly=y2 * (y2>=0)
  library(scales)
  mycolor=alpha("blue", alpha = 0.1)
  polygon(c(y2poly, rev(y1poly)), c(all_prob, rev(all_prob)), 
          col = mycolor, lty = 0)
  n_d_1=which(y1<=1)
  n_d_2=which(y2>=0)
  lines(y1[n_d_1], all_prob[n_d_1], col="blue", lty=2)      
  lines(y2[n_d_2], all_prob[n_d_2], col="blue", lty=2) 
  dev.off()
}




