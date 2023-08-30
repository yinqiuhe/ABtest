#This is the R codes for testing the computation under a simplied multiple-mediator setting

#the R codes to generate original data of (S, eM, eY) when there are two mediators
#(without intercept )
gendata_r_mult_simplified <- function(n, num_med){
  mydata_tmp = matrix(NA, n, (2+num_med) );
  
  #1. exposure, binary
  mydata_tmp[,1] = rbinom(n, 1, 0.5); 
  
  #2. errors of M (mediator) models
  for( indj in 1:num_med){
    mydata_tmp[, (1+indj)] = rnorm(n, sd=0.5); 
  }
  
  #2. error Y model
  mydata_tmp[,(2+num_med)] = rnorm(n, sd=0.5);  
  
  return(mydata_tmp)
}

#function that generates additional covariates in model
gendata_r_covariates <- function(n, k){
  covariate_data = matrix(NA, n, k);
  
  covariate_data[ ,1] = rep(1,n);
  
  if(k>1){
    half_ind = floor((k-1)/2);
    # half binary and half continuous
    for (indk in 2:(half_ind+1)){
      covariate_data[,indk] = rbinom(n, 1, 0.5); 
    }
    for (indk in (half_ind+2):k){
      covariate_data[,indk] = rnorm(n, sd=0.5); 
    }
  }
  return(covariate_data)
}

summary_power1<- function(all_res_C, alpha){
  est_power = mean(all_res_C[1, ] < alpha) #adaptive MaxP
  est_power_non = mean(all_res_C[6 + N.boot.out,  ]< alpha ) #two-sided MaxP test with nonparametric bootstrap
  tmp_pval_asym = pnorm( all_res_C[5,] )
  pval_asym = 2 * apply( rbind( tmp_pval_asym, (1-tmp_pval_asym) ), 2, FUN=min )
  est_power_asym = mean( pval_asym  < alpha )
  all_power_res = c(est_power, est_power_non, est_power_asym)
  
  return(all_power_res)
}

summary_power2<- function(all_res_C, alpha){
  est_power = mean(all_res_C[1, ] < alpha) #adaptive MaxP
  est_power_non = mean(all_res_C[6 + N.boot.out,  ]< alpha ) #two-sided MaxP test with nonparametric bootstrap
  tmp_pval_asym = pnorm( all_res_C[ (7+N.boot.out) , ] )
  pval_asym = 2 * apply( rbind( tmp_pval_asym, (1-tmp_pval_asym) ), 2, FUN=min )
  est_power_asym = mean( pval_asym  < alpha )
  all_power_res = c(est_power, est_power_non, est_power_asym)
  
  return(all_power_res)
}

summary_power_m_out_n<-function( m_out_n_res, num_res, alpha){
  all_power_m_out_n_min = matrix(NA, 1, num_res*2)
  all_power_m_out_n_prod = matrix(NA, 1, num_res*2)
  ind = 1;
  for(ind1 in 1:2){
    for(ind2 in 1:num_res){
      all_power_m_out_n_min[ind] = mean( m_out_n_res[ind2, ind1, ] < alpha )
      ind = ind + 1 } }
  ind = 1;
  for(ind1 in 3:4){
    for(ind2 in 1:num_res){
      all_power_m_out_n_prod[ind] = mean( m_out_n_res[ind2, ind1, ] < alpha )
      ind = ind + 1} }
  
  return(rbind(all_power_m_out_n_min, all_power_m_out_n_prod))
}

plot_p_val<- function(all_res_C, figtext, all_power_res, sel_ind){
  clab=1.3
  cmain=1.5
  caxis=1.2
  pdf(file=figtext)
  est_power_legend = paste0("type I error/power: min ", all_power_res[1])
  est_power_legend_non = paste0("two side class boot: (blue) ", all_power_res[2])
  est_power_legend_asym = paste0("min asym (green) ", all_power_res[3])
  
  candidate_pval = c(5,  7+N.boot.out)
  tmp_pval_asym = pnorm( all_res_C[candidate_pval[sel_ind],])
  pval_asym = 2 * apply( rbind( tmp_pval_asym, (1-tmp_pval_asym) ), 2, FUN=min )
  
  candiate_name = c("MaxP test", "Product Coefficient Test")
  candiate_legend = c("JST", "Sobel Test" )
  
  plot(seq(0, 1, 0.01),quantile(all_res_C[1, ], probs = seq(0, 1, 0.01)), pch = 1, cex = 0.8, cex.lab=clab, cex.axis=caxis, cex.main=cmain, col="red", ylab = "Sample Quantiles", xlab = "Theoretical Quantiles", main = candiate_name[sel_ind], xlim=c(0,1), ylim=c(0,1))
  points(seq(0, 1, 0.01),quantile(all_res_C[(6 + N.boot.out),  ], probs = seq(0, 1, 0.01)), pch = 2, cex = 0.8, col="blue") #two-sided nonparametric bootstrap of MaxP
  points(seq(0,1,0.01), quantile(pval_asym, probs = seq(0, 1, 0.01)), col="green", pch = 0, cex = 0.8) 
  abline(0,1, col="orange")
  legend("topleft", legend = c( paste0("Adaptive-B: ", all_power_res[1]), paste0("Percentile-B:", all_power_res[2]), paste0(candiate_legend[sel_ind],": ",all_power_res[3])), pt.cex = c(1.3,1.3,1.3), cex = 1.3, pch=c(1, 2, 0), col=c("red","blue","green"), text.col = "black")
  
  #legend("topleft", legend = c("Adaptive-B", "Percentile-B", candiate_legend[sel_ind] ), pt.cex = c(1.3,1.3,1.3), cex = 1.4, pch=c(1, 2, 0), col=c("red","blue","green"), text.col = "black")
  
  #legend("topleft", legend = c(est_power_legend, est_power_legend_non, est_power_legend_asym))
  dev.off()
}




plot_p_val2<- function(all_res_C, figtext, all_power_res, sel_ind, MT_res, power_MT){
  clab=1.3
  cmain=1.5
  caxis=1.2
  pdf(file=figtext)
  est_power_legend = paste0("type I error/power: min ", all_power_res[1])
  est_power_legend_non = paste0("two side class boot: (blue) ", all_power_res[2])
  est_power_legend_asym = paste0("min asym (green) ", all_power_res[3])
  
  candidate_pval = c(5,  7+N.boot.out)
  tmp_pval_asym = pnorm( all_res_C[candidate_pval[sel_ind],])
  pval_asym = 2 * apply( rbind( tmp_pval_asym, (1-tmp_pval_asym) ), 2, FUN=min )
  
  candiate_name = c("MaxP test", "Product Coefficient Test")
  candiate_legend = c("JST", "Sobel Test" )
  
  all_colors=c("firebrick3","dodgerblue","forestgreen","darkorange")
  plot(seq(0, 1, 0.01),quantile(all_res_C[1, ], probs = seq(0, 1, 0.01)), pch = 1, cex = 0.8, cex.lab=clab, cex.axis=caxis, cex.main=cmain, font.main=2, col=all_colors[1], ylab = "Sample Quantiles", xlab = "Theoretical Quantiles", main = candiate_name[sel_ind], xlim=c(0,1), ylim=c(0,1))
  points(seq(0, 1, 0.01),quantile(all_res_C[(6 + N.boot.out),  ], probs = seq(0, 1, 0.01)), pch = 2, cex = 0.8, col=all_colors[2]) #two-sided nonparametric bootstrap of MaxP
  points(seq(0,1,0.01), quantile(pval_asym, probs = seq(0, 1, 0.01)), col=all_colors[3], pch = 0, cex = 0.8) 
  points( seq(0,1,0.01), quantile(MT_res, probs = seq(0, 1, 0.01)), col=all_colors[4], pch=4, cex = 0.8)
  abline(0,1, col="gold")
  legend("topleft", legend = c( paste0("Adaptive-B: ", all_power_res[1]), paste0("Percentile-B: ", all_power_res[2]), paste0(candiate_legend[sel_ind],": ",all_power_res[3]), paste0("MT-Comp: ",power_MT)), pt.cex = rep(1.4,4), pt.lwd=rep(1.4,4), cex = 1.4, pch=c(1, 2, 0, 4), col=all_colors, text.col = "black")
  
  #legend("topleft", legend = c("Adaptive-B", "Percentile-B", candiate_legend[sel_ind] ), pt.cex = c(1.3,1.3,1.3), cex = 1.4, pch=c(1, 2, 0), col=c("red","blue","green"), text.col = "black")
  
  #legend("topleft", legend = c(est_power_legend, est_power_legend_non, est_power_legend_asym))
  dev.off()
}





plot_p_val3<- function(all_res_C, figtext, all_power_res, sel_ind, MT_res, power_MT, package_res, power_package, m_out_n_res, power_m_out_n, power_num, N.boot.out){
  power_m_out_n_sel = power_m_out_n[sel_ind, ]
  num_power_res = length(power_m_out_n_sel) # 2 types of p-values times number of m choices
  clab=1.3
  cmain=1.5
  caxis=1.2
  pdf(file=figtext, width=8.5, height=6.3)
  # est_power_legend = paste0("type I error/power: min ", all_power_res[1])
  # est_power_legend_non = paste0("two side class boot: (blue) ", all_power_res[2])
  # est_power_legend_asym = paste0("min asym (green) ", all_power_res[3])
  # est_power_1 = paste0(  )
  
  par(xpd=TRUE, mar=c(5, 4, 4, 2)+c(0.2,11,0.2,0.2))
  
  candidate_pval = c(5,  7+N.boot.out) #select JST or product
  tmp_pval_asym = pnorm( all_res_C[candidate_pval[sel_ind],])
  pval_asym = 2 * apply( rbind( tmp_pval_asym, (1-tmp_pval_asym) ), 2, FUN=min )
  
  candiate_name = c("JST", "Product Coefficients Test")
  candiate_legend = c("JST", "Sobel" )
  
  # all_colors=c("firebrick3","dodgerblue","forestgreen","darkorange")
  all_colors=c("firebrick3","dodgerblue","forestgreen","darkorange", "red", "blue", "deeppink", "purple", "brown4", "gold", "snow4", "lightpink", "green")
  all_colors = all_colors[ 1:( 5 + num_power_res) ]
  plot(seq(0, 1, 0.01),quantile(all_res_C[1, ], probs = seq(0, 1, 0.01)), pch = 1, cex = 0.8, cex.lab=clab, cex.axis=caxis, cex.main=cmain, font.main=2, col=all_colors[1], ylab = "Sample Quantiles", xlab = "Theoretical Quantiles", main = candiate_name[sel_ind], xlim=c(0,1), ylim=c(0,1))
  points( seq(0, 1, 0.01),quantile(all_res_C[(6 + N.boot.out),  ], probs = seq(0, 1, 0.01)), pch = 2, cex = 0.8, col=all_colors[2] ) #two-sided nonparametric bootstrap of MaxP
  points( seq(0,1,0.01), quantile(pval_asym, probs = seq(0, 1, 0.01)), col=all_colors[3], pch = 0, cex = 0.8 ) 
  points( seq(0,1,0.01), quantile(MT_res, probs = seq(0, 1, 0.01)), col=all_colors[4], pch=4, cex = 0.8 )
  points( seq(0,1,0.01), quantile(package_res, probs = seq(0, 1, 0.01)), col=all_colors[5], pch=5, cex = 0.8 )
  
  color_ind = 6
  for(ind1 in (sel_ind*2-1):(sel_ind*2)){
    for(ind2 in 1:(num_power_res/2)){
      points( seq(0,1,0.01), quantile( m_out_n_res[ind2, ind1, ] , probs = seq(0, 1, 0.01)), col=all_colors[color_ind], pch=6, cex = 0.8 )
      color_ind = color_ind + 1
    }
  }
  all_names_m = NULL
  for( ind_res in 1:(num_power_res/2)){
    all_names_m = c(all_names_m,  paste0("PB m: ", round(power_num[ind_res],digits=2),": ", round(power_m_out_n_sel[ind_res],digits=3)) ) }
  for( ind_res in (num_power_res/2 + 1):(num_power_res) ){
    all_names_m = c(all_names_m,  paste0("CPB m: ",round(power_num[ind_res-num_power_res/2],digits=2),": ", round(power_m_out_n_sel[ind_res],digits=3)) ) }
  
  lines(c(-0.04,1.04),c(-0.04,1.04), col="black", lty = 2)
  
  legend(-0.65,0.8, legend = c( paste0("AB: ", round(all_power_res[1], digits=3)), 
                               paste0("PB: ", round(all_power_res[2], digits=3)), 
                               paste0(candiate_legend[sel_ind],": ",round(all_power_res[3], digits=3)), 
                               paste0("MT-Comp: ", round(power_MT, digits=3)),
                               paste0("package:", round(power_package, digits=3)), 
                               all_names_m ), pt.cex = 1, pt.lwd=1, cex = 1, pch=c(1, 2, 0, 4, 5, rep(6, num_power_res)), col=all_colors, text.col = all_colors)
  dev.off()
  #par(mar=c(5, 4, 4, 2) + 0.1)
  
  # legend(-0.8,0.8, legend = c( paste0("AB: ", all_power_res[1]), 
  #                               paste0("PB: ", all_power_res[2]), 
  #                               paste0(candiate_legend[sel_ind],": ",all_power_res[3]), 
  #                               paste0("MT-Comp: ",power_MT),
  #                               paste0("package:", power_package), 
  #                               all_names_m ), pt.cex = 1, pt.lwd=1, cex = 1, pch=c(1, 2, 0, 4), col=all_colors, text.col = all_colors)
  # 
  #legend("topleft", legend = c("Adaptive-B", "Percentile-B", candiate_legend[sel_ind] ), pt.cex = c(1.3,1.3,1.3), cex = 1.4, pch=c(1, 2, 0), col=c("red","blue","green"), text.col = "black")
  
  #legend("topleft", legend = c(est_power_legend, est_power_legend_non, est_power_legend_asym))

}




plot_lambdaval<- function(all_res_C, figtext_lambda_1, figtext_lambda_2, figtext_rej_prop){
  ##lambda 1 chosen
  pdf(file = figtext_lambda_1)
  hist(all_res_C[3, ], main = figtext_lambda_1, freq = FALSE)
  dev.off()
  print(summary(factor(all_res_C[3, ])))
  ##lambda 2 chosen
  pdf(file = figtext_lambda_2)
  hist(all_res_C[4, ], main =  figtext_lambda_2, freq = FALSE)
  dev.off()
  print(summary(factor(all_res_C[4, ])))
  ##rejection proportion
  
  pdf(file = figtext_rej_prop)
  hist( all_res_C[2, ] , main = figtext_rej_prop, freq = FALSE)
  dev.off()
  print(summary(all_res_C[2, ])) ##summary statistics
  print(summary(factor(all_res_C[2, ])))
}








#generate data for using mediation package
gene_data_second <- function(n, gamma_coef, beta_coef, other_gamma_coef, other_beta_coef, direct_effect){
  num_gamma = length(other_gamma_coef);
  num_beta = length(other_beta_coef); #the intercept, may be modified to a vector to incorporate other covariates later
  mydata_tmp = gendata_r(n, num_gamma, num_beta)
  
  my_data1 = matrix(NA, n, num_gamma+2)   #data for the model M ~ S
  my_data2 = matrix(NA, n, num_beta+3)    #data for the model Y ~ M
  
  S = mydata_tmp[,1]
  X = mydata_tmp[ , (2: (2+num_gamma-1))]
  errM = mydata_tmp[, 2+num_gamma ]
  my_data1[ , 1] = S                      #exposure: S
  my_data1[ , (2: (2+num_gamma-1)) ] = X  #confounders
  M = S * gamma_coef + X %*% other_gamma_coef + errM
  my_data1[ , 2+num_gamma] = M
  
  my_data2[ , 1] = M
  my_data2[ , (2: (2+num_beta-1)) ] = X  #same confounders as above
  errY = mydata_tmp[, 3+num_gamma+num_beta ]
  Y = M * beta_coef + X %*% other_beta_coef + direct_effect * S +  errY
  my_data2[ , (2+num_beta) ] = S
  my_data2[ , (3+num_beta) ] = Y
  
  return( list(my_data1, my_data2) )
}


test_mediate<-function(n, gamma_coef, beta_coef, other_gamma_coef, other_beta_coef, direct_effect){
  #generate data
  data_tmp = gene_data_second(n, gamma_coef, beta_coef, other_gamma_coef, other_beta_coef, direct_effect); 
  my_data1 = data_tmp[[1]];
  my_data2 = data_tmp[[2]];
  
  #fit linear regression model
  my_data1 = as.data.frame(my_data1)
  colnames(my_data1) = c("S", "X0", "X1", "X2", "M")
  my_data2 = as.data.frame(my_data2)
  colnames(my_data2) = c("M", "X0", "X1", "X2", "S", "Y")
  fit.sm <- lm( M ~ 0 + S + X0 + X1 + X2 , data = my_data1)
  fit.my <- lm( Y ~ 0 + M + X0 + X1 + X2 + S, data = my_data2)
  
  med_res <- summary(mediate(fit.sm, fit.my, treat ="S", mediator ="M", boot=TRUE))
  #return p-value
  return(med_res$d1.p)
}






# plot_p_val<- function(all_res_C, figtext, all_power_res, sel_ind){
#   clab=1.3
#   cmain=1.5
#   caxis=1.2
#   pdf(file=figtext)
#   est_power_legend = paste0("type I error/power: min ", all_power_res[1])
#   est_power_legend_non = paste0("two side class boot: (blue) ", all_power_res[2])
#   est_power_legend_asym = paste0("min asym (green) ", all_power_res[3])
#   
#   candidate_pval = c(5,  7+N.boot.out)
#   tmp_pval_asym = pnorm( all_res_C[candidate_pval[sel_ind],])
#   pval_asym = 2 * apply( rbind( tmp_pval_asym, (1-tmp_pval_asym) ), 2, FUN=min )
#   
#   candiate_name = c("MaxP test", "Product Coefficient Test")
#   
#   plot(seq(0, 1, 0.01),quantile(all_res_C[1, ], probs = seq(0, 1, 0.01)), pch = 20, cex.lab=clab, cex.axis=caxis, cex.main=cmain, col="red", ylab = "Sample Quantiles", xlab = "Theoretical Quantiles", main = candiate_name[sel_ind], xlim=c(0,1), ylim=c(0,1))
#   points(seq(0, 1, 0.01),quantile(all_res_C[(6 + N.boot.out),  ], probs = seq(0, 1, 0.01)), pch = 17, cex = 0.8, col="blue") #two-sided nonparametric bootstrap of MaxP
#   points(seq(0,1,0.01), quantile(pval_asym, probs = seq(0, 1, 0.01)), cex = 0.8, col="green", pch = 15) 
#   abline(0,1, col="orange")
#   legend("topleft", legend = c("Adaptive bootstrap", "Nonparametric bootsrap","Asymptotic JS"), pt.cex = c(1.8,1.3,1.3), cex = 1.4, pch=c(20,17,15), col=c("red","blue","green"), text.col = "black")
#   
#   #legend("topleft", legend = c(est_power_legend, est_power_legend_non, est_power_legend_asym))
#   dev.off()
# }









