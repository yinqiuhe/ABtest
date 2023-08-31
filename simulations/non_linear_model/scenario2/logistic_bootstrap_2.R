source('compute_stat2.R', echo=TRUE)

library(parallel)
library(abind)

#bootstrap of logistic regression

#data generation
#0. parameters
n <- 500
alpha_coef <- 0
beta_coef <- 0

#repeat simulation
B_num <- 500
n_repeat <- 1000
num_cores <- 8 #number of cores for parallel computing
set.seed(312)
all_seeds = floor(runif(n_repeat) * 10^5)

########################################
#classical and adaptive bootstrap
########################################
lambda1 <- 1.9; lambda2 <- 3.3
lambda_alpha <- lambda1*sqrt(n)/log(n); lambda_beta <- lambda2*sqrt(n)/log(n)

all_boot_pval <- do.call(abind, list( mclapply(mc.cores = num_cores, all_seeds, function(seed.num) {
  set.seed(seed.num)
  
  my_data <- data_generate_model(n, alpha_coef, beta_coef)
  tmp_res <- compute_stat(my_data)
  
  test_stat <- tmp_res$me
  alpha_residual <- tmp_res$alpha_residual
  beta_residual <- tmp_res$beta_residual
  z_alpha <- tmp_res$z_alpha
  z_beta <- tmp_res$z_beta
  g_alpha_fit_values <- tmp_res$g_alpha_fit_values

  one_boot_res <- one_ab_bootstrap(my_data, B_num, test_stat, 
                                   alpha_residual, beta_residual,
                                   z_alpha, z_beta, 
                                   lambda_alpha, lambda_beta,
                                   g_alpha_fit_values)
  #Result 1: classical bootstrap
  tmp_p_class <- mean( one_boot_res[1, ] > 0) 
  p_val_classical_boot <- 2 * min(tmp_p_class, 1-tmp_p_class)
  
  #Result 2: adaptive bootstrap
  # tmp_p_ab <- mean( one_boot_res[2, ] > test_stat) 
  tmp_p_ab <- mean( one_boot_res[2, ] > 0) 
  p_val_ab_boot <- 2 * min(tmp_p_ab, 1-tmp_p_ab)
  
  return( c(p_val_classical_boot, p_val_ab_boot)) }), along = 2))

##################################################
#plot qq-plots
namestxt=c("AB","B")
ptcex=1.1
clab=2.1
cmain=2.2
caxis=2

cex_legend = 1.96
ptlwd_legend = 1.96
ptcex_legend = 1.96

figtext = paste0("qq_plot.pdf")
pdf(file=figtext,width=7, height=7)
par( mar = c(5, 4.6, 0.3, 0.2) + 0.1 )
plot(seq(0, 1, 0.01),quantile(all_boot_pval[1, ], probs = seq(0, 1, 0.01)), 
     pch = 1, cex = 0.8, cex.lab=clab, cex.axis=caxis, 
     ylab = "Sample Quantiles", xlab = "Theoretical Quantiles", type = "n",
     xlim=c(0,1), ylim=c(0,1))
points(seq(0, 1, 0.01), quantile(all_boot_pval[1,  ], probs = seq(0, 1, 0.01)), pch = 2, col="blue", lwd=1.2, cex = ptcex )
points(seq(0, 1, 0.01), quantile(all_boot_pval[2,  ], probs = seq(0, 1, 0.01)), pch = 1, col="red", lwd=1.2, cex = ptcex )
abline( 0, 1, col="orange")
legend("bottomright", namestxt, col = c("red", "blue"), pch = c(1,2), 
       cex = cex_legend, pt.cex = ptcex_legend, pt.lwd=ptlwd_legend)
dev.off()



saveRDS(all_boot_pval, file = paste0("logit_boot_pval.Rds"))

##################################################
hist(all_boot_pval[1,])
hist(all_boot_pval[2,])
#rejection rates
paste(  mean(all_boot_pval[1,] < 0.05), "&",mean(all_boot_pval[2,] < 0.05) )
paste(  mean(all_boot_pval[1,] < 0.1), "&",mean(all_boot_pval[2,] < 0.1) )

rej_rates = rbind(cbind(mean(all_boot_pval[1,] < 0.05), mean(all_boot_pval[2,] < 0.05)), 
                  cbind(mean(all_boot_pval[1,] < 0.1), mean(all_boot_pval[2,] < 0.1)))
saveRDS(rej_rates, file = paste0("logit_power.Rds"))


ratio_coef
alpha_coef
beta_coef



