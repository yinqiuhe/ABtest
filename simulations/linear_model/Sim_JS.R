
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(abind)
library(DBmypackage6two)

source('Rsource_fun1.R') 

#path of the results
#all_res_C_1 #results when lambda=2

n = 500
omega = 0.05
N_boot_out = 500; #bootstrap number
N_boot_in = 0;
num_cores = 8 #number of cores for parallel computing
N_rep = 500 #the number of Monte Carlo replication
R_seed = 1234
set.seed(R_seed)
all_seeds = floor(runif(N_rep) * 10^5)

#coefficients in model
other_coef = c(1,1,1) #only intercept if 1, (intercept + binary + continuous)
direct_effect = 1
coef_X_alpha = other_coef
coef_X_beta = other_coef
num_other_alpha = length(coef_X_alpha)
num_other_beta = length(coef_X_beta) + 1 #additional 1 because of the direct effect S
num_med = 1;  #set number of mediators
alpha_coef = 0; beta_coef = 0; 
lambda1all = 2; lambda2all = lambda1all;

########################################################
#test the maximum statistic
all_res_C = do.call(abind, list( mclapply(mc.cores = num_cores, all_seeds, function(mc_ind) {
  set.seed(mc_ind);

  #generate original data of (S, eM1, ..., eMm, eY)
  mydata_tmp = gendata_r_mult_simplified(n, num_med);

  #generate additional covariates including the intercept
  other_covariates = gendata_r_covariates(n, length(other_coef));

  #compute data from the model S, M1, ..., Mm, Y
  mydata_original = DBmypackage6two::gendata_model_X_C( n,  alpha_coef, beta_coef,  mydata_tmp, other_covariates,
                                                        coef_X_alpha, coef_X_beta, direct_effect, num_med)

  DBmypackage6two::One_MC_DB_max_X_noProj_C(mydata_original, other_covariates, n, omega,
                                            num_med, lambda1all, lambda2all, N_boot_out, N_boot_in)} ),
  along = 2))

########################################################
main_title_a = "alpha: ("; main_title_b = "beta: ("
for (inj in 1:length(alpha_coef)){  
  main_title_a = paste(main_title_a,alpha_coef[inj]);
  main_title_b = paste(main_title_b, beta_coef[inj])
}
main_title = paste(main_title_a, ") \n",  main_title_b, ")" ); cmain=2
plot(seq(0, 1, 0.01),quantile(all_res_C[1, ], probs = seq(0, 1, 0.01)), 
     pch = 1, cex = 0.8, cex.main=cmain, #cex.lab=clab, cex.axis=caxis,  
     col="red", ylab = "Sample Quantiles", xlab = "Theoretical Quantiles", 
     main = main_title, xlim=c(0,1), ylim=c(0,1))
points(seq(0, 1, 0.01),quantile(all_res_C[6,  ], probs = seq(0, 1, 0.01)), pch = 2, cex = 0.8, col="blue")
tmp_pval_asym_prod = pnorm( all_res_C[5,  ] ) #maximum statistic
pval_asym = 2 * apply(cbind(tmp_pval_asym_prod, (1- tmp_pval_asym_prod) ), 1 ,min)
points(seq(0, 1, 0.01), quantile(pval_asym, probs = seq(0, 1, 0.01)), pch = 1, cex = 0.8, col="green")
abline( 0, 1, col="orange")

all_power = c(mean(all_res_C[1, ]<omega), mean(all_res_C[6, ]<omega), mean(pval_asym<omega))
names(all_power) <-c("AB-P", "B-P", "Sobel")
all_power
saveRDS(all_res_C, file = paste0(fig_folder, "all_res.Rds"))

paste(  mean(all_res_C[1,] < 0.05), "&",mean(all_res_C[6,] < 0.05) )
paste(  mean(all_res_C[1,] < 0.1), "&",mean(all_res_C[6,] < 0.1) )






