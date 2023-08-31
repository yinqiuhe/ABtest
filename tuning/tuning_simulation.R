#This file presents codes of conducting double bootstrap in simulated data
#Methods are described in Section "Double Bootstrap for Choosing the Tuning Parameter" of Supplementary Material


################################################
#Load packages and functions
################################################
#add simulations with fixed lambda
library(Rcpp)
library(RcppArmadillo)
olibrary(DBmypackage6two)
library(qqplotr) #for plotting results
library(BoutrosLab.plotting.general) #for plotting results
source('Rsource_fun1.R') 
source('qqplot_with_band.R')

################################################
#Set parameters for simulating data
################################################
n = 200 #sample size
omega = 0.05 #significance level
N_boot_out = 500; #outer bootstrap size
N_boot_in = 500; #inner bootstrap size
direct_effect = 1 #direct effect tau
other_coef = 1 #coefficients of intercept 
coef_X_alpha = other_coef #mediator model coefficients
coef_X_beta = other_coef #outcome model coefficients
num_other_alpha = length(coef_X_alpha) #number of adjusted covariates (including intercept) in the mediator model
num_other_beta = length(coef_X_beta) + 1 #number of adjusted covariates (including intercept) in the outcome model #add 1 because of the direct effect S
alpha_coef = 0.5; #alphaS value
beta_coef = 0;  #betaM value
mc_ind = 123
set.seed(mc_ind)


######set tuning parameter######
lambda1all = 0; lambda2all = 0; #value 0 indicates that AB reduces to classical bootstrap
lambda_n1_alpha = lambda1all*sqrt(n)/log(n); lambda_n2_alpha = lambda2all*sqrt(n)/log(n)
lambda_n1_beta = lambda1all*sqrt(n)/log(n); lambda_n2_beta = lambda2all*sqrt(n)/log(n)
lambda_n1_both = lambda1all*sqrt(n)/log(n); lambda_n2_both = lambda2all*sqrt(n)/log(n)
################################################


################################################
######1. Simulate data######
################################################
num_med = 1 #number of mediators
mydata_tmp = gendata_r_mult_simplified(n, num_med); 
other_covariates = gendata_r_covariates(n, length(other_coef));
mydata_original = DBmypackage6two::gendata_model_X_C( n,  alpha_coef, beta_coef,  mydata_tmp, other_covariates,
                                                      coef_X_alpha, coef_X_beta, direct_effect, num_med)

################################################
######2. Run double bootstrap ######
################################################
res_tmp = DBmypackage6two::boot_tuning_par_two_C(mydata_original, other_covariates, 
                                                 n,  omega,  num_med,
                                                 lambda_n1_alpha, lambda_n2_alpha, 
                                                 lambda_n1_beta, lambda_n2_beta, 
                                                 lambda_n1_both, lambda_n2_both, 
                                                 N_boot_out, N_boot_in,
                                                 num_other_alpha, num_other_beta)

################################################
######3. Plot p-values obtained by double bootstrap ######
################################################
omega1=0.05 #plot p-values of 1-omega1 level confidence band

#3.1 p-values formed by processed data with \hat{\alpha}=0  (P_{alpha}*(lambda))
alpha_p_val_boot = res_tmp$Boot_res_alpha_0
qqplot_with_band(alpha_p_val_boot, omega1)

#3.2 p-values formed by processed data with \hat{\beta}=0  (P_{beta}*(lambda)) 
beta_p_val_boot = res_tmp$Boot_res_beta_0 
qqplot_with_band(beta_p_val_boot, omega1)

#3.3 p-values formed by processed data with  \hat{\alpha}=0 and \hat{\beta}=0  (P_{beta}*(lambda))  P_{alpha,beta}*(lambda)
both_p_val_boot = res_tmp$Boot_res_both_0
qqplot_with_band(both_p_val_boot, omega1)

################################################
######4. Interpretation of results
################################################
#When generating data with alpha_coef = 0.5 and beta_coef=0 and set lambda=0, we expect
#3.1 processed data with \hat{\alpha}=0 is very conservative 
#3.2 processed data with \hat{\beta}=0 is close to uniform distribution (as alpha_coef is non-zero)
#3.3 processed data with \hat{\alpha}=0 and \hat{\beta}=0 is very conservative  




