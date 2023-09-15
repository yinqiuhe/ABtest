#This file presents codes of conducting double bootstrap in real data
#Methods are described in Section "Double Bootstrap for Choosing the Tuning Parameter" of Supplementary Material
#Due to privacy contraint, we illustrate the procedure with the noisy data (raw data with added noise)

####################################################################################
#1. read data and set up
########################################################
exposure_data <- readRDS("dataset_noisy/exposure_data.Rds")
outcome_data <- readRDS("dataset_noisy/outcome_data.Rds")
mediator_data <- readRDS("dataset_noisy/mediator_data.Rds")
other_covariates <- readRDS("dataset_noisy/other_covariates.Rds")
other_covariates <- other_covariates[,c("FOLIOCC","Intercept")]

library(DBmypackage6two)
source('functions/pre_functions_mult.R')
source('functions/qqplot_with_band.R')
library(mediation)

all_med = colnames(mediator_data)
omega = 0.05
N_boot_out = 400;
N_boot_in = 400; 
Num_mediator= 149

n = 382
num_med = 1; #number of tested mediators
mc_ind = 1234
set.seed(mc_ind);

#indices of selected mediators
#sel_ind = c(3, 40, 41, 53, 132, 133, 136, 2, 14, 50, 52, 67, 74, 78, 95)
#names of all selected mediators 
#all_med[sel_ind]

j=136 #LAURIC.ACID #a rejected mediator
# j=10 #FA.12.0-OH (F.12) #a non-rejected mediator
data_model <- data_combine_function0( j,  exposure_data, outcome_data, mediator_data, other_covariates)
n = nrow(data_model)
mydata_original = cbind(data_model$exposure, data_model$mediator, data_model$outcome)
other_covariates = cbind(rep(1, n), data_model$SEX, data_model$age)

num_other_alpha = 3
num_other_beta = 3 + 1 #additional 1 because of the direct effect S
########################################################

########################################################
#2. process data and apply double bootstrap 
########################################################
#set tuning parameter
lambda1all <- lambda2all <- 0 * sqrt(n)/log(n); 
lambda_n1_alpha <- lambda_n2_alpha <- lambda_n1_beta <- lambda_n2_beta <- lambda_n1_both <- lambda_n2_both <- lambda_n1_original<- lambda_n2_original <- lambda1all

#run double bootstrap
res_tmp = DBmypackage6two::boot_tuning_par_all_C(mydata_original, other_covariates, 
                                                 n,  omega,  num_med,
                                                 lambda_n1_alpha, lambda_n2_alpha, 
                                                 lambda_n1_beta, lambda_n2_beta, 
                                                 lambda_n1_both, lambda_n2_both, 
                                                 lambda_n1_original, lambda_n2_original, 
                                                 N_boot_out, N_boot_in,
                                                 num_other_alpha, num_other_beta)

########################################################
#3. plot results
########################################################
omega1=0.05 #plot p-values of 1-omega1 level confidence band

#3.0 p-values estimated by applying double bootstrap to observed data
plot_boot_sample = res_tmp$Boot_res_original
qqplot_with_band_separate(plot_boot_sample, omega1 )

#3.1 p-values estimated by applying double bootstrap to processed data with \hat{\alpha}=0  (P_{alpha}*(lambda))
alpha_p_val_boot = res_tmp$Boot_res_alpha_0
qqplot_with_band(alpha_p_val_boot, omega1)

#3.2 p-values estimated by applying double bootstrap to  processed data with \hat{\beta}=0  (P_{beta}*(lambda)) 
beta_p_val_boot = res_tmp$Boot_res_beta_0 
qqplot_with_band(beta_p_val_boot, omega1)

#3.3 p-values estimated by applying double bootstrap to  processed data with  \hat{\alpha}=0 and \hat{\beta}=0  (P_{beta}*(lambda))  P_{alpha,beta}*(lambda)
both_p_val_boot = res_tmp$Boot_res_both_0
qqplot_with_band(both_p_val_boot, omega1)


################################################
######4. Interpretation of results
################################################
#For rejected mediators, it is likely that alpha_coef is nonzero and beta_coef is nonzero,
#setting lambda=0, we expect
#3.0 double bootstrap to original data gives qq-plots of p-values bend downward if under alternative hypothesis, and vice versa
#3.1 processed data with \hat{\alpha}=0 is  close to uniform distribution (as beta_coef is non-zero)
#3.2 processed data with \hat{\beta}=0 is close to uniform distribution (as alpha_coef is non-zero)
#3.3 processed data with \hat{\alpha}=0 and \hat{\beta}=0 is very conservative  




