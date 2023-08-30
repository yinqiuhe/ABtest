################################################################
#Step 0: Preparation
################################################################
### (i) Load packages
library(dplyr)
library(Rcpp)
library(RcppArmadillo)
library(mediation)
library(DBmypackage6two)
library(abind)
library(parallel)

# setwd("~/Library/CloudStorage/Dropbox/work/work_13_mediation/Mediation Test Data/revision summary")
#functions for running the test in Huang 2019
source("MT_Comp.R")
Rcpp::sourceCpp('DBC7_MT_Comp.cpp')
source('pre_functions_mult.R')  

### (ii) Read data
exposure_data <- readRDS("exposure_data.Rds")
outcome_data <- readRDS("outcome_data.Rds")
mediator_data <- readRDS("mediator_data.Rds")
other_covariates <- readRDS("other_covariates.Rds")


### (iii) Set Parameters
N.boot.out = 10^5
N.boot.in = 0
# n=nrow(exposure_data)-3#delete three subjects with missing SEX
num_alpha = 1 + 2 #number of baseline covariates: intercept, SEX, age
num_beta = 1 + 2  #number of baseline covariates: intercept, SEX, age
lambda1all = 2 #c(1.5, 2)
lambda2all = lambda1all
Num_mediator= 149
omega=0.05
N_cores=36


################################################################
#Step 1: Obtain p-values of Joint Significance Tests
################################################################
#(adaptive, classic, asymptotic)
all_res_max <- mclapply(1:Num_mediator, function(j){
  set.seed(j); 
  #### process data (include only j-th mediator)
  data_model <- data_combine_function( j,  exposure_data, outcome_data, mediator_data, other_covariates)
  n_all=nrow(data_model)   #sample size
  # mydata <- pre_project_data(data_model) #project data for inference  
  pvals_0 <- p_val_max(n_all, data_model, N.boot.out,  N.boot.in, omega, lambda1all, lambda2all)
  
  return( pvals_0 )
}, 
mc.cores = N_cores
)
# saveRDS(all_res_max, file = 'realdata_max.Rds')


################################################################
#Step 2: Obtain p-values of Product of Coefficients Tests
################################################################
#(adaptive, classic, asymptotic)
all_res_prod <- mclapply(1:Num_mediator, function(j){
  set.seed(j); 
  #### process data (include only j-th mediator)
  data_model <- data_combine_function( j,  exposure_data, outcome_data, mediator_data, other_covariates)
  n_all=nrow(data_model)   #sample size
  # mydata <- pre_project_data(data_model) #project data for inference  
  pvals_0 <- p_val_prod(n_all, data_model, N.boot.out,  N.boot.in, omega, lambda1all, lambda2all)
  
  return( pvals_0 )
}, 
mc.cores = N_cores
)
# saveRDS(all_res_prod, file = 'realdata_prod.Rds')


##################################################################################
#Step 3: Obtain p-values of other methods (MT-Huang and Causal Mediation Package)
##################################################################################
all_res_other <- mclapply(1:Num_mediator, function(j){
  set.seed(j); 
  #### process data (include only j-th mediator)
  data_model <- data_combine_function( j,  exposure_data, outcome_data, mediator_data, other_covariates)
  n_all=nrow(data_model)   #sample size
  mydata <- pre_project_data(data_model) #project data for inference  
  
  other_cova1 <- data_model[,c(3,4,5)]
  mydata_original_tmp <- data_model[, c(1, 6, 2)] #read data
  
  sim_num <- 10^3
  pvals_0 <- p_val_other(n_all, mydata, mydata_original_tmp, other_cova1, sim_num, num_alpha, num_beta, N.boot.out,  N.boot.in, omega, lambda1all, lambda2all)
  
  return( pvals_0 )
}, 
mc.cores = N_cores
)
# saveRDS(all_res_other, file = 'realdata_other.Rds')



##################################################################################
#Step 4: Processing output results
##################################################################################
num_test = 8
res_all=matrix(NA, Num_mediator, num_test)

#Estimated p-valuse for 149 mediators and 8 tests
for (j in 1:Num_mediator){
  #results of maximum statistic
  res_all[j,1:3] = all_res_max[[j]]
  
  #results of product statistic
  res_all[j,4:6] = all_res_prod[[j]]
  
  #results of other tests
  res_all[j,7:8] = all_res_other[[j]]
}

colnametext = c("J_AB","J_B", "J_Asym", "P_AB","P_B", "P_Sobel", "MT", "CMA")
colnames(res_all) = colnametext
res_all




