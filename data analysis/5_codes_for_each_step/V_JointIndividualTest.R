################################################################
#Step 0: Preparation
################################################################
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


################################################################
#Step 1: Read Data
################################################################
exposure_data <- readRDS("exposure_data.Rds")
outcome_data <- readRDS("outcome_data.Rds")
mediator_data <- readRDS("mediator_data.Rds")
other_covariates <- readRDS("other_covariates.Rds")

all_sel_ind = c(3, 40, 41, 53, 132, 133, 136, 2, 14, 50, 52, 67, 74, 78, 95)
#link exposure, outcome, X, and selected mediators
data_model_tmp <- merge(exposure_data, outcome_data, by="FOLIOCC")  
data_model <- merge(data_model_tmp, other_covariates, by="FOLIOCC") 
mediator_data <- mediator_data[ ,c(1,all_sel_ind+1)]  #selected mediators, +1 as index of mediators start from 2
data_model <- merge(data_model, mediator_data, by="FOLIOCC") 
colnames(data_model)[1:6] <- c("subject", "exposure", "outcome", "SEX", "age", "Intercept")
data_model = na.omit(data_model) #remove NA values



################################################################
#Step 2: Conduct Tests
################################################################
n = nrow(data_model)
omega <- 0.05
N.boot.out <- 10^4 #bootstrap number
num_med <- length(all_sel_ind) #number of selected mediators

individual_test <- function(ind_j) {
  set.seed(1234)
  
  ind_choose <- ind_j + 6
  other_covariates <- data_model[,-c(1, 2, 3, ind_choose)]
  mydata_original_tmp <- data_model[, c(2, ind_choose, 3)] #read data
  num_med_tested <- 1
  
  ############################################################
  #Product type methods
  ############################################################
  lambda1all <- lambda2all <- 2 #c(1.5, 2)
  N.boot.in <- 0
  res_Prod <- DBmypackage6two::One_MC_DB_X_noProj_C(as.matrix(mydata_original_tmp), 
                                                    as.matrix(other_covariates), 
                                                    n, omega, num_med_tested, lambda1all, 
                                                    lambda2all, N.boot.out, N.boot.in)
  #results from asymptotic normality
  tmp_pval_asym_prod = pnorm(res_Prod[(N.boot.out+7)])
  pval_asym_prod = 2 * min(tmp_pval_asym_prod, (1- tmp_pval_asym_prod) ) #two-sided p-value
  
  ############################################################
  #Maximum type methods
  ############################################################
  lambda1all = 2 
  lambda2all = lambda1all
  res_maxP <- DBmypackage6two::One_MC_DB_max_X_noProj_C(as.matrix(mydata_original_tmp), 
                                                        as.matrix(other_covariates), 
                                                        n, omega, num_med_tested, lambda1all, 
                                                        lambda2all, N.boot.out, N.boot.in)
  
  #results from asymptotic normality
  tmp_pval_asym_maxP = pnorm(res_maxP[5])
  pval_asym_maxP = 2 * min( tmp_pval_asym_maxP, (1-tmp_pval_asym_maxP) )
  
  #####################################
  #compute p-values by other tests
  #####################################
  
  ################################################################################
  #Comparison 1:  Huang's method: 
  #compute projected data: mydata_proj
  colnames(mydata_original_tmp)[2] <- "mediator"
  mydata_proj <- pre_project_data_multi(as.matrix(mydata_original_tmp), 
                                        as.matrix(other_covariates))
  num_alpha = 1 + 2 + num_med - 1 # intercept, SEX, age, and other mediators
  num_beta = 1 + 2  + num_med - 1 # intercept, SEX, age, and other mediators
  res_MT <- My_MT_Comp(mydata_proj, n, num_alpha, num_beta) 
  
  AB = computeTest_C(mydata_proj, n, num_alpha, num_beta)
  mediation_effect<- AB[1]*AB[2]
  
  ################################################################################
  #Comparison 2: the mediation package method
  sim_num <- 10^3 #change to a larger number for reproducing the results
  res_package <- test_mediate(n, mydata_original_tmp, other_covariates, sim_num)
  
  return(c(res_Prod[1], res_Prod[6], pval_asym_prod, 
           res_maxP[1], res_maxP[6], pval_asym_maxP, 
           res_MT, res_package, mediation_effect)) }





num_cores=num_med
all_boot_pval <- do.call(abind, list( mclapply(mc.cores = num_cores, 1:num_med, function(ind_j){individual_test(ind_j)}), along = 2))



################################################################
#Step 3: Processing and Save Outputs
################################################################
all_boot_pval <- t(all_boot_pval)

sel_med_names <- colnames(mediator_data)[-1]
rownames(all_boot_pval) <- sel_med_names

all_mediation_effect <- all_boot_pval[,9]
all_boot_pval <- all_boot_pval[, 1:8]
colnames(all_boot_pval) <- c("PoC-AB", "PoC-B", "PoC-Sobel", "JS-AB", "JS-B", "JS-MaxP", "MT", "CMA")


all_boot_pval[,c("JS-AB", "JS-MaxP", "PoC-AB", "PoC-B", "PoC-Sobel",  "CMA")]


# saveRDS(all_boot_pval, file = "all_boot_pval.Rds")










