#This file tests the joint mediation effect of 15 selected lipids
################################################################
#Step 0: preparation
################################################################
library(DBmypackage6two)


################################################################
#Step 1: read data
################################################################
exposure_data <- readRDS("exposure_data.Rds")
outcome_data <- readRDS("outcome_data.Rds")
mediator_data <- readRDS("mediator_data.Rds")
other_covariates <- readRDS("other_covariates.Rds")

#link exposure, outcome, X, and selected mediators
data_model_tmp <- merge(exposure_data, outcome_data, by="FOLIOCC")  
data_model <- merge(data_model_tmp, other_covariates, by="FOLIOCC") 

all_sel_ind = c(3, 40, 41, 53, 132, 133, 136, 2, 14, 50, 52, 67, 74, 78, 95) #indexes of selected lipids from Step 1 screeing
mediator_data <- mediator_data[ ,c(1,all_sel_ind+1)]  #selected mediators, +1 as index of mediators start from 2
data_model <- merge(data_model, mediator_data, by="FOLIOCC") 

colnames(data_model)[1:6] <- c("subject", "exposure", "outcome", "SEX", "age", "Intercept")
data_model = na.omit(data_model) #remove NA values

num_med <- length(all_sel_ind) #number of selected mediators
mydata_original = data_model[, c(2, 7:(6+num_med), 3)] #exposure, selected mediators, outcome

other_covariates = data_model[,c("SEX", "age", "Intercept")]
mydata_original = data_model[, c(2, 7:(6+num_med), 3)] #read data
n = nrow(mydata_original)
mydata_original = as.matrix(mydata_original)
other_covariates = as.matrix(other_covariates)


################################################################
#Step 2: Conduct the Global Test
################################################################
N.boot.out <- 10^4 #10^4

set.seed(1234)
lambda1all <- 2; lambda2all <- lambda1all
omega=0.05
res_Prod <- DBmypackage6two::One_MC_DB_X_noProj_C(mydata_original, other_covariates, 
                                                  n, omega, num_med, lambda1all, 
                                                  lambda2all, N.boot.out, 0)

res_Prod[1] 


################################################################
#Step 3: Comparison with the global test in Wei and Song (2022)
################################################################

source("LR.R")
#read data to compare with other methods
k = ncol(mydata_original)
X_gen = mydata_original[, 1] #exposure
M_gen = mydata_original[, 2:(k-1)] #mediators
Y_gen = as.matrix(mydata_original[, k]) #outcome
Z_gen = other_covariates
#compare by the codes in Wei
pvals_other_tmp=cluster_pval(X=X_gen,M=M_gen,Y=Y_gen,Z=Z_gen,total_iter=1000,tol=1e-6,total_iter_lambda=10,m=10000)
pvals_other = c(pvals_other_tmp$p_val_LR, pvals_other_tmp$p_val_PTN,  pvals_other_tmp$p_val_PTNP )
pvals_other


