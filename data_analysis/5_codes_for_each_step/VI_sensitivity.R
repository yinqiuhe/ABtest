################################################################
#Step 0: Preparation
################################################################
library(mediation)
library(abind)
library(parallel)
source('pre_functions_mult.R')  

#read data
exposure_data <- readRDS("exposure_data.Rds") 
outcome_data <- readRDS("outcome_data.Rds")
mediator_data <- readRDS("mediator_data.Rds")
other_covariates <- readRDS("other_covariates.Rds")


################################################################
#Step 1: conduct sensitivity analysis
################################################################
Num_mediator_sens <- (dim(mediator_data)[2]-1)

N_cores = 8
all_res_sensitivity <- mclapply(1:Num_mediator, function(j){
  set.seed(j); 
  #### process data (include only j-th mediator)
  data_model <- data_combine_function( j,  exposure_data, outcome_data, mediator_data, other_covariates)
  #n_all=nrow(data_model)   #sample size
  med.fit <- lm(mediator ~ exposure + age + SEX, data = data_model)
  out.fit <- lm(outcome ~ mediator + exposure + age + SEX, data = data_model)
  med.out <- mediate(med.fit, out.fit, treat = "exposure", mediator = "mediator", robustSE = TRUE, sims = 100)
  # med.out_list[[j]] = med.out
  
  sens.out <- medsens(med.out, rho.by = 0.1, effect.type = "indirect", sims = 100)
  # sens.out_list[[j]]=sens.out
  return( list("ind" = j, "med" = med.out, "sens" = sens.out) )
}, 
mc.cores = N_cores
)

# saveRDS(all_res_sensitivity, file = "realdata_sens2.Rds")


################################################################
#Step 2: processing outputs
################################################################

res_sens <- all_res_sensitivity

all_med_names <- colnames(mediator_data)[-1]
all_sel_ind = c(3, 40, 41, 53, 132, 133, 136, 2, 14, 50, 52, 67, 74, 78, 95)
sel_med_names <- all_med_names[all_sel_ind]
Num_mediator= Num_mediator_sens

all_var_corr <- matrix(NA,Num_mediator,4)
alpha_val <- matrix(NA, Num_mediator, 2) #the first column is alpha.hat, and the second column is its p-value
beta_val <- matrix(NA, Num_mediator, 2) #the first column is beta.hat, and the second column is its p-value
med_val <- rep(NA, Num_mediator) #alpha.hat * beta.hat
#residual correlation by each mediator
for(j in 1:Num_mediator){
  data_model <- data_combine_function( j,  exposure_data, outcome_data, mediator_data, other_covariates)
  med.fit <- lm(mediator ~ exposure + age + SEX + 1, data = data_model)
  out.fit <- lm(outcome ~ mediator + exposure + age + SEX + 1, data = data_model)
  
  corval <- cor(med.fit$residuals, out.fit$residuals)
  var_out <- var(out.fit$residuals)
  var_med <- var(med.fit$residuals)
  
  all_var_corr[j,1] = corval
  all_var_corr[j,2] = var_out
  all_var_corr[j,3] = var_med
  
  med.fit.summary <- summary(med.fit)
  out.fit.summary <- summary(out.fit)
  alpha_val[j, ] <- med.fit.summary$coefficients[2, c(1,4)]
  beta_val[j, ] <- out.fit.summary$coefficients[2, c(1,4)]
  med_val[j] = alpha_val[j, 1] * beta_val[j, 1]
}

min_rho_zero <- rep(NA, Num_mediator)
for(j in 1:Num_mediator){
  # med.out <- res_sens[[j]]$med
  sens.out <- res_sens[[j]]$sens
  #keep the minimum rho that has confounded results
  min_rho_zero[j] = sens.out$err.cr.d
}

sens_res <- data.frame(sel_med_names, med_val[all_sel_ind],  alpha_val[all_sel_ind,], beta_val[all_sel_ind,], min_rho_zero[all_sel_ind],  all_var_corr[all_sel_ind,1])

greeks = c(alpha='\u03b1', tau='\u03c4', sigma='\u03c3', sigmaSq='\u03c3\u00B2', beta='\u03b2', gamma='\u03b3')

colnames(sens_res) = c("", "NIE", "alpha", "p_a", "beta", "p_b", "rho", "rho_m")

sens_res[,1] = c("X16.0.LYSO.PC_3", "FA 7:0-OH._1", "FA 7:0-OH._2", "FA 18:3_1", 
                 "GLY", "GLY-H2O", "LAURIC.ACID", "X16.0.LYSO.PC_2",  "FA.16:0-OH", 
                 "LPC 16:1_3",  "LPC 18:2", "URSODIOL_2", "DECENOYL", "FA 12:0-DiC", "FA 18:0-DiC")

sens_res
