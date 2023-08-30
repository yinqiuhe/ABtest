#######################################################
#the preprossing projection function
pre_project_data_multi <-function(mydata_original_tmp, other_covariates){ 
  n = nrow(mydata_original_tmp) #sample size
  #for outcome ~ mediator model, remove the exposure (direct effect), and intercept, SEX and age and other mediators
  mediator_tilde_O <- lm( mydata_original_tmp[,"mediator"] ~ 0 + mydata_original_tmp[,"exposure"] + other_covariates )$residuals
  outcome_tilde <- lm( mydata_original_tmp[,"outcome"] ~ 0 + mydata_original_tmp[,"exposure"] + other_covariates)$residuals
  
  #for mediator ~ exposure model, remove the intercept, and SEX and age and other mediators
  exposure_tilde <- lm( mydata_original_tmp[,"exposure"] ~ 0 + other_covariates )$residuals
  mediator_tilde_E <- lm( mydata_original_tmp[,"mediator"] ~ 0 + other_covariates )$residuals
  
  my_data_model <- as.matrix(cbind(outcome_tilde, exposure_tilde, mediator_tilde_O, mediator_tilde_E))
  # colnames(my_data_model) <- c("outcome","exposure","mediator_O","mediator_E")
  
  return(my_data_model)
}


#######################################################
#the preprossing projection function
pre_project_data <-function(data_model){
  n = nrow(data_model) #sample size
  #for outcome ~ mediator model, remove the intercept and exposure (direct effect), and SEX and age
  mediator_tilde_O <- lm( data_model$mediator ~ 0 + data_model$Intercept +  data_model$exposure + data_model$SEX + data_model$age )$residuals
  outcome_tilde <- lm( data_model$outcome ~ 0 + data_model$Intercept +  data_model$exposure + data_model$SEX + data_model$age )$residuals

  #for mediator ~ exposure model, remove the intercept, and SEX and age
  exposure_tilde <- lm( data_model$exposure ~ 0 + data_model$Intercept + data_model$SEX + data_model$age )$residuals
  mediator_tilde_E <- lm(data_model$mediator~ 0 + data_model$Intercept + data_model$SEX + data_model$age )$residuals

  my_data_model <- as.matrix(cbind(outcome_tilde, exposure_tilde, mediator_tilde_O, mediator_tilde_E))
  # colnames(my_data_model) <- c("outcome","exposure","mediator_O","mediator_E")

  return(my_data_model)
}

#######################################################
#compute max statistic p-value function 
p_val_max <- function(n, data_model, N.boot.out,  N.boot.in, omega, lambda1all, lambda2all){
  #data_model: is the original data
  #mydata: is the original data data_model after projection
  # res_maxP <- One_MC_DoubleBoot_C(mydata, n, num_gamma, num_beta, alpha, lambda1all, lambda2all, N.boot.out, N.boot.in) 
  
  mydata_original <- as.matrix(data_model[,c("exposure", "mediator", "outcome")]) #original data without projection
  other_covariates <- as.matrix(data_model[, c("Intercept", "SEX", "age")])
  n = nrow(mydata_original)
  num_med <- 1
  res_maxP <- DBmypackage5two::One_MC_DB_max_X_noProj_C(mydata_original, other_covariates, n, omega, 
                                                 num_med, lambda1all, lambda2all, N.boot.out, N.boot.in)
  
  #results from asymptotic normality
  tmp_pval_asym_maxP = pnorm(res_maxP[5])
  pval_asym_maxP = 2 * min( tmp_pval_asym_maxP, (1-tmp_pval_asym_maxP) )
  
  all_pval <- c(res_maxP[1],res_maxP[6], pval_asym_maxP)
  names(all_pval) <- c("maxP","maxP_class", "maxP_Asym")
  
  return(all_pval)
}

#######################################################
#compute product statistic p-value function 
p_val_prod <- function(n, data_model, N.boot.out,  N.boot.in, omega, lambda1all, lambda2all){
  #data_model: is the original data
  #mydata: is the original data data_model after projection
  # res_Prod <- One_MC_DoubleBoot_C(mydata, n, num_gamma, num_beta, alpha, lambda1all, lambda2all, N.boot.out, N.boot.in) 
  
  mydata_original <- as.matrix(data_model[,c("exposure", "mediator", "outcome")]) #original data without projection
  other_covariates <- as.matrix(data_model[, c("Intercept", "SEX", "age")])
  n = nrow(mydata_original)
  num_med <- 1
  res_Prod <- DBmypackage5two::One_MC_DB_X_noProj_C(mydata_original, other_covariates, n, omega, 
                                                 num_med, lambda1all, lambda2all, N.boot.out, N.boot.in)
  
  #results from asymptotic normality
  tmp_pval_asym_prod = pnorm(res_Prod[(N.boot.out+7)])
  pval_asym_prod = 2 * min(tmp_pval_asym_prod, (1- tmp_pval_asym_prod) ) #two-sided p-value
  
  all_pval <- c(res_Prod[1], res_Prod[6], pval_asym_prod)
  names(all_pval) <- c("Prod","Prod_class", "Prod_Asym")
  
  return(all_pval)
}

#####################################
#compute p-value for other statistics
p_val_other <- function(n, mydata, data_original, other_covariates, sim_num, num_gamma, num_beta, N.boot.out,  N.boot.in, omega, lambda1all, lambda2all){
  
  ################################################################################
  #Compare 1:  Huang's method: 
  res_MT <- My_MT_Comp(mydata, n, num_gamma, num_beta) 
  
  ################################################################################
  #Compare 2: the mediation package method
  res_package <- test_mediate(n, data_original, other_covariates, sim_num)
  
  all_pval <- c(res_MT, res_package)
  names(all_pval) <- c("MT", "Package")
  
  return(all_pval)   
}


#####################################
#data combine
data_combine_function <- function( j,  exposure_data, outcome_data, mediator_data, other_covariates){
  data_model <- merge(exposure_data, outcome_data, by="FOLIOCC")
  data_model <- merge(data_model, other_covariates, by="FOLIOCC")
  mediator_data <- mediator_data[, c(1, j+1)]
  data_model <- merge(data_model, mediator_data, by="FOLIOCC")
  colnames(data_model) <- c("subject", "exposure", "outcome", "SEX", "age", "Intercept", "mediator")
  data_model <- data_model[,-1]  #delete "subject label" 
  missing_ind <- which(is.na((data_model$SEX))) #indexes of data missing SEX information
  data_model <- data_model[-missing_ind, ] #delete data with missing SEX

  return(data_model)
}



#####################################
#data combine
data_combine_function0 <- function( j,  exposure_data, outcome_data, mediator_data, other_covariates){
  data_model <- merge(exposure_data, outcome_data, by="FOLIOCC")
  data_model <- merge(data_model, other_covariates, by="FOLIOCC")
  mediator_data <- mediator_data[, c(1, j+1)]
  data_model <- merge(data_model, mediator_data, by="FOLIOCC")
  colnames(data_model) <- c("subject", "exposure", "outcome", "Intercept", "mediator")
  data_model <- data_model[,-1]  #delete "subject label" 
  data_model <- na.omit(data_model)
  return(data_model)
}
