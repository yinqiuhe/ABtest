################################################################
#logistic functions
################################################################
logit_transform_fun <- function(x){
  l_val <- exp(x)/(1+exp(x))
  return(l_val)
}
logit_deriv <- function(x){
  l_val <- exp(x)/(1+exp(x))^2
  return(l_val)
}
logit_inverse_fun <- function(x){
  return( log(x/(1-x)) )
}

################################################################
#fit alpha logistic regression, beta linear regression, and ME
################################################################
compute_stat <- function(boot_data){
  ########################################
  #1. fit logistic: alpha coefficient
  ########################################
  ########################################
  #equivalent closed-form
  # S_boot <- boot_data$exposure
  # M_boot <- boot_data$mediator
  lm_alpha_tmp <- glm(mediator ~ exposure + covariate1, family=binomial(link='logit'),data=boot_data)
  lm_alpha <- summary(lm_alpha_tmp)
  result_alpha <- lm_alpha$coefficients
  alpha_hat <- result_alpha["exposure", "Estimate"]
  alpha_int_hat <- result_alpha["(Intercept)", "Estimate"]
  z_alpha <- result_alpha["exposure", "z value"]
  
  g_alpha_hat <- logit_transform_fun( alpha_hat + alpha_int_hat)
  g_alpha_hat_0 <- logit_transform_fun( alpha_int_hat)
  d_g_alpha_hat <- g_alpha_hat - g_alpha_hat_0
  g_alpha_fit_values <- lm_alpha_tmp$fitted.values
  
  alpha_residual <- boot_data$mediator - lm_alpha_tmp$fitted.values
  
  ########################################
  #2. fit linear: beta coefficient
  ########################################
  lm_beta <- summary(lm( outcome ~ mediator + exposure + covariate1, data=boot_data))
  result_beta <- lm_beta$coefficients
  beta_hat <- result_beta["mediator", "Estimate"]
  z_beta <- result_beta["mediator", "t value"]
  beta_residual <- lm_beta$residuals 
  
  #get projected M after bootstrap
  M_proj_boot <- lm(mediator ~ exposure + covariate1, data=boot_data)$residuals
  
  ########################################
  #3. compute test statistic
  ########################################
  n <- nrow(boot_data)
  mediation_estimate <- sqrt(n) * as.numeric(beta_hat * d_g_alpha_hat)
  
  return(list(me=mediation_estimate, alpha_residual = alpha_residual, 
              beta_residual = beta_residual, z_alpha = z_alpha, z_beta = z_beta,
              M_proj_boot = M_proj_boot, g_alpha_fit_values = g_alpha_fit_values,
              alpha_hat = alpha_hat, alpha_int_hat = alpha_int_hat))
}

################################################################
#compute classical bootstrap
################################################################
one_bootstrap <- function(my_data, B_num){
  n <- nrow(my_data)
  one_boot_res <- sapply( 1:B_num, function(b, my_data, n){
    boot_index <- sample(n, replace=T)
    boot_data <- my_data[boot_index,]
    tmp_estiamte <- compute_stat(boot_data)
    mediation_estiamte <- tmp_estiamte$me
    return(mediation_estiamte)
  }, my_data = my_data, n = n )
  
  return(one_boot_res)
}

################################################################
#compute the local bootstrap at (0,0)
################################################################
compute_local_stat <- function(boot_data, M_proj_boot, alpha_residual, beta_residual,
                               alpha_hat, alpha_int_hat, g_alpha_fit_values){

  W_alpha_t <- c( logit_deriv(alpha_hat + alpha_int_hat),  0, 0 )
  
  D_boot_mat <- as.matrix(boot_data[, c("exposure", "covariate0", "covariate1")])
  V_alpha <- t(D_boot_mat) %*% sweep( D_boot_mat, MARGIN = 1, g_alpha_fit_values * (1-g_alpha_fit_values), '*' )
  
  alpha_stat_res <- W_alpha_t %*% solve(V_alpha) %*% t(D_boot_mat) %*% alpha_residual
  beta_stat_res <- sum( M_proj_boot * beta_residual ) / sum(M_proj_boot^2)
  
  n <- length(alpha_residual)
  residual_stat <- sqrt(n) * alpha_stat_res * beta_stat_res
  return(residual_stat)
}

################################################################
#compute adaptive bootstrap
################################################################
one_ab_bootstrap <- function(my_data, B_num, test_stat, 
                             alpha_residual, beta_residual,
                             z_alpha, z_beta, 
                             lambda_alpha, lambda_beta,
                             g_alpha_fit_values){
  n <- nrow(my_data)
  one_boot_res <- sapply( 1:B_num, function(b, my_data, n){
    #1. bootstrap indexes and then data
    boot_index <- sample(n, replace=T)
    boot_data <- my_data[boot_index,]
    boot_res_alpha <- alpha_residual[boot_index]
    boot_res_beta <- beta_residual[boot_index]
    boot_g_alpha_fit_values <- g_alpha_fit_values[boot_index]
    
    #2. compute classical bootstrap statistic
    tmp_boot_res <- compute_stat(boot_data)
    mediation_est_boot <- tmp_boot_res$me #classical bootstrap estimate
    z_alpha_boot <- tmp_boot_res$z_alpha
    z_beta_boot <- tmp_boot_res$z_beta
    M_proj_boot <- tmp_boot_res$M_proj_boot
    # g_alpha_fit_values <- tmp_boot_res$g_alpha_fit_values
    alpha_hat_boot <- tmp_boot_res$alpha_hat
    alpha_int_hat_boot <- tmp_boot_res$alpha_int_hat
    
    #3. depends on passing threshold or not
    t_alpha <- ( abs(z_alpha) <= lambda_alpha ) & ( abs(z_alpha_boot) <= lambda_alpha)
    t_beta <- ( abs(z_beta) <= lambda_beta ) & ( abs(z_beta_boot) <= lambda_beta )
    
    #4. compute local expansion bootstrap
    if( t_alpha & t_beta ){
      S_boot <- boot_data$exposure; M_boot <- boot_data$mediator
      ab_stat <- compute_local_stat(boot_data, M_proj_boot, boot_res_alpha, boot_res_beta,
                                    alpha_hat_boot, alpha_int_hat_boot, 
                                    boot_g_alpha_fit_values) - test_stat
    }else{
      
      ab_stat <- mediation_est_boot 
    }
    
    return( c(mediation_est_boot, ab_stat) )
  }, my_data = my_data, n = n )
  
  return(one_boot_res)
}

################################################################
#generate data
################################################################
data_generate_model<-function(n, alpha_coef, beta_coef){
  
  v0 <- 1
  alpha_int <- -1; beta_int <- -1
  alpha_cov <- 1; beta_cov <- 1
  tauS <- 1
  
  #1. generate exposure
  S <- rbinom(n, 1, 0.5)
  X <- rbinom(n, 1, 0.5)
  
  #2. generate mediator
  mu_vec <- S * alpha_coef + alpha_int + alpha_cov * X
  mu_logit <- logit_transform_fun(mu_vec)
  M <- rbinom(n, size = 1, prob = mu_logit)
  
  #3. generate outcome
  sigma_sd = 0.5
  mean_y <- M * beta_coef + beta_int + beta_cov * X + tauS * S
  Y <- mean_y + rnorm(n, 0, sigma_sd)
  
  my_data = data.frame(exposure=S, mediator=M, outcome=Y, covariate1 = X, covariate0 = rep(1,n))
  
  return(my_data)
}






