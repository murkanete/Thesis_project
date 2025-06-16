## The following code is adapted from ISTAT implementation (Deliu et al., 2025)

library(fastDummies)
library(MASS)

GMSE = function(sample_data, register_data, covariates, domain, target, C, imp_model) {
  
  X_design = model.matrix(as.formula(paste("~", paste(covariates, collapse = " + "))), 
                     data = register_data)
  
  if (C == 2) {
  
    # estimate the coefficients and p.hat if needed
    #mod_est = glm(formula, data = sample_data, family = "binomial")
    #p.hat = predict(mod_est, newdata = register_data, "response")
    p.hat = register_data[, "pp"]
    
    pi_incl = register_data$pi
    
    # Sigma_matrix
    #Sigma_matrix = lapply(p.hat, function(p) p*(1-p))
    
    # F_term
    F_terms = X_design[,]*p.hat*(1-p.hat)
    
    # list of covariate values to get matrix U
    #G_matrix = list()
    #for(i in 1:N){
    #  G_matrix[[i]] = X_design[i,]
    #}
    
    # A matrix
    A_matrix = t(X_design) %*% (X_design * as.vector(p.hat * (1 - p.hat) * pi_incl))
    invA = ginv(A_matrix)
    
    # evaluate the big sum with Ui (it includes invA, G and Sigma)
    #big_sum = c(0)
    #for(i in 1:N){

      #U_i = invA%*%G_matrix[[i]]
      #sum_i = U_i%*%Sigma_matrix[[i]]%*%t(U_i)*pi_incl[i]
      
      #big_sum = big_sum + sum_i 
      
    #}
    
    big_sum <- (invA %*% t(X_design)) %*% ((X_design %*% invA) * as.vector(p.hat * (1 - p.hat) * pi_incl))
    
    stochastic_noise = p.hat*(1-p.hat)
    
    if(domain == "all"){
      gamma = rep(1, N)
      #theta.hat = sum(p.hat[gamma==1])
      GMSE = gamma%*%F_terms%*%big_sum%*%t(F_terms)%*%gamma
      GMSE_stoc = GMSE + sum(gamma*stochastic_noise)
      colnames(GMSE) = levels(register_data[[target]])[2]
      colnames(GMSE_stoc) = levels(register_data[[target]])[2]
      #CV = sqrt(GMSE)/theta.hat
    } else {
      GMSE = GMSE_stoc = theta.hat = NULL
      for(d in levels(register_data[[domain]])){
        gamma = (register_data[[domain]]==d)*1
        #theta.hat = rbind(theta.hat, sum(p.hat[gamma==1]))
        GMSE_d = c()
        GMSE_d_stoc = c()
        GMSE_d = gamma%*%F_terms%*%big_sum%*%t(F_terms)%*%gamma
        GMSE_d_stoc = GMSE_d + sum(gamma*stochastic_noise)
        GMSE = rbind(GMSE, GMSE_d)
        GMSE_stoc = rbind(GMSE_stoc, GMSE_d_stoc)
      }
      rownames(GMSE) = levels(register_data[[domain]])
      colnames(GMSE) = levels(register_data[[target]])[2]
      
      rownames(GMSE_stoc) = levels(register_data[[domain]])
      colnames(GMSE_stoc) = levels(register_data[[target]])[2]
      
      #CV = sqrt(GMSE)/theta.hat
    }
    
  } else if (C > 2) {
    
    
    # estimate the coefficients and p.hat if needed
    #register_data$y.true_ref <- relevel(register_data$onderwijs, ref = 1)
    #sample_data$y.true_ref = relevel(sample_data$onderwijs, ref = 1)
    #formula = as.formula(paste("y.true_ref ~", paste(covs, collapse = " + ")))
    #mod_est <- multinom(formula, data = sample_data)
    #beta.hat = coef(mod_est)
    #p.hat = predict(mod_est, newdata = register_data, "probs")
    p.hat = as.matrix(register_data[, paste0("p", 1:C)])
    
    beta_hat = coef(imp_model)

    J = length(beta_hat[1,])
    K = length(beta_hat[,1]) + 1
    H = J*K
    
    pi_incl = register_data$pi
    
    # Sigma_matrix
    Sigma_matrix = list()
    for(i in 1:N){
      Sigma_matrix[[i]] = -p.hat[i,]%*%t(p.hat[i,])
      diag(Sigma_matrix[[i]]) = p.hat[i,]*(1-p.hat[i,])
    }
    
    F_terms = list()
    for(k in 1:K){
      F_term_1 = c()
      for(k1 in 1:K){
        F_k = -X_design*p.hat[,k]*p.hat[,k1]
        F_term_1 = cbind(F_term_1, F_k)
      }
      F_term_1[,((k-1)*J+1):(k*J)] = X_design[,]*p.hat[,k]*(1-p.hat[,k])
      F_terms[[k]] = F_term_1
    }
    
    
    
    # now let's compute the G terms of Eq. (16)
    G_matrix = list()
    for(i in 1:N){
      Gi_term = NULL
      for(k in 1:K){
        der.gi_term = rep(0, H) # error adjusted 23/02/2024; previously this line was outside the for loop
        der.gi_term[(J*(k-1)+1):(J*k)] = X_design[i,]
        Gi_term = cbind(Gi_term, der.gi_term)
      }
      G_matrix[[i]] = Gi_term
    }
    
    # A with lambda_i = pi_i
    A_matrix = NULL
    for(k in 1:K){
      A_row_blocks = NULL
      for(j in 1:K){
        A_matrix_blocks = matrix(0, nrow = J, ncol = J)
        for(i in 1:N){
          A_matrix_blocks = -pi_incl[i]*X_design[i,]%*%t(X_design[i,])*Sigma_matrix[[i]][j,k] + A_matrix_blocks
        }
        
        A_row_blocks = rbind(A_row_blocks, A_matrix_blocks)
      }
      
      A_matrix = cbind(A_matrix, A_row_blocks)
    }
    
    #invA = solve(A_matrix)
    invA = ginv(A_matrix)
    
    # evaluate the big sum with Ui (it includes invA, G and Sigma)
    big_sum = c(0)
    for(i in 1:N){
      
      # this is according to Piero version
      U_i = invA%*%G_matrix[[i]]
      sum_i = U_i%*%Sigma_matrix[[i]]%*%t(U_i)*pi_incl[i]
      
      big_sum = big_sum + sum_i 
    }
    
    stochastic_noise = apply(p.hat, 2, function(p.hat) p.hat*(1-p.hat))
    
    if(domain == "all"){
      gamma = rep(1, N)
      #theta.hat = apply(p.hat[gamma==1,], 2, sum)
      GMSE = c()
      GMSE_stoc = c()
      for(k in 1:K){
        GMSE[k] = gamma%*%F_terms[[k]]%*%big_sum%*%t(F_terms[[k]])%*%gamma
        GMSE_stoc[k] = GMSE[k] + sum(gamma*stochastic_noise[, k])
      }
      #CV = sqrt(GMSE)/theta.hat
    } else {
      GMSE = GMSE_stoc = theta.hat = NULL
      for(d in levels(register_data[[domain]])){
        gamma = (register_data[[domain]]==d)*1
        #theta.hat = rbind(theta.hat, apply(p.hat[gamma==1,], 2, sum))
        GMSE_d = c()
        GMSE_stoc_d = c()
        for(k in 1:K){
          GMSE_d[k] = gamma%*%F_terms[[k]]%*%big_sum%*%t(F_terms[[k]])%*%gamma
          GMSE_stoc_d[k] = GMSE_d[k] + sum(gamma*stochastic_noise[, k])
        }
        GMSE = rbind(GMSE, GMSE_d)
        GMSE_stoc = rbind(GMSE_stoc, GMSE_stoc_d)
      }
      rownames(GMSE) = levels(register_data[[domain]])
      colnames(GMSE) = levels(register_data[[target]])
      rownames(GMSE_stoc) = levels(register_data[[domain]])
      colnames(GMSE_stoc) = levels(register_data[[target]])
      #CV = sqrt(GMSE)/theta.hat
    }
  }
  
  return(list(GMSE = GMSE, GMSE_stoc = GMSE_stoc))
  
}
