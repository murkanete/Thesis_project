#install.packages("tidyverse")
library(tidyverse)

PEM_srs = function(sample_data, 
                    register_data,
                    C, 
                    target, 
                    domain,
                    matrix = c("domain-specific", "non-specific")) {
  
  target_imp = paste0(target, "_imp")
  
  if (target_imp %in% colnames(sample_data)) {
    sample_data[[target_imp]] = NULL
  }
  
  sample_data = sample_data %>%
    semi_join(register_data, by = "id")
  
  if (C ==2 ) {
    ids = sample_data$id
    imp_col = register_data[register_data$id %in% ids, c("pp",target_imp, "id")]
    sample_data = merge(sample_data, imp_col, by = "id")
    w = 1/sample_data$pi
    sample_data$w = w
    
    N = nrow(register_data)
    n = nrow(sample_data)
    
    domain_levels = levels(register_data[[domain]])
    target_levels = levels(register_data[[target]])
    
    
    bias = variance = MSE = p11 = p00 = NULL
    
    for (d in domain_levels) {
      
      if (matrix == "domain-specific") {
        sample_data_select = sample_data[sample_data[[domain]]==d, ]
        register_data_select = register_data[register_data[[domain]]==d, ]
        
        class1_samp = sample_data_select[sample_data_select[[target]]==target_levels[2], ]
        class0_samp = sample_data_select[sample_data_select[[target]]==target_levels[1], ]
        
        phat11 = sum(class1_samp$pp)/nrow(class1_samp)
        phat00 = sum(1-class0_samp$pp)/nrow(class0_samp)
        
        Nhat1 = sum(register_data_select[, "pp"])
        Nhat0 = sum(1-register_data_select[, "pp"])
        
        bias_d = Nhat1*(phat11-1) + Nhat0*(1-phat00)
        variance_d = (N/n)*(Nhat1*(phat11*(1-phat11)) +
                              Nhat0*(phat00*(1-phat00)))
        
        MSE_d = bias_d**2 + variance_d
        
        bias = rbind(bias, bias_d)
        variance = rbind(variance, variance_d)
        MSE = rbind(MSE, MSE_d)
        p11 = rbind(phat11, p11)
        p00 = rbind(phat00, p00)
        
        
      } else if (matrix == "non-specific") {
        
        sample_data_select = sample_data[sample_data[[domain]]==d, ]
        register_data_select = register_data[register_data[[domain]]==d, ]
        
        class1_samp = sample_data[sample_data[[target]]==target_levels[2], ]
        class0_samp = sample_data[sample_data[[target]]==target_levels[1], ]
        
        phat11 = sum(class1_samp$pp)/nrow(class1_samp)
        phat00 = sum(1-class0_samp$pp)/nrow(class0_samp)
        
        Nhat1 = sum(register_data_select[, "pp"])
        Nhat0 = sum(1-register_data_select[, "pp"])
        
        bias_d = Nhat1*(phat11-1) + Nhat0*(1-phat00)
        variance_d = (N/n)*(Nhat1*(phat11*(1-phat11)) +
                              Nhat0*(phat00*(1-phat00)))
        
        MSE_d = bias_d**2 + variance_d
        
        bias = rbind(bias, bias_d)
        variance = rbind(variance, variance_d)
        MSE = rbind(MSE, MSE_d)
        p11 = rbind(phat11, p11)
        p00 = rbind(phat00, p00)
      }
      
      
    }
    
    rownames(bias) = domain_levels
    colnames(bias) = target_levels[2]
    rownames(variance) = domain_levels
    colnames(variance) = target_levels[2]
    rownames(MSE) = domain_levels
    colnames(MSE) = target_levels[2]
    rownames(p11) = domain_levels
    colnames(p11) = target_levels[2]
    rownames(p00) = domain_levels
    colnames(p00) = target_levels[1]
    
    return(list(MSE = MSE,
                bias = bias,
                variance = variance,
                p11 = p11,
                p00 = p00))
    
    
  } else if (C > 2) {
    
    pcols = paste0("p", 1:C)
    
    ids = sample_data$id
    imp_col = register_data[register_data$id %in% ids, c(pcols, "id")]
    sample_data = merge(sample_data, imp_col, by = "id")
    w = 1/sample_data$pi
    sample_data$w = w
    
    N = nrow(register_data)
    n = nrow(sample_data)
    
    domain_levels = levels(register_data[[domain]])
    target_levels = levels(register_data[[target]])
    
    
    bias = variance = Nhat = NULL
    P_list = vector("list", length(domain_levels))
    it=0
    
    for (d in domain_levels) {
      
      it=it+1
      sample_data_select = sample_data[sample_data[[domain]]==d, ]
      register_data_select = register_data[register_data[[domain]]==d, ]
      
      P = matrix(data=NA, nrow=C, ncol=C)
      Nhat_d = numeric(C)
      
      #compute P from predicted probabilities
      
      for (g in 1:C) {
        
        classk_samp = sample_data_select[sample_data_select[[target]]==target_levels[g], ]
        Nhat_d[g] = sum(register_data_select[, pcols[g]])
        
        for (k in 1:C) {
          
          P[g,k] = sum(classk_samp[, pcols[k]])/nrow(classk_samp)
          
          
        }
      }
      
      P_list[[it]] = as.matrix(P)
      # compute bias and variance from P
      
      bias_d = variance_d = numeric(C)
      
      for (c in 1:C) {
        probs = P[,c]
        bias_d[c] = (probs[c]-1)*Nhat_d[c] + sum(probs[-c]*Nhat_d[-c])
        variance_d[c] = (N/n)*sum(probs*(1-probs)*Nhat_d)
        
      }
      
      
      
      bias = rbind(bias, bias_d)
      variance = rbind(variance, variance_d)
      Nhat = rbind(Nhat, Nhat_d)
      
    }
    
    rownames(bias) = domain_levels
    colnames(bias) = target_levels
    rownames(variance) = domain_levels
    colnames(variance) = target_levels
    MSE = bias**2 + variance
    
    return(list(bias=bias,
                variance=variance,
                MSE=MSE,
                P=P_list,
                Nhat=Nhat))
    
  }
}


PEM_strat_binom = function(sample_data, 
                          register_data,
                          C, 
                          target, 
                          domain,
                          strat_var) {
  
  target_imp = paste0(target, "_imp")
  
  if (target_imp %in% colnames(sample_data)) {
    sample_data[[target_imp]] = NULL
  }
  
  sample_data = sample_data %>%
    semi_join(register_data, by = "id")
  
  ids = sample_data$id
  imp_col = register_data[register_data$id %in% ids, c("pp",target_imp, "id")]
  sample_data = merge(sample_data, imp_col, by = "id")
  w = 1/sample_data$pi
  sample_data$w = w
  
  
  domain_levels = levels(register_data[[domain]])
  target_levels = levels(register_data[[target]])
  
  
  bias = variance = MSE = p11 = p00 = NULL
  
  for (d in domain_levels) {
    
    sample_data_select = sample_data[sample_data[[domain]]==d, ]
    register_data_select = register_data[register_data[[domain]]==d, ]
      
    class1_samp = sample_data_select[sample_data_select[[target]]==target_levels[2], ]
    class0_samp = sample_data_select[sample_data_select[[target]]==target_levels[1], ]
      
    phat11 = sum(class1_samp$w*class1_samp$pp)/sum(class1_samp$w)
    phat00 = sum(class0_samp$w*(1-class0_samp$pp))/sum(class0_samp$w)
      
    Nhat1 = sum(register_data_select[, "pp"])
    Nhat0 = sum(1-register_data_select[, "pp"])
      
    bias_d = Nhat1*(phat11-1) + Nhat0*(1-phat00)
    
    strat_levels = levels(sample_data[[strat_var]])
    variance_s = numeric(length(strat_levels))
    
    for (stratum in seq_along(strat_levels)) {
      
      
      Nh = summary(register_data[[strat_var]])
      nh = summary(sample_data[[strat_var]])
      
      #class1_samp_strat = class1_samp[class1_samp[[strat_var]] == strat_levels[stratum], ]
      #class0_samp_strat = class0_samp[class0_samp[[strat_var]] == strat_levels[stratum], ]
      
      
      #phat11_strat = sum(class1_samp_strat$pp)/nrow(class1_samp_strat)
      #phat00_strat = sum(1-class0_samp_strat$pp)/nrow(class0_samp_strat)
      
      Nhat1_strat = sum(register_data_select[register_data_select[[strat_var]] == strat_levels[stratum], "pp"])
      Nhat0_strat = sum(1-register_data_select[register_data_select[[strat_var]] == strat_levels[stratum], "pp"])
      
      variance_s[stratum] = (Nh[stratum]/nh[stratum])*(Nhat1_strat*(phat11*(1-phat11)) +
                            Nhat0_strat*(phat00*(1-phat00)))
      
    }
    
    variance_d = sum(variance_s)
      
    MSE_d = bias_d**2 + variance_d
      
    bias = rbind(bias, bias_d)
    variance = rbind(variance, variance_d)
    MSE = rbind(MSE, MSE_d)
    p11 = rbind(phat11, p11)
    p00 = rbind(phat00, p00)
      
  }
  
  rownames(bias) = domain_levels
  colnames(bias) = target_levels[2]
  rownames(variance) = domain_levels
  colnames(variance) = target_levels[2]
  rownames(MSE) = domain_levels
  colnames(MSE) = target_levels[2]
  rownames(p11) = domain_levels
  colnames(p11) = target_levels[2]
  rownames(p00) = domain_levels
  colnames(p00) = target_levels[1]
  
  return(list(MSE = MSE,
              bias = bias,
              variance = variance,
              p11 = p11,
              p00 = p00))
  
  
}

