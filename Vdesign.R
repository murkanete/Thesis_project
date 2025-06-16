#install.packages("data.table")
#install.packages("survey")
#install.packages("svyVGAM")
library(data.table)
library(survey)
library(svyVGAM)

### helper functions

Vmatrix_blocks = function(sandwich) { # function to extract blocks from the sandwich estimator
  
  rows = rownames(sandwich)
  cols = colnames(sandwich)
  
  row_classes = unique(gsub(".*:(\\d+)$", "\\1", rows))
  unique_covariates = unique(sub("(:\\d+)$", "", rows))
  
  out = matrix(data=NA, nrow=length(row_classes)*length(unique_covariates), 
               ncol=length(row_classes)*length(unique_covariates))
  blocks = list()
  
  for (d in 1:length(row_classes)) {
    for (D in 1:length(row_classes)) {
      selected_rows = grep(paste0(":", d, "$"), rows)
      selected_cols = grep(paste0(":", D, "$"), cols)
      
      V_d_D = sandwich[selected_rows, selected_cols, drop = FALSE]
      blocks = append(blocks, list(V_d_D))
      
    }
  }
  return(blocks)
}

delta_function = function(p, c) { # function to compute delta vector
  delta_c = -p[c] * p  
  delta_c[c] = p[c] * (1 - p[c])  
  
  delta_c = delta_c[-1] 
  return(delta_c)
}

### Vdesign

Vdesign = function(sample_data, register_data, covariates, domain, target, C, design = c("SRS", "Stratified")) {
  
  register_data = as.data.table(register_data)
  
  if (C == 2) {
    # first variance term
    register_data[, v1 := pp * (1 - pp)]
    
    # second variance term
    
    ## sandwich
    n = nrow(register_data)
    N = nrow(sample_data)
    
    if (design == "SRS") {
      rownames(sample_data) = NULL
      sample_data$id = 1:nrow(sample_data)
      sample_data$w = 1/sample_data$pi
      survey = svydesign(id ~ 1, weights = ~w, data = sample_data)
      model = svyglm(as.formula(paste(target, "~", paste(covariates, collapse = " + "))),
                     family = 'binomial', design = survey)
      sandwich = vcov(model)
    } else if(design == "Stratified") {
      rownames(sample_data) = NULL
      sample_data$id = 1:nrow(sample_data)
      sample_data$w = 1/sample_data$pi
      survey = svydesign(id ~ 1, strata = ~strat_var, weights = ~w, data = sample_data)
      model = svyglm(as.formula(paste(target, "~", paste(covariates, collapse = " + "))),
                     family = 'binomial', design = survey)
      sandwich = vcov(model)
    }
    
    ## create dummies
    register_select = model.matrix(as.formula(paste("~", paste(covariates, collapse = " + "))), data = register_data)
    dummy_covs = colnames(register_select)
    register_select = as.data.table(register_select)
    
    register_select[, (domain) := register_data[[domain]]]
    register_select[, v1 := register_data$v1]
    
    if (domain == "all") {
      
      i_cov_values = as.matrix(register_select[, dummy_covs, with = FALSE])
      
      u_tot = colSums(as.numeric(register_select[, mget(v1)]) * i_cov_values[1:N_domain,])
      
      covpcpc = t(as.numeric(u_tot)) %*% sandwich %*% as.numeric(u_tot)
      v1_comp = sum(register_select[, mget(v1)])
      
      Vdesign_out = as.numeric(v1_comp) + covpcpc
      
      
    } else {
      domain_levels = levels(register_data[[domain]])
      
      subsets = vector("list", length(domain_levels))
      for (i in 1:length(domain_levels)) {
        subsets[[i]] = register_select[get(domain) == domain_levels[i], ]
      }

      Vdesign_out = NULL
      
      for (d in 1:length(subsets)) {
        N_domain = nrow(subsets[[d]])
        register_domain = subsets[[d]] 
        
        # covariate matrix
        i_cov_values = as.matrix(register_domain[, dummy_covs, with = FALSE])
        
        # u-tot vector
        u_tot = colSums(as.numeric(register_domain$v1) * i_cov_values[1:N_domain,])
        
        # second and first variance term
        Vdesign_out_d = c()
        covpcpc_d = t(as.numeric(u_tot)) %*% sandwich %*% as.numeric(u_tot)
        v1_comp = register_data[get(domain) == domain_levels[d], .(v1 = sum(v1))]
        Vdesign_out_d = as.numeric(v1_comp) + covpcpc_d
        Vdesign_out = rbind(Vdesign_out, Vdesign_out_d)
        
      }
      rownames(Vdesign_out) = domain_levels
      colnames(Vdesign_out) = levels(register_data[[target]])[2]
    }
    
  } else if (C > 2) {
    
    # first variance term
    for (c in 1:C) eval(parse(text = sprintf('register_data[ , pp%d := p%d * (1 - p%d)]', c, c, c)))
    
    # second variance term
    
    ## sandwich
    rownames(sample_data) = NULL
    sample_data$id = 1:nrow(sample_data)
    survey = svydesign(id ~ 1, fpc = rep(N, n), data = sample_data)
    model = svy_vglm(as.formula(paste(target, "~", paste(covariates, collapse = " + "))),
                     family = multinomial(refLevel = 1), design = survey)
    sandwich = vcov(model)
    vblocks = Vmatrix_blocks(sandwich)
    
    register_select = model.matrix(as.formula(paste("~", paste(covariates, collapse = " + "))), data = register_data)
    dummy_covs = colnames(register_select)
    register_select = cbind(register_select, register_data[, mget(paste0("p", 1:C))])
    register_select = cbind(register_select, register_data[, mget(paste0("pp", 1:C))])
    register_select = as.data.table(register_select)
    
    if (domain == "all") {
      
      # covariate matrix
      i_cov_values = as.matrix(register_select[, dummy_covs, with = FALSE])
      # predicted probabilities
      i_p_values = as.matrix(register_select[, paste0("p", 1:C), with = FALSE])
      
      # delta values
      delta_ip = sapply(1:C, function(c) {
        apply(i_p_values, 1, function(p) delta_function(p = as.numeric(p), c = c))
      }, simplify = "array")
      
      # u-tot vector
      u_tot = vector("list", C)
      for (c in 1:C) {
        d_list = vector("list", C-1)
        for (d in 1:(C-1)) {
          d_list[[d]] = colSums(delta_ip[d, 1:nrow(register_select), c] * i_cov_values[1:nrow(register_select), ])
        }
        u_tot[[c]] = d_list
      }
      
      # second and first variance term
      Vdesign_out = numeric(C)
      for (c in 1:C) {
        covpcpc = numeric((C-1)*(C-1))
        for (b in 1:(C-1)) {
          for (bhat in 1:(C-1)) {
            iteration = (b - 1) * (C-1) + bhat
            sum_bbhat = t(as.numeric(u_tot[[c]][[b]])) %*% vblocks[[iteration]] %*% as.numeric(u_tot[[c]][[bhat]])
            covpcpc[iteration] = sum_bbhat
          }
        }
        Vdesign_out[c] = sum(covpcpc) + sum(register_select[,mget(paste0("pp", c))])
      }
      names(Vdesign_out) = levels(register_data[[target]])
      
    } else {
      register_select[, (domain) := register_data[[domain]]]
      domain_levels = levels(register_data[[domain]])
      
      subsets = vector("list", length(domain_levels))
      for (i in 1:length(domain_levels)) {
        subsets[[i]] = register_select[get(domain) == domain_levels[i], ]
      }
      
      Vdesign_out = NULL
      
      for (d in 1:length(subsets)) {

        N_domain = nrow(subsets[[d]])
        register_domain = subsets[[d]] 
        
        # covariate matrix
        i_cov_values = as.matrix(register_domain[, dummy_covs, with = FALSE])
        # predicted probabilities
        i_p_values = as.matrix(register_domain[, paste0("p", 1:C), with = FALSE])
        
        # delta values
        delta_ip = sapply(1:C, function(c) {
          apply(i_p_values, 1, function(p) delta_function(p = as.numeric(p), c = c))
        }, simplify = "array")
        
        # u-tot vector
        u_tot = vector("list", C)
        for (c in 1:C) {
          d_list = vector("list", C-1)
          for (delta in 1:(C-1)) {
            d_list[[delta]] = colSums(delta_ip[delta, 1:N_domain, c] * i_cov_values[1:N_domain, ])
          }
          u_tot[[c]] = d_list
        }
        
        # second and first variance term
        Vdesign_out_d = numeric(C)
        for (c in 1:C) {
          covpcpc = numeric((C-1)*(C-1))
          for (b in 1:(C-1)) {
            for (bhat in 1:(C-1)) {
              iteration = (b-1)*(C-1)+bhat
              sum_bbhat = t(as.numeric(u_tot[[c]][[b]])) %*% vblocks[[iteration]] %*% as.numeric(u_tot[[c]][[bhat]])
              covpcpc[iteration] = sum_bbhat
            }
          }
          Vdesign_out_d[c] = sum(covpcpc) + sum(register_domain[,mget(paste0("pp", c))])
        }
        Vdesign_out = rbind(Vdesign_out, Vdesign_out_d)
      }
      rownames(Vdesign_out) = domain_levels
      colnames(Vdesign_out) = levels(register_data[[target]])
    }
  }
  return(Vdesign_out)
  }
