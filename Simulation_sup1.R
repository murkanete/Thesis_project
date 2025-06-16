##########################################################################################

# The following script performs the simulation study at each experimental condition
# in superpopulation 1. The results from each condition are saved to disc and combined
# at the end of the script.

# The Samplonia datafile and function scripts should be in the same working directory.

##########################################################################################


#libraries and helper functions

#install.packages("nnet")
library(nnet)
source("PEM_functions.R") 
source("GMSE.R")
source("Vdesign.R")

###########  Define the superpopulation at baseline settings   ##########


##### Transform the original Samplonia file

pop <- read.csv2(file = 'samplonie.csv', header = TRUE,
                 stringsAsFactors = TRUE)

pop <- pop[pop$leeftijd >= 15, ] #select people over 15


#turn income into categorical variable
pop$inkomen3 <- cut(pop$inkomen, 
                    breaks = c(0, 350, 1400, 5000), 
                    include.lowest = TRUE, 
                    labels = c("Low", "Middle", "High"))


# large N 
pop <- pop[rep(1:nrow(pop), times = 125), ]

# binomial target variable with balanced classes
c = 2

pop$onderwijs <- ifelse(pop$onderwijs %in% c('Geen','Basisonderwijs', 'Vmbo','Havo/vwo','Mbo'), 1, 2)
pop$onderwijs <- factor(pop$onderwijs, labels = c('Geen wo','Wo'))


# oversample Wo class to create balanced classes 
set.seed(4789)
n_to_sample = sum(pop$onderwijs == "Geen wo") - sum(pop$onderwijs == "Wo") # nr extra Wo samples needed
wo_subset = pop[pop$onderwijs == "Wo", ]
extra_samples = wo_subset[sample(nrow(wo_subset), n_to_sample, replace = TRUE), ] # sample with replacement
pop_balanced = rbind(pop, extra_samples)
rm(extra_samples, wo_subset)

N_large = nrow(pop_balanced) # N = 122 250
n_large = round(0.05*N_large)  # n = 6112


# define inclusion probabilities
pop_balanced$pi = n_large/N_large
pop_balanced$id = 1:nrow(pop_balanced)



######## Define the superpopulation 


# log-scale the numerical variable
pop_balanced$leeftijd = scale(log(pop_balanced$leeftijd+0.0001))

covs_1 = c("geslacht", "leeftijd", "inkomen3") # geslacht will be the domain variable
formula_1 = as.formula(paste("onderwijs ~", paste(covs_1, collapse = " + "))) 
superpopulation_1 = glm(formula_1, data = pop_balanced, family='binomial')



##########################################################################################

                           ######## Simulation  #########


##########################################################################################

                               # Baseline condition #


##########################################################################################

set.seed(4789)
b_samp = 100
B_pop = 100

# define the items to store for each finite population; so storing the model distribution
bias_model = variance_model = theta_model = theta_hat_model = MSE_model = CSRM_MSE_model = 
  vdesign_model = CSRM_vdesign_model = gmse_model = CSRM_gmse_model =
  pem_mse_model = CSRM_pem_mse_model = pem_bias_model = pem_variance_model = CSRM_pem_bias_model = CSRM_pem_variance_model = 
  P_model = CSRM_bias_model = CSRM_variance_model =
  vector("list", B_pop)

for (B in 1:B_pop) {
  
  print(paste(B, "Baseline"))
  
  # generate y in register; model variability
  sim_register = pop_balanced
  p.true = predict(superpopulation_1, newdata = sim_register, type="response")
  sim_register$onderwijs = factor(sapply(p.true, function(p) rbinom(1, 1, p)),
                                  labels = c("Geen wo", "Wo"))
  
  theta_model[[B]] = unclass(as.matrix(
    table(sim_register$geslacht, sim_register$onderwijs)[, 2, drop=F])) # matrix object for later data manipulation
  
  # store results from each sample drawn from the finite population; sampling distribution
  difference_sampling = theta_hat_sampling = vdesign_sampling = CSRM_vdesign_sampling =
    gmse_sampling = CSRM_gmse_sampling =
    pem_mse_sampling = pem_bias_sampling = pem_variance_sampling = CSRM_pem_mse_sampling = 
    CSRM_pem_bias_sampling = 
    CSRM_pem_variance_sampling = P_sampling =
    vector("list", b_samp)
  
  
  for (b in 1:b_samp) { # sampling variability
    
    print(b)
    
    # draw an SRS without replacement
    S = sample(N_large, n_large)
    samp = sim_register[S,]
    
    # train an imputation model
    imp_model = glm(formula_1, data = samp, family='binomial') # parameter variability
    
    # imputation
    predicted_probs = predict(imp_model, newdata = sim_register, type = "response")
    sim_register$onderwijs_imp = factor(sapply(predicted_probs, function(p) rbinom(1, 1, p)),
                                        labels = c("Geen wo", "Wo")) # imputation variability
    sim_register_b = sim_register
    sim_register_b$pp = predicted_probs
    
    
    # for benchmarks
    theta_hat_sampling[[b]] = unclass(as.matrix(
      table(sim_register$geslacht, sim_register$onderwijs_imp)[, 2, drop=F]))
    difference_sampling[[b]] = theta_hat_sampling[[b]] - theta_model[[B]]
    
    
    # accuracy estimators
    vdesign_sampling[[b]] = Vdesign(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_1,
      domain = "geslacht",
      target = "onderwijs",
      C = c,
      design="SRS"
    )
    
    CSRM_vdesign_sampling[[b]] = sqrt(vdesign_sampling[[b]])/theta_hat_sampling[[b]]
    
    gmse_results = GMSE(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_1,
      domain = "geslacht",
      target = "onderwijs",
      C = c,
      imp_model = imp_model
    )
    
    gmse_sampling[[b]] = gmse_results$GMSE_stoc
    CSRM_gmse_sampling[[b]] = sqrt(gmse_sampling[[b]])/theta_hat_sampling[[b]]
    
    PEM_results = PEM_srs(
      sample_data = samp,
      register_data = sim_register_b,
      C = c,
      target = "onderwijs",
      domain = "geslacht",
      matrix = "domain-specific"
    )
    
    
    pem_mse_sampling[[b]] = PEM_results$MSE
    pem_bias_sampling[[b]] = PEM_results$bias
    pem_variance_sampling[[b]] = PEM_results$variance
    P_sampling[[b]] = cbind(PEM_results$p00, PEM_results$p11)
    
    CSRM_pem_mse_sampling[[b]] = sqrt(pem_mse_sampling[[b]])/theta_hat_sampling[[b]]
    CSRM_pem_bias_sampling[[b]] = pem_bias_sampling[[b]]/theta_hat_sampling[[b]]
    CSRM_pem_variance_sampling[[b]] = sqrt(pem_variance_sampling[[b]])/theta_hat_sampling[[b]]
    
  }
  
  # benchmark across b samples
  theta_hat_model[[B]] = apply(simplify2array(theta_hat_sampling), c(1,2), mean)
  difference_sampling = simplify2array(difference_sampling)
  bias_model[[B]] = apply(difference_sampling, c(1, 2), mean)
  CSRM_bias_model[[B]] = bias_model[[B]]/theta_model[[B]]
  variance_model[[B]] = apply(difference_sampling, c(1,2), var)
  CSRM_variance_model[[B]] = sqrt(variance_model[[B]])/theta_model[[B]]
  MSE_model[[B]] = bias_model[[B]]**2 + variance_model[[B]]
  CSRM_MSE_model[[B]] = sqrt(MSE_model[[B]])/theta_model[[B]]
  
  # accuracy estimators across b samples
  vdesign_model[[B]] = apply(simplify2array(vdesign_sampling), c(1, 2), mean)
  CSRM_vdesign_model[[B]] = apply(simplify2array(CSRM_vdesign_sampling), c(1, 2), mean)
  
  gmse_model[[B]] = apply(simplify2array(gmse_sampling), c(1,2), mean)
  CSRM_gmse_model[[B]] = apply(simplify2array(CSRM_gmse_sampling), c(1,2), mean)
  
  
  pem_mse_model[[B]] = apply(simplify2array(pem_mse_sampling), c(1,2), mean)
  CSRM_pem_mse_model[[B]] = apply(simplify2array(CSRM_pem_mse_sampling), c(1,2), mean)
  pem_bias_model[[B]] = apply(simplify2array(pem_bias_sampling), c(1,2), mean)
  pem_variance_model[[B]] = apply(simplify2array(pem_variance_sampling), c(1,2), mean)
  CSRM_pem_bias_model[[B]] = apply(simplify2array(CSRM_pem_bias_sampling), c(1,2), mean)
  CSRM_pem_variance_model[[B]] = apply(simplify2array(CSRM_pem_variance_sampling), c(1,2), mean)
  
  P_model[[B]] = apply(simplify2array(P_sampling), c(1,2), mean)
  
}

results_model_names = c(
  "vdesign_model", "CSRM_vdesign_model", "gmse_model", "CSRM_gmse_model",
  "pem_mse_model", "CSRM_pem_mse_model", "pem_bias_model", "pem_variance_model", "CSRM_pem_bias_model",
  "CSRM_pem_variance_model", "bias_model", "variance_model", "theta_model", 
  "theta_hat_model", "MSE_model", "CSRM_MSE_model",
  "CSRM_bias_model", "CSRM_variance_model"
)

results_sampling_names = c(
  "theta_hat_sampling", "vdesign_sampling", "CSRM_vdesign_sampling",
  "gmse_sampling", "CSRM_gmse_sampling",
  "pem_mse_sampling", "pem_bias_sampling", "pem_variance_sampling", 
  "CSRM_pem_mse_sampling",
  "CSRM_pem_bias_sampling", 
  "CSRM_pem_variance_sampling"
)

results_dfs = function(results_model_names, results_sampling_names) {
  
  list_of_results_model = mget(results_model_names, envir = .GlobalEnv)
  
  results_model = list()
  
  for (i in seq_along(list_of_results_model)) {
    result = list_of_results_model[[i]]
    result_name = results_model_names[i]
    
    result_df = bind_rows(lapply(seq_along(result), function(j) {
      df = as.data.frame(result[[j]])
      df$Iteration = j
      df$Estimate = result_name
      df$Domain = rownames(df)
      rownames(df) = NULL
      return(df)
    }))
    
    results_model[[i]] = result_df
  }
  
  results_model_df = bind_rows(results_model)  # model distribution
  results_joint_df = results_model_df %>% group_by(Estimate, Domain) %>%  # joint distribution
    summarise(joint = mean(Wo))
  
  
  
  list_results_sampling = mget(results_sampling_names, envir = .GlobalEnv)
  results_sampling = list()
  
  for (i in seq_along(list_results_sampling)) {
    result = list_results_sampling[[i]]
    result_name = results_sampling_names[i]
    
    result_df = bind_rows(lapply(seq_along(result), function(j) {
      df = as.data.frame(result[[j]])
      df$Iteration = j
      df$Estimate = result_name
      df$Domain = rownames(df)
      rownames(df) = NULL
      return(df)
    }))
    
    results_sampling[[i]] = result_df
  }
  
  results_sampling_df = bind_rows(results_sampling)  # sampling distribution; last iteration of B
  
  return(list(
    results_sampling_df = results_sampling_df,
    results_model_df = results_model_df,
    results_joint_df = results_joint_df
  ))
  
}

results = results_dfs(results_model_names, results_sampling_names)


# write results to disk
save(P_model, P_sampling, file = "sup1_baseline_P.RData")
write.csv(results$results_sampling_df, "sup1_baseline_sampling_distribution.csv")
write.csv(results$results_model_df, "sup1_baseline_model_distribution.csv")
write.csv(results$results_joint_df, "sup1_baseline_joint_distribution.csv")


##########################################################################################

                             # Small sample condition #


##########################################################################################


n_small = round(0.01*N_large)
pop_balanced$pi = n_small/N_large

set.seed(4789)

# define the items to store for each finite population; so storing the model distribution
bias_model = variance_model = theta_model = theta_hat_model = MSE_model = CSRM_MSE_model = 
  vdesign_model = CSRM_vdesign_model = gmse_model = CSRM_gmse_model =
  pem_mse_model = CSRM_pem_mse_model = pem_bias_model = pem_variance_model = CSRM_pem_bias_model = CSRM_pem_variance_model = 
  P_model = CSRM_bias_model = CSRM_variance_model =
  vector("list", B_pop)

for (B in 1:B_pop) {
  
  print(paste(B, "Small sample"))
  
  # generate y in register; model variability
  sim_register = pop_balanced
  p.true = predict(superpopulation_1, newdata = sim_register, type="response")
  sim_register$onderwijs = factor(sapply(p.true, function(p) rbinom(1, 1, p)),
                                  labels = c("Geen wo", "Wo"))
  
  theta_model[[B]] = unclass(as.matrix(
    table(sim_register$geslacht, sim_register$onderwijs)[, 2, drop=F])) # matrix object for later data manipulation
  
  # store results from each sample drawn from the finite population; sampling distribution
  difference_sampling = theta_hat_sampling = vdesign_sampling = CSRM_vdesign_sampling =
    gmse_sampling = CSRM_gmse_sampling =
    pem_mse_sampling = pem_bias_sampling = pem_variance_sampling = CSRM_pem_mse_sampling = 
    CSRM_pem_bias_sampling = 
    CSRM_pem_variance_sampling = P_sampling =
    vector("list", b_samp)
  
  
  for (b in 1:b_samp) { # sampling variability
    
    print(b)
    
    # draw an SRS without replacement
    S = sample(N_large, n_small)
    samp = sim_register[S,]
    
    # train an imputation model
    imp_model = glm(formula_1, data = samp, family='binomial') # parameter variability
    
    # imputation
    predicted_probs = predict(imp_model, newdata = sim_register, type = "response")
    sim_register$onderwijs_imp = factor(sapply(predicted_probs, function(p) rbinom(1, 1, p)),
                                        labels = c("Geen wo", "Wo")) # imputation variability
    sim_register_b = sim_register
    sim_register_b$pp = predicted_probs
    
    
    # for benchmarks
    theta_hat_sampling[[b]] = unclass(as.matrix(
      table(sim_register$geslacht, sim_register$onderwijs_imp)[, 2, drop=F]))
    difference_sampling[[b]] = theta_hat_sampling[[b]] - theta_model[[B]]
    
    
    # accuracy estimators
    vdesign_sampling[[b]] = Vdesign(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_1,
      domain = "geslacht",
      target = "onderwijs",
      C = c,
      design="SRS"
    )
    
    CSRM_vdesign_sampling[[b]] = sqrt(vdesign_sampling[[b]])/theta_hat_sampling[[b]]
    
    gmse_results = GMSE(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_1,
      domain = "geslacht",
      target = "onderwijs",
      C = c,
      imp_model = imp_model
    )
    
    gmse_sampling[[b]] = gmse_results$GMSE_stoc
    CSRM_gmse_sampling[[b]] = sqrt(gmse_sampling[[b]])/theta_hat_sampling[[b]]
    
    PEM_results = PEM_srs(
      sample_data = samp,
      register_data = sim_register_b,
      C = c,
      target = "onderwijs",
      domain = "geslacht",
      matrix = "domain-specific"
    )
    
    
    pem_mse_sampling[[b]] = PEM_results$MSE
    pem_bias_sampling[[b]] = PEM_results$bias
    pem_variance_sampling[[b]] = PEM_results$variance
    P_sampling[[b]] = cbind(PEM_results$p00, PEM_results$p11)
    
    CSRM_pem_mse_sampling[[b]] = sqrt(pem_mse_sampling[[b]])/theta_hat_sampling[[b]]
    CSRM_pem_bias_sampling[[b]] = pem_bias_sampling[[b]]/theta_hat_sampling[[b]]
    CSRM_pem_variance_sampling[[b]] = sqrt(pem_variance_sampling[[b]])/theta_hat_sampling[[b]]
    
  }
  
  # benchmark across b samples
  theta_hat_model[[B]] = apply(simplify2array(theta_hat_sampling), c(1,2), mean)
  difference_sampling = simplify2array(difference_sampling)
  bias_model[[B]] = apply(difference_sampling, c(1, 2), mean)
  CSRM_bias_model[[B]] = bias_model[[B]]/theta_model[[B]]
  variance_model[[B]] = apply(difference_sampling, c(1,2), var)
  CSRM_variance_model[[B]] = sqrt(variance_model[[B]])/theta_model[[B]]
  MSE_model[[B]] = bias_model[[B]]**2 + variance_model[[B]]
  CSRM_MSE_model[[B]] = sqrt(MSE_model[[B]])/theta_model[[B]]
  
  # accuracy estimators across b samples
  vdesign_model[[B]] = apply(simplify2array(vdesign_sampling), c(1, 2), mean)
  CSRM_vdesign_model[[B]] = apply(simplify2array(CSRM_vdesign_sampling), c(1, 2), mean)
  
  gmse_model[[B]] = apply(simplify2array(gmse_sampling), c(1,2), mean)
  CSRM_gmse_model[[B]] = apply(simplify2array(CSRM_gmse_sampling), c(1,2), mean)
  
  
  pem_mse_model[[B]] = apply(simplify2array(pem_mse_sampling), c(1,2), mean)
  CSRM_pem_mse_model[[B]] = apply(simplify2array(CSRM_pem_mse_sampling), c(1,2), mean)
  pem_bias_model[[B]] = apply(simplify2array(pem_bias_sampling), c(1,2), mean)
  pem_variance_model[[B]] = apply(simplify2array(pem_variance_sampling), c(1,2), mean)
  CSRM_pem_bias_model[[B]] = apply(simplify2array(CSRM_pem_bias_sampling), c(1,2), mean)
  CSRM_pem_variance_model[[B]] = apply(simplify2array(CSRM_pem_variance_sampling), c(1,2), mean)
  
  P_model[[B]] = apply(simplify2array(P_sampling), c(1,2), mean)
  
}

results = results_dfs(results_model_names, results_sampling_names)


# write results to disk
save(P_model, P_sampling, file = "sup1_small_sample_P.RData")
write.csv(results$results_sampling_df, "sup1_small_sample_sampling_distribution.csv")
write.csv(results$results_model_df, "sup1_small_sample_model_distribution.csv")
write.csv(results$results_joint_df, "sup1_small_sample_joint_distribution.csv")


##########################################################################################

                              # Small population condition #


##########################################################################################


set.seed(4789)

pop <- read.csv2(file = 'samplonie.csv', header = TRUE,
                 stringsAsFactors = TRUE)

pop <- pop[pop$leeftijd >= 15, ] #select people over 15


#turn income into categorical variable
pop$inkomen3 <- cut(pop$inkomen, 
                    breaks = c(0, 350, 1400, 5000), 
                    include.lowest = TRUE, 
                    labels = c("Low", "Middle", "High"))


# large N 
pop <- pop[rep(1:nrow(pop), times = 50), ]

# binomial target variable with balanced classes
c = 2

pop$onderwijs <- ifelse(pop$onderwijs %in% c('Geen','Basisonderwijs', 'Vmbo','Havo/vwo','Mbo'), 1, 2)
pop$onderwijs <- factor(pop$onderwijs, labels = c('Geen wo','Wo'))

#prop.table(table(pop$onderwijs)) # the classes are moderately imbalanced in the original data; Wo will be oversampled to result in 50/50
#prop.table(table(pop$geslacht, pop$onderwijs), margin = 1) 

# oversample Wo class to create balanced classes 
n_to_sample = sum(pop$onderwijs == "Geen wo") - sum(pop$onderwijs == "Wo") # nr extra Wo samples needed
wo_subset = pop[pop$onderwijs == "Wo", ]
extra_samples = wo_subset[sample(nrow(wo_subset), n_to_sample, replace = TRUE), ] # sample with replacement
pop_balanced_small = rbind(pop, extra_samples)
rm(extra_samples, wo_subset)

N_small = nrow(pop_balanced_small) 
n_large_2 = round(0.05*N_small) 


pop_balanced_small$pi = n_large_2/N_small
pop_balanced_small$id = 1:nrow(pop_balanced_small)

# log-scale the numerical variable
pop_balanced_small$leeftijd = scale(log(pop_balanced_small$leeftijd+0.0001))


superpopulation_smallpop = glm(formula_1, data = pop_balanced_small, family='binomial')


# define the items to store for each finite population; so storing the model distribution
bias_model = variance_model = theta_model = theta_hat_model = MSE_model = CSRM_MSE_model = 
  vdesign_model = CSRM_vdesign_model = gmse_model = CSRM_gmse_model =
  pem_mse_model = CSRM_pem_mse_model = pem_bias_model = pem_variance_model = CSRM_pem_bias_model = CSRM_pem_variance_model = 
  P_model = CSRM_bias_model = CSRM_variance_model =
  vector("list", B_pop)

for (B in 1:B_pop) {
  
  print(paste(B, "Small population"))
  
  # generate y in register; model variability
  sim_register = pop_balanced_small
  p.true = predict(superpopulation_smallpop, newdata = sim_register, type="response")
  sim_register$onderwijs = factor(sapply(p.true, function(p) rbinom(1, 1, p)),
                                  labels = c("Geen wo", "Wo"))
  
  theta_model[[B]] = unclass(as.matrix(
    table(sim_register$geslacht, sim_register$onderwijs)[, 2, drop=F])) # matrix object for later data manipulation
  
  # store results from each sample drawn from the finite population; sampling distribution
  difference_sampling = theta_hat_sampling = vdesign_sampling = CSRM_vdesign_sampling =
    gmse_sampling = CSRM_gmse_sampling =
    pem_mse_sampling = pem_bias_sampling = pem_variance_sampling = CSRM_pem_mse_sampling = 
    CSRM_pem_bias_sampling = 
    CSRM_pem_variance_sampling = P_sampling =
    vector("list", b_samp)
  
  
  for (b in 1:b_samp) { # sampling variability
    
    print(b)
    
    # draw an SRS without replacement
    S = sample(N_small, n_large_2)
    samp = sim_register[S,]
    
    # train an imputation model
    imp_model = glm(formula_1, data = samp, family='binomial') # parameter variability
    
    # imputation
    predicted_probs = predict(imp_model, newdata = sim_register, type = "response")
    sim_register$onderwijs_imp = factor(sapply(predicted_probs, function(p) rbinom(1, 1, p)),
                                        labels = c("Geen wo", "Wo")) # imputation variability
    sim_register_b = sim_register
    sim_register_b$pp = predicted_probs
    
    
    # for benchmarks
    theta_hat_sampling[[b]] = unclass(as.matrix(
      table(sim_register$geslacht, sim_register$onderwijs_imp)[, 2, drop=F]))
    difference_sampling[[b]] = theta_hat_sampling[[b]] - theta_model[[B]]
    
    
    # accuracy estimators
    vdesign_sampling[[b]] = Vdesign(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_1,
      domain = "geslacht",
      target = "onderwijs",
      C = c,
      design="SRS"
    )
    
    CSRM_vdesign_sampling[[b]] = sqrt(vdesign_sampling[[b]])/theta_hat_sampling[[b]]
    
    gmse_results = GMSE(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_1,
      domain = "geslacht",
      target = "onderwijs",
      C = c,
      imp_model = imp_model
    )
    
    gmse_sampling[[b]] = gmse_results$GMSE_stoc
    CSRM_gmse_sampling[[b]] = sqrt(gmse_sampling[[b]])/theta_hat_sampling[[b]]
    
    PEM_results = PEM_srs(
      sample_data = samp,
      register_data = sim_register_b,
      C = c,
      target = "onderwijs",
      domain = "geslacht",
      matrix = "domain-specific"
    )
    
    
    pem_mse_sampling[[b]] = PEM_results$MSE
    pem_bias_sampling[[b]] = PEM_results$bias
    pem_variance_sampling[[b]] = PEM_results$variance
    P_sampling[[b]] = cbind(PEM_results$p00, PEM_results$p11)
    
    CSRM_pem_mse_sampling[[b]] = sqrt(pem_mse_sampling[[b]])/theta_hat_sampling[[b]]
    CSRM_pem_bias_sampling[[b]] = pem_bias_sampling[[b]]/theta_hat_sampling[[b]]
    CSRM_pem_variance_sampling[[b]] = sqrt(pem_variance_sampling[[b]])/theta_hat_sampling[[b]]
    
  }
  
  # benchmark across b samples
  theta_hat_model[[B]] = apply(simplify2array(theta_hat_sampling), c(1,2), mean)
  difference_sampling = simplify2array(difference_sampling)
  bias_model[[B]] = apply(difference_sampling, c(1, 2), mean)
  CSRM_bias_model[[B]] = bias_model[[B]]/theta_model[[B]]
  variance_model[[B]] = apply(difference_sampling, c(1,2), var)
  CSRM_variance_model[[B]] = sqrt(variance_model[[B]])/theta_model[[B]]
  MSE_model[[B]] = bias_model[[B]]**2 + variance_model[[B]]
  CSRM_MSE_model[[B]] = sqrt(MSE_model[[B]])/theta_model[[B]]
  
  # accuracy estimators across b samples
  vdesign_model[[B]] = apply(simplify2array(vdesign_sampling), c(1, 2), mean)
  CSRM_vdesign_model[[B]] = apply(simplify2array(CSRM_vdesign_sampling), c(1, 2), mean)
  
  gmse_model[[B]] = apply(simplify2array(gmse_sampling), c(1,2), mean)
  CSRM_gmse_model[[B]] = apply(simplify2array(CSRM_gmse_sampling), c(1,2), mean)
  
  
  pem_mse_model[[B]] = apply(simplify2array(pem_mse_sampling), c(1,2), mean)
  CSRM_pem_mse_model[[B]] = apply(simplify2array(CSRM_pem_mse_sampling), c(1,2), mean)
  pem_bias_model[[B]] = apply(simplify2array(pem_bias_sampling), c(1,2), mean)
  pem_variance_model[[B]] = apply(simplify2array(pem_variance_sampling), c(1,2), mean)
  CSRM_pem_bias_model[[B]] = apply(simplify2array(CSRM_pem_bias_sampling), c(1,2), mean)
  CSRM_pem_variance_model[[B]] = apply(simplify2array(CSRM_pem_variance_sampling), c(1,2), mean)
  
  P_model[[B]] = apply(simplify2array(P_sampling), c(1,2), mean)
  
}


results = results_dfs(results_model_names, results_sampling_names)

# write results to disk
save(P_model, P_sampling, file = "sup1_small_population_P.RData")
write.csv(results$results_sampling_df, "sup1_small_population_sampling_distribution.csv")
write.csv(results$results_model_df, "sup1_small_population_model_distribution.csv")
write.csv(results$results_joint_df, "sup1_small_population_joint_distribution.csv")


##########################################################################################

                               # Overparameterisation condition #


##########################################################################################

set.seed(4789)

pop_balanced$pi = n_large/N_large

covs_imp = c("geslacht", "leeftijd", "inkomen3", "gemeente", "random")
formula_imp = as.formula(paste("onderwijs ~", paste(covs_imp, collapse = " + "))) 

assign_random <- function(edu) {
  if (edu == "Wo") {
    sample(c("A", "B", "C", "D"), size = 1, prob = c(0.25, 0.25, 0.25, 0.25))  
  } else {
    sample(c("A", "B", "C", "D"), size = 1, prob = c(0.25, 0.25, 0.25, 0.25))  
  }
}

pop_balanced$random <- sapply(pop_balanced$onderwijs, assign_random)
pop_balanced$random = as.factor(pop_balanced$random)

# define the items to store for each finite population; so storing the model distribution
bias_model = variance_model = theta_model = theta_hat_model = MSE_model = CSRM_MSE_model = 
  vdesign_model = CSRM_vdesign_model = gmse_model = CSRM_gmse_model =
  pem_mse_model = CSRM_pem_mse_model = pem_bias_model = pem_variance_model = CSRM_pem_bias_model = CSRM_pem_variance_model = 
  P_model = CSRM_bias_model = CSRM_variance_model =
  vector("list", B_pop)

for (B in 1:B_pop) {
  
  print(paste(B, "Overparameterisation"))
  
  # generate y in register; model variability
  sim_register = pop_balanced
  p.true = predict(superpopulation_1, newdata = sim_register, type="response")
  sim_register$onderwijs = factor(sapply(p.true, function(p) rbinom(1, 1, p)),
                                  labels = c("Geen wo", "Wo"))
  
  theta_model[[B]] = unclass(as.matrix(
    table(sim_register$geslacht, sim_register$onderwijs)[, 2, drop=F])) # matrix object for later data manipulation
  
  # store results from each sample drawn from the finite population; sampling distribution
  difference_sampling = theta_hat_sampling = vdesign_sampling = CSRM_vdesign_sampling =
    gmse_sampling = CSRM_gmse_sampling =
    pem_mse_sampling = pem_bias_sampling = pem_variance_sampling = CSRM_pem_mse_sampling = 
    CSRM_pem_bias_sampling = 
    CSRM_pem_variance_sampling = P_sampling =
    vector("list", b_samp)
  
  
  for (b in 1:b_samp) { # sampling variability
    
    print(b)
    
    # draw an SRS without replacement
    S = sample(N_large, n_large)
    samp = sim_register[S,]
    
    # train an imputation model
    imp_model = glm(formula_imp, data = samp, family='binomial') # parameter variability
    
    # imputation
    predicted_probs = predict(imp_model, newdata = sim_register, type = "response")
    sim_register$onderwijs_imp = factor(sapply(predicted_probs, function(p) rbinom(1, 1, p)),
                                        labels = c("Geen wo", "Wo")) # imputation variability
    sim_register_b = sim_register
    sim_register_b$pp = predicted_probs
    
    
    # for benchmarks
    theta_hat_sampling[[b]] = unclass(as.matrix(
      table(sim_register$geslacht, sim_register$onderwijs_imp)[, 2, drop=F]))
    difference_sampling[[b]] = theta_hat_sampling[[b]] - theta_model[[B]]
    
    
    # accuracy estimators
    vdesign_sampling[[b]] = Vdesign(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_imp,
      domain = "geslacht",
      target = "onderwijs",
      C = c,
      design="SRS"
    )
    
    CSRM_vdesign_sampling[[b]] = sqrt(vdesign_sampling[[b]])/theta_hat_sampling[[b]]
    
    gmse_results = GMSE(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_imp,
      domain = "geslacht",
      target = "onderwijs",
      C = c,
      imp_model = imp_model
    )
    
    gmse_sampling[[b]] = gmse_results$GMSE_stoc
    CSRM_gmse_sampling[[b]] = sqrt(gmse_sampling[[b]])/theta_hat_sampling[[b]]
    
    PEM_results = PEM_srs(
      sample_data = samp,
      register_data = sim_register_b,
      C = c,
      target = "onderwijs",
      domain = "geslacht",
      matrix = "domain-specific"
    )
    
    
    pem_mse_sampling[[b]] = PEM_results$MSE
    pem_bias_sampling[[b]] = PEM_results$bias
    pem_variance_sampling[[b]] = PEM_results$variance
    P_sampling[[b]] = cbind(PEM_results$p00, PEM_results$p11)
    
    CSRM_pem_mse_sampling[[b]] = sqrt(pem_mse_sampling[[b]])/theta_hat_sampling[[b]]
    CSRM_pem_bias_sampling[[b]] = pem_bias_sampling[[b]]/theta_hat_sampling[[b]]
    CSRM_pem_variance_sampling[[b]] = sqrt(pem_variance_sampling[[b]])/theta_hat_sampling[[b]]
    
  }
  
  # benchmark across b samples
  theta_hat_model[[B]] = apply(simplify2array(theta_hat_sampling), c(1,2), mean)
  difference_sampling = simplify2array(difference_sampling)
  bias_model[[B]] = apply(difference_sampling, c(1, 2), mean)
  CSRM_bias_model[[B]] = bias_model[[B]]/theta_model[[B]]
  variance_model[[B]] = apply(difference_sampling, c(1,2), var)
  CSRM_variance_model[[B]] = sqrt(variance_model[[B]])/theta_model[[B]]
  MSE_model[[B]] = bias_model[[B]]**2 + variance_model[[B]]
  CSRM_MSE_model[[B]] = sqrt(MSE_model[[B]])/theta_model[[B]]
  
  # accuracy estimators across b samples
  vdesign_model[[B]] = apply(simplify2array(vdesign_sampling), c(1, 2), mean)
  CSRM_vdesign_model[[B]] = apply(simplify2array(CSRM_vdesign_sampling), c(1, 2), mean)
  
  gmse_model[[B]] = apply(simplify2array(gmse_sampling), c(1,2), mean)
  CSRM_gmse_model[[B]] = apply(simplify2array(CSRM_gmse_sampling), c(1,2), mean)
  
  
  pem_mse_model[[B]] = apply(simplify2array(pem_mse_sampling), c(1,2), mean)
  CSRM_pem_mse_model[[B]] = apply(simplify2array(CSRM_pem_mse_sampling), c(1,2), mean)
  pem_bias_model[[B]] = apply(simplify2array(pem_bias_sampling), c(1,2), mean)
  pem_variance_model[[B]] = apply(simplify2array(pem_variance_sampling), c(1,2), mean)
  CSRM_pem_bias_model[[B]] = apply(simplify2array(CSRM_pem_bias_sampling), c(1,2), mean)
  CSRM_pem_variance_model[[B]] = apply(simplify2array(CSRM_pem_variance_sampling), c(1,2), mean)
  
  P_model[[B]] = apply(simplify2array(P_sampling), c(1,2), mean)
  
}

results = results_dfs(results_model_names, results_sampling_names)

# write results to disk
save(P_model, P_sampling, file = "sup1_overparameterisation_P.RData")
write.csv(results$results_sampling_df, "sup1_overparameterisation_sampling_distribution.csv")
write.csv(results$results_model_df, "sup1_overparameterisation_model_distribution.csv")
write.csv(results$results_joint_df, "sup1_overparameterisation_joint_distribution.csv")


##########################################################################################

                             # Non-ignorable sampling condition #


##########################################################################################
set.seed(4789)

Nh = summary(pop_balanced$inkomen3) # oversample people with middle and higher incomes (minority groups)

nh = c(round(0.2*n_large), round(0.3*n_large), round(0.5*n_large)) # divide n unequally between the strata
pi = nh/Nh
pop_balanced$pi = ifelse(pop_balanced$inkomen3 == "High", pi[3],
                         ifelse(pop_balanced$inkomen3 == "Middle", pi[2],
                                pi[1]))


# define the items to store for each finite population; so storing the model distribution
bias_model = variance_model = theta_model = theta_hat_model = MSE_model = CSRM_MSE_model = 
  vdesign_model = CSRM_vdesign_model = gmse_model = CSRM_gmse_model =
  pem_mse_model = CSRM_pem_mse_model = pem_bias_model = pem_variance_model = CSRM_pem_bias_model = CSRM_pem_variance_model = 
  P_model = CSRM_bias_model = CSRM_variance_model =
  vector("list", B_pop)

for (B in 1:B_pop) {
  
  print(paste(B, "Non-ignorable"))
  
  # generate y in register; model variability
  sim_register = pop_balanced
  p.true = predict(superpopulation_1, newdata = sim_register, type="response")
  sim_register$onderwijs = factor(sapply(p.true, function(p) rbinom(1, 1, p)),
                                  labels = c("Geen wo", "Wo"))
  
  theta_model[[B]] = unclass(as.matrix(
    table(sim_register$geslacht, sim_register$onderwijs)[, 2, drop=F])) # matrix object for later data manipulation
  
  # store results from each sample drawn from the finite population; sampling distribution
  difference_sampling = theta_hat_sampling = vdesign_sampling = CSRM_vdesign_sampling =
    gmse_sampling = CSRM_gmse_sampling =
    pem_mse_sampling = pem_bias_sampling = pem_variance_sampling = CSRM_pem_mse_sampling = 
    CSRM_pem_bias_sampling = 
    CSRM_pem_variance_sampling = P_sampling =
    vector("list", b_samp)
  
  
  for (b in 1:b_samp) { # sampling variability
    
    print(b)
    
    # draw an SRS without replacement in each strata
    samp = vector("list", length(Nh))
    for (stratum in seq_along(Nh)) {
      strat_levels = levels(sim_register$inkomen3)
      sim_subset = sim_register[sim_register$inkomen3 == strat_levels[stratum], ]
      S = sample(Nh[stratum], nh[stratum])
      samp[[stratum]] = sim_subset[S, ]
    }
    
    samp = bind_rows(samp)
    
    ### unweighted modelling or assuming non-informativeness; input for GMSE
    
    # train an imputation model
    imp_model = glm(formula_1, data = samp, family='binomial') # input for GMSE
    
    # imputation
    predicted_probs = predict(imp_model, newdata = sim_register, type = "response")
    sim_register$onderwijs_imp = factor(sapply(predicted_probs, function(p) rbinom(1, 1, p)),
                                        labels = c("Geen wo", "Wo")) 
    sim_register_b = sim_register
    sim_register_b$pp = predicted_probs
    
    gmse_results = GMSE(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_1,
      domain = "geslacht",
      target = "onderwijs",
      C = c,
      imp_model = imp_model
    )
    
    
    ### weighted modelling
    
    # train an imputation model
    imp_model = glm(formula_1, data = samp, family='binomial', weights = 1/samp$pi) # input for Vdesign and MSEboot
    
    # imputation
    predicted_probs = predict(imp_model, newdata = sim_register, type = "response")
    sim_register$onderwijs_imp = factor(sapply(predicted_probs, function(p) rbinom(1, 1, p)),
                                        labels = c("Geen wo", "Wo")) 
    sim_register_b = sim_register
    sim_register_b$pp = predicted_probs
    
    # for benchmarks
    theta_hat_sampling[[b]] = unclass(as.matrix(
      table(sim_register$geslacht, sim_register$onderwijs_imp)[, 2, drop=F]))
    difference_sampling[[b]] = theta_hat_sampling[[b]] - theta_model[[B]]
    
    
    # accuracy estimators
    
    gmse_sampling[[b]] = gmse_results$GMSE_stoc
    CSRM_gmse_sampling[[b]] = sqrt(gmse_sampling[[b]])/theta_hat_sampling[[b]]
    
    
    samp$strat_var = samp$inkomen3
    
    vdesign_sampling[[b]] = Vdesign(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_1,
      domain = "geslacht",
      target = "onderwijs",
      C = c,
      design="SRS"
    )
    
    CSRM_vdesign_sampling[[b]] = sqrt(vdesign_sampling[[b]])/theta_hat_sampling[[b]]
    
    PEM_results = PEM_strat_binom(
      sample_data = samp,
      register_data = sim_register_b,
      C = c,
      target = "onderwijs",
      domain = "geslacht",
      strat_var = "inkomen3"
    )
    
    
    pem_mse_sampling[[b]] = PEM_results$MSE
    pem_bias_sampling[[b]] = PEM_results$bias
    pem_variance_sampling[[b]] = PEM_results$variance
    P_sampling[[b]] = cbind(PEM_results$p00, PEM_results$p11)
    
    CSRM_pem_mse_sampling[[b]] = sqrt(pem_mse_sampling[[b]])/theta_hat_sampling[[b]]
    CSRM_pem_bias_sampling[[b]] = pem_bias_sampling[[b]]/theta_hat_sampling[[b]]
    CSRM_pem_variance_sampling[[b]] = sqrt(pem_variance_sampling[[b]])/theta_hat_sampling[[b]]
    
  }
  
  # benchmark across b samples
  theta_hat_model[[B]] = apply(simplify2array(theta_hat_sampling), c(1,2), mean)
  difference_sampling = simplify2array(difference_sampling)
  bias_model[[B]] = apply(difference_sampling, c(1, 2), mean)
  CSRM_bias_model[[B]] = bias_model[[B]]/theta_model[[B]]
  variance_model[[B]] = apply(difference_sampling, c(1,2), var)
  CSRM_variance_model[[B]] = sqrt(variance_model[[B]])/theta_model[[B]]
  MSE_model[[B]] = bias_model[[B]]**2 + variance_model[[B]]
  CSRM_MSE_model[[B]] = sqrt(MSE_model[[B]])/theta_model[[B]]
  
  # accuracy estimators across b samples
  vdesign_model[[B]] = apply(simplify2array(vdesign_sampling), c(1, 2), mean)
  CSRM_vdesign_model[[B]] = apply(simplify2array(CSRM_vdesign_sampling), c(1, 2), mean)
  
  gmse_model[[B]] = apply(simplify2array(gmse_sampling), c(1,2), mean)
  CSRM_gmse_model[[B]] = apply(simplify2array(CSRM_gmse_sampling), c(1,2), mean)
  
  
  pem_mse_model[[B]] = apply(simplify2array(pem_mse_sampling), c(1,2), mean)
  CSRM_pem_mse_model[[B]] = apply(simplify2array(CSRM_pem_mse_sampling), c(1,2), mean)
  pem_bias_model[[B]] = apply(simplify2array(pem_bias_sampling), c(1,2), mean)
  pem_variance_model[[B]] = apply(simplify2array(pem_variance_sampling), c(1,2), mean)
  CSRM_pem_bias_model[[B]] = apply(simplify2array(CSRM_pem_bias_sampling), c(1,2), mean)
  CSRM_pem_variance_model[[B]] = apply(simplify2array(CSRM_pem_variance_sampling), c(1,2), mean)
  
  P_model[[B]] = apply(simplify2array(P_sampling), c(1,2), mean)
  
}

results = results_dfs(results_model_names, results_sampling_names)


# write results to disk
save(P_model, P_sampling, file = "sup1_nonignorable_P.RData")
write.csv(results$results_sampling_df, "sup1_nonignorable_sampling_distribution.csv")
write.csv(results$results_model_df, "sup1_nonignorable_model_distribution.csv")
write.csv(results$results_joint_df, "sup1_nonignorable_joint_distribution.csv")


##########################################################################################

                             # Design bias condition #


##########################################################################################
set.seed(4789)

pop_balanced$pi = n_large/N_large

# independently generate domain variable
pop_balanced$domain_random = sample(c("A", "B"), nrow(pop_balanced), replace = TRUE)
pop_balanced$domain_random = as.factor(pop_balanced$domain_random)


# define the items to store for each finite population; so storing the model distribution
bias_model = variance_model = theta_model = theta_hat_model = MSE_model = CSRM_MSE_model = 
  vdesign_model = CSRM_vdesign_model = gmse_model = CSRM_gmse_model =
  pem_mse_model = CSRM_pem_mse_model = pem_bias_model = pem_variance_model = CSRM_pem_bias_model = CSRM_pem_variance_model = 
  P_model = CSRM_bias_model = CSRM_variance_model =
  vector("list", B_pop)

for (B in 1:B_pop) {
  
  print(paste(B, "Design bias"))
  
  # generate y in register; model variability
  sim_register = pop_balanced
  p.true = predict(superpopulation_1, newdata = sim_register, type="response")
  sim_register$onderwijs = factor(sapply(p.true, function(p) rbinom(1, 1, p)),
                                  labels = c("Geen wo", "Wo"))
  
  theta_model[[B]] = unclass(as.matrix(
    table(sim_register$domain_random, sim_register$onderwijs)[, 2, drop=F])) # matrix object for later data manipulation
  
  # store results from each sample drawn from the finite population; sampling distribution
  difference_sampling = theta_hat_sampling = vdesign_sampling = CSRM_vdesign_sampling =
    gmse_sampling = CSRM_gmse_sampling =
    pem_mse_sampling = pem_bias_sampling = pem_variance_sampling = CSRM_pem_mse_sampling = 
    CSRM_pem_bias_sampling = 
    CSRM_pem_variance_sampling = P_sampling =
    vector("list", b_samp)
  
  
  for (b in 1:b_samp) { # sampling variability
    
    print(b)
    
    # draw an SRS without replacement
    S = sample(N_large, n_large)
    samp = sim_register[S,]
    
    # train an imputation model
    imp_model = glm(formula_1, data = samp, family='binomial') # parameter variability
    
    # imputation
    predicted_probs = predict(imp_model, newdata = sim_register, type = "response")
    sim_register$onderwijs_imp = factor(sapply(predicted_probs, function(p) rbinom(1, 1, p)),
                                        labels = c("Geen wo", "Wo")) # imputation variability
    sim_register_b = sim_register
    sim_register_b$pp = predicted_probs
    
    
    # for benchmarks
    theta_hat_sampling[[b]] = unclass(as.matrix(
      table(sim_register$domain_random, sim_register$onderwijs_imp)[, 2, drop=F]))
    difference_sampling[[b]] = theta_hat_sampling[[b]] - theta_model[[B]]
    
    
    # accuracy estimators
    vdesign_sampling[[b]] = Vdesign(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_1,
      domain = "domain_random",
      target = "onderwijs",
      C = c,
      design="SRS"
    )
    
    CSRM_vdesign_sampling[[b]] = sqrt(vdesign_sampling[[b]])/theta_hat_sampling[[b]]
    
    gmse_results = GMSE(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_1,
      domain = "domain_random",
      target = "onderwijs",
      C = c,
      imp_model = imp_model
    )
    
    gmse_sampling[[b]] = gmse_results$GMSE_stoc
    CSRM_gmse_sampling[[b]] = sqrt(gmse_sampling[[b]])/theta_hat_sampling[[b]]
    
    PEM_results = PEM_srs(
      sample_data = samp,
      register_data = sim_register_b,
      C = c,
      target = "onderwijs",
      domain = "domain_random",
      matrix = "non-specific"
    )
    
    
    pem_mse_sampling[[b]] = PEM_results$MSE
    pem_bias_sampling[[b]] = PEM_results$bias
    pem_variance_sampling[[b]] = PEM_results$variance
    P_sampling[[b]] = cbind(PEM_results$p00, PEM_results$p11)
    
    CSRM_pem_mse_sampling[[b]] = sqrt(pem_mse_sampling[[b]])/theta_hat_sampling[[b]]
    CSRM_pem_bias_sampling[[b]] = pem_bias_sampling[[b]]/theta_hat_sampling[[b]]
    CSRM_pem_variance_sampling[[b]] = sqrt(pem_variance_sampling[[b]])/theta_hat_sampling[[b]]
    
  }
  
  # benchmark across b samples
  theta_hat_model[[B]] = apply(simplify2array(theta_hat_sampling), c(1,2), mean)
  difference_sampling = simplify2array(difference_sampling)
  bias_model[[B]] = apply(difference_sampling, c(1, 2), mean)
  CSRM_bias_model[[B]] = bias_model[[B]]/theta_model[[B]]
  variance_model[[B]] = apply(difference_sampling, c(1,2), var)
  CSRM_variance_model[[B]] = sqrt(variance_model[[B]])/theta_model[[B]]
  MSE_model[[B]] = bias_model[[B]]**2 + variance_model[[B]]
  CSRM_MSE_model[[B]] = sqrt(MSE_model[[B]])/theta_model[[B]]
  
  # accuracy estimators across b samples
  vdesign_model[[B]] = apply(simplify2array(vdesign_sampling), c(1, 2), mean)
  CSRM_vdesign_model[[B]] = apply(simplify2array(CSRM_vdesign_sampling), c(1, 2), mean)
  
  gmse_model[[B]] = apply(simplify2array(gmse_sampling), c(1,2), mean)
  CSRM_gmse_model[[B]] = apply(simplify2array(CSRM_gmse_sampling), c(1,2), mean)
  
  
  pem_mse_model[[B]] = apply(simplify2array(pem_mse_sampling), c(1,2), mean)
  CSRM_pem_mse_model[[B]] = apply(simplify2array(CSRM_pem_mse_sampling), c(1,2), mean)
  pem_bias_model[[B]] = apply(simplify2array(pem_bias_sampling), c(1,2), mean)
  pem_variance_model[[B]] = apply(simplify2array(pem_variance_sampling), c(1,2), mean)
  CSRM_pem_bias_model[[B]] = apply(simplify2array(CSRM_pem_bias_sampling), c(1,2), mean)
  CSRM_pem_variance_model[[B]] = apply(simplify2array(CSRM_pem_variance_sampling), c(1,2), mean)
  
  P_model[[B]] = apply(simplify2array(P_sampling), c(1,2), mean)
  
}

results = results_dfs(results_model_names, results_sampling_names)


# write results to disk
save(P_model, P_sampling, file = "sup1_designbias_P.RData")
write.csv(results$results_sampling_df, "sup1_designbias_sampling_distribution.csv")
write.csv(results$results_model_df, "sup1_designbias_model_distribution.csv")
write.csv(results$results_joint_df, "sup1_designbias_joint_distribution.csv")


##########################################################################################

                                # Model bias (small) condition #


##########################################################################################
set.seed(4789)

covs_imp = c("leeftijd", "inkomen3")
formula_imp = as.formula(paste("onderwijs ~", paste(covs_imp, collapse = " + "))) 


# define the items to store for each finite population; so storing the model distribution
bias_model = variance_model = theta_model = theta_hat_model = MSE_model = CSRM_MSE_model = 
  vdesign_model = CSRM_vdesign_model = gmse_model = CSRM_gmse_model =
  pem_mse_model = CSRM_pem_mse_model = pem_bias_model = pem_variance_model = CSRM_pem_bias_model = CSRM_pem_variance_model = 
  P_model = CSRM_bias_model = CSRM_variance_model =
  vector("list", B_pop)

for (B in 1:B_pop) {
  
  print(paste(B, "Model bias (small)"))
  
  # generate y in register; model variability
  sim_register = pop_balanced
  p.true = predict(superpopulation_1, newdata = sim_register, type="response")
  sim_register$onderwijs = factor(sapply(p.true, function(p) rbinom(1, 1, p)),
                                  labels = c("Geen wo", "Wo"))
  
  theta_model[[B]] = unclass(as.matrix(
    table(sim_register$geslacht, sim_register$onderwijs)[, 2, drop=F])) # matrix object for later data manipulation
  
  # store results from each sample drawn from the finite population; sampling distribution
  difference_sampling = theta_hat_sampling = vdesign_sampling = CSRM_vdesign_sampling =
    gmse_sampling = CSRM_gmse_sampling =
    pem_mse_sampling = pem_bias_sampling = pem_variance_sampling = CSRM_pem_mse_sampling = 
    CSRM_pem_bias_sampling = 
    CSRM_pem_variance_sampling = P_sampling =
    vector("list", b_samp)
  
  
  for (b in 1:b_samp) { # sampling variability
    
    print(b)
    
    # draw an SRS without replacement
    S = sample(N_large, n_large)
    samp = sim_register[S,]
    
    # train an imputation model
    imp_model = glm(formula_imp, data = samp, family='binomial') # parameter variability
    
    # imputation
    predicted_probs = predict(imp_model, newdata = sim_register, type = "response")
    sim_register$onderwijs_imp = factor(sapply(predicted_probs, function(p) rbinom(1, 1, p)),
                                        labels = c("Geen wo", "Wo")) # imputation variability
    sim_register_b = sim_register
    sim_register_b$pp = predicted_probs
    
    
    # for benchmarks
    theta_hat_sampling[[b]] = unclass(as.matrix(
      table(sim_register$geslacht, sim_register$onderwijs_imp)[, 2, drop=F]))
    difference_sampling[[b]] = theta_hat_sampling[[b]] - theta_model[[B]]
    
    
    # accuracy estimators
    vdesign_sampling[[b]] = Vdesign(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_imp,
      domain = "geslacht",
      target = "onderwijs",
      C = c,
      design="SRS"
    )
    
    CSRM_vdesign_sampling[[b]] = sqrt(vdesign_sampling[[b]])/theta_hat_sampling[[b]]
    
    gmse_results = GMSE(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_imp,
      domain = "geslacht",
      target = "onderwijs",
      C = c,
      imp_model = imp_model
    )
    
    gmse_sampling[[b]] = gmse_results$GMSE_stoc
    CSRM_gmse_sampling[[b]] = sqrt(gmse_sampling[[b]])/theta_hat_sampling[[b]]
    
    PEM_results = PEM_srs(
      sample_data = samp,
      register_data = sim_register_b,
      C = c,
      target = "onderwijs",
      domain = "geslacht",
      matrix = "non-specific"
    )
    
    
    pem_mse_sampling[[b]] = PEM_results$MSE
    pem_bias_sampling[[b]] = PEM_results$bias
    pem_variance_sampling[[b]] = PEM_results$variance
    P_sampling[[b]] = cbind(PEM_results$p00, PEM_results$p11)
    
    CSRM_pem_mse_sampling[[b]] = sqrt(pem_mse_sampling[[b]])/theta_hat_sampling[[b]]
    CSRM_pem_bias_sampling[[b]] = pem_bias_sampling[[b]]/theta_hat_sampling[[b]]
    CSRM_pem_variance_sampling[[b]] = sqrt(pem_variance_sampling[[b]])/theta_hat_sampling[[b]]
    
  }
  
  # benchmark across b samples
  theta_hat_model[[B]] = apply(simplify2array(theta_hat_sampling), c(1,2), mean)
  difference_sampling = simplify2array(difference_sampling)
  bias_model[[B]] = apply(difference_sampling, c(1, 2), mean)
  CSRM_bias_model[[B]] = bias_model[[B]]/theta_model[[B]]
  variance_model[[B]] = apply(difference_sampling, c(1,2), var)
  CSRM_variance_model[[B]] = sqrt(variance_model[[B]])/theta_model[[B]]
  MSE_model[[B]] = bias_model[[B]]**2 + variance_model[[B]]
  CSRM_MSE_model[[B]] = sqrt(MSE_model[[B]])/theta_model[[B]]
  
  # accuracy estimators across b samples
  vdesign_model[[B]] = apply(simplify2array(vdesign_sampling), c(1, 2), mean)
  CSRM_vdesign_model[[B]] = apply(simplify2array(CSRM_vdesign_sampling), c(1, 2), mean)
  
  gmse_model[[B]] = apply(simplify2array(gmse_sampling), c(1,2), mean)
  CSRM_gmse_model[[B]] = apply(simplify2array(CSRM_gmse_sampling), c(1,2), mean)
  
  
  pem_mse_model[[B]] = apply(simplify2array(pem_mse_sampling), c(1,2), mean)
  CSRM_pem_mse_model[[B]] = apply(simplify2array(CSRM_pem_mse_sampling), c(1,2), mean)
  pem_bias_model[[B]] = apply(simplify2array(pem_bias_sampling), c(1,2), mean)
  pem_variance_model[[B]] = apply(simplify2array(pem_variance_sampling), c(1,2), mean)
  CSRM_pem_bias_model[[B]] = apply(simplify2array(CSRM_pem_bias_sampling), c(1,2), mean)
  CSRM_pem_variance_model[[B]] = apply(simplify2array(CSRM_pem_variance_sampling), c(1,2), mean)
  
  P_model[[B]] = apply(simplify2array(P_sampling), c(1,2), mean)
  
}

results = results_dfs(results_model_names, results_sampling_names)

# write results to disk
save(P_model, P_sampling, file = "sup1_modelbias_small_P.RData")
write.csv(results$results_sampling_df, "sup1_modelbias_small_sampling_distribution.csv")
write.csv(results$results_model_df, "sup1_modelbias_small_model_distribution.csv")
write.csv(results$results_joint_df, "sup1_modelbias_small_joint_distribution.csv")

##########################################################################################

                              # Model bias (large) condition #


##########################################################################################
set.seed(4789)

assign_gender <- function(gen) {
  if (gen == "Wo") {
    sample(c("Man", "Vrouw"), size = 1, prob = c(0.6, 0.4))  
  } else {
    sample(c("Man", "Vrouw"), size = 1, prob = c(0.4, 0.6))  
  }
}

# Reassign income3 based on education
pop_balanced$geslacht_syn <- sapply(pop_balanced$onderwijs, assign_gender)
pop_balanced$geslacht_syn = as.factor(pop_balanced$geslacht_syn)

covs_syn = c("geslacht_syn", "leeftijd", "inkomen3") # geslacht will be the domain variable
formula_syn = as.formula(paste("onderwijs ~", paste(covs_syn, collapse = " + "))) 
superpopulation_syn = glm(formula_syn, data = pop_balanced, family='binomial')

# define the items to store for each finite population; so storing the model distribution
bias_model = variance_model = theta_model = theta_hat_model = MSE_model = CSRM_MSE_model = 
  vdesign_model = CSRM_vdesign_model = gmse_model = CSRM_gmse_model =
  pem_mse_model = CSRM_pem_mse_model = pem_bias_model = pem_variance_model = CSRM_pem_bias_model = CSRM_pem_variance_model = 
  P_model = CSRM_bias_model = CSRM_variance_model =
  vector("list", B_pop)

for (B in 1:B_pop) {
  
  print(paste(B, "Model bias (large)"))
  
  # generate y in register; model variability
  sim_register = pop_balanced
  p.true = predict(superpopulation_syn, newdata = sim_register, type="response")
  sim_register$onderwijs = factor(sapply(p.true, function(p) rbinom(1, 1, p)),
                                  labels = c("Geen wo", "Wo"))
  
  theta_model[[B]] = unclass(as.matrix(
    table(sim_register$geslacht_syn, sim_register$onderwijs)[, 2, drop=F])) # matrix object for later data manipulation
  
  # store results from each sample drawn from the finite population; sampling distribution
  difference_sampling = theta_hat_sampling = vdesign_sampling = CSRM_vdesign_sampling =
    gmse_sampling = CSRM_gmse_sampling =
    pem_mse_sampling = pem_bias_sampling = pem_variance_sampling = CSRM_pem_mse_sampling = 
    CSRM_pem_bias_sampling = 
    CSRM_pem_variance_sampling = P_sampling =
    vector("list", b_samp)
  
  
  for (b in 1:b_samp) { # sampling variability
    
    print(b)
    
    # draw an SRS without replacement
    S = sample(N_large, n_large)
    samp = sim_register[S,]
    
    # train an imputation model
    imp_model = glm(formula_imp, data = samp, family='binomial') # parameter variability
    
    # imputation
    predicted_probs = predict(imp_model, newdata = sim_register, type = "response")
    sim_register$onderwijs_imp = factor(sapply(predicted_probs, function(p) rbinom(1, 1, p)),
                                        labels = c("Geen wo", "Wo")) # imputation variability
    sim_register_b = sim_register
    sim_register_b$pp = predicted_probs
    
    
    # for benchmarks
    theta_hat_sampling[[b]] = unclass(as.matrix(
      table(sim_register$geslacht_syn, sim_register$onderwijs_imp)[, 2, drop=F]))
    difference_sampling[[b]] = theta_hat_sampling[[b]] - theta_model[[B]]
    
    
    # accuracy estimators
    vdesign_sampling[[b]] = Vdesign(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_imp,
      domain = "geslacht_syn",
      target = "onderwijs",
      C = c,
      design="SRS"
    )
    
    CSRM_vdesign_sampling[[b]] = sqrt(vdesign_sampling[[b]])/theta_hat_sampling[[b]]
    
    gmse_results = GMSE(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_imp,
      domain = "geslacht_syn",
      target = "onderwijs",
      C = c,
      imp_model = imp_model
    )
    
    gmse_sampling[[b]] = gmse_results$GMSE_stoc
    CSRM_gmse_sampling[[b]] = sqrt(gmse_sampling[[b]])/theta_hat_sampling[[b]]
    
    PEM_results = PEM_srs(
      sample_data = samp,
      register_data = sim_register_b,
      C = c,
      target = "onderwijs",
      domain = "geslacht_syn",
      matrix = "non-specific"
    )
    
    
    pem_mse_sampling[[b]] = PEM_results$MSE
    pem_bias_sampling[[b]] = PEM_results$bias
    pem_variance_sampling[[b]] = PEM_results$variance
    P_sampling[[b]] = cbind(PEM_results$p00, PEM_results$p11)
    
    CSRM_pem_mse_sampling[[b]] = sqrt(pem_mse_sampling[[b]])/theta_hat_sampling[[b]]
    CSRM_pem_bias_sampling[[b]] = pem_bias_sampling[[b]]/theta_hat_sampling[[b]]
    CSRM_pem_variance_sampling[[b]] = sqrt(pem_variance_sampling[[b]])/theta_hat_sampling[[b]]
    
  }
  
  # benchmark across b samples
  theta_hat_model[[B]] = apply(simplify2array(theta_hat_sampling), c(1,2), mean)
  difference_sampling = simplify2array(difference_sampling)
  bias_model[[B]] = apply(difference_sampling, c(1, 2), mean)
  CSRM_bias_model[[B]] = bias_model[[B]]/theta_model[[B]]
  variance_model[[B]] = apply(difference_sampling, c(1,2), var)
  CSRM_variance_model[[B]] = sqrt(variance_model[[B]])/theta_model[[B]]
  MSE_model[[B]] = bias_model[[B]]**2 + variance_model[[B]]
  CSRM_MSE_model[[B]] = sqrt(MSE_model[[B]])/theta_model[[B]]
  
  # accuracy estimators across b samples
  vdesign_model[[B]] = apply(simplify2array(vdesign_sampling), c(1, 2), mean)
  CSRM_vdesign_model[[B]] = apply(simplify2array(CSRM_vdesign_sampling), c(1, 2), mean)
  
  gmse_model[[B]] = apply(simplify2array(gmse_sampling), c(1,2), mean)
  CSRM_gmse_model[[B]] = apply(simplify2array(CSRM_gmse_sampling), c(1,2), mean)
  
  
  pem_mse_model[[B]] = apply(simplify2array(pem_mse_sampling), c(1,2), mean)
  CSRM_pem_mse_model[[B]] = apply(simplify2array(CSRM_pem_mse_sampling), c(1,2), mean)
  pem_bias_model[[B]] = apply(simplify2array(pem_bias_sampling), c(1,2), mean)
  pem_variance_model[[B]] = apply(simplify2array(pem_variance_sampling), c(1,2), mean)
  CSRM_pem_bias_model[[B]] = apply(simplify2array(CSRM_pem_bias_sampling), c(1,2), mean)
  CSRM_pem_variance_model[[B]] = apply(simplify2array(CSRM_pem_variance_sampling), c(1,2), mean)
  
  P_model[[B]] = apply(simplify2array(P_sampling), c(1,2), mean)
  
}

results = results_dfs(results_model_names, results_sampling_names)

# write results to disk
save(P_model, P_sampling, file = "sup1_modelbias_large_P.RData")
write.csv(results$results_sampling_df, "sup1_modelbias_large_sampling_distribution.csv")
write.csv(results$results_model_df, "sup1_modelbias_large_model_distribution.csv")
write.csv(results$results_joint_df, "sup1_modelbias_large_joint_distribution.csv")

##########################################################################################

                      # Model bias (large) X small sample condition #


##########################################################################################
set.seed(4789)
pop_balanced$pi = n_small/N_large

# define the items to store for each finite population; so storing the model distribution
bias_model = variance_model = theta_model = theta_hat_model = MSE_model = CSRM_MSE_model = 
  vdesign_model = CSRM_vdesign_model = gmse_model = CSRM_gmse_model =
  pem_mse_model = CSRM_pem_mse_model = pem_bias_model = pem_variance_model = CSRM_pem_bias_model = CSRM_pem_variance_model = 
  P_model = CSRM_bias_model = CSRM_variance_model =
  vector("list", B_pop)

for (B in 1:B_pop) {
  
  print(paste(B, "Model bias (large) X small sample"))
  
  # generate y in register; model variability
  sim_register = pop_balanced
  p.true = predict(superpopulation_syn, newdata = sim_register, type="response")
  sim_register$onderwijs = factor(sapply(p.true, function(p) rbinom(1, 1, p)),
                                  labels = c("Geen wo", "Wo"))
  
  theta_model[[B]] = unclass(as.matrix(
    table(sim_register$geslacht_syn, sim_register$onderwijs)[, 2, drop=F])) # matrix object for later data manipulation
  
  # store results from each sample drawn from the finite population; sampling distribution
  difference_sampling = theta_hat_sampling = vdesign_sampling = CSRM_vdesign_sampling =
    gmse_sampling = CSRM_gmse_sampling =
    pem_mse_sampling = pem_bias_sampling = pem_variance_sampling = CSRM_pem_mse_sampling = 
    CSRM_pem_bias_sampling = 
    CSRM_pem_variance_sampling = P_sampling =
    vector("list", b_samp)
  
  
  for (b in 1:b_samp) { # sampling variability
    
    print(b)
    
    # draw an SRS without replacement
    S = sample(N_large, n_small)
    samp = sim_register[S,]
    
    # train an imputation model
    imp_model = glm(formula_imp, data = samp, family='binomial') # parameter variability
    
    # imputation
    predicted_probs = predict(imp_model, newdata = sim_register, type = "response")
    sim_register$onderwijs_imp = factor(sapply(predicted_probs, function(p) rbinom(1, 1, p)),
                                        labels = c("Geen wo", "Wo")) # imputation variability
    sim_register_b = sim_register
    sim_register_b$pp = predicted_probs
    
    
    # for benchmarks
    theta_hat_sampling[[b]] = unclass(as.matrix(
      table(sim_register$geslacht_syn, sim_register$onderwijs_imp)[, 2, drop=F]))
    difference_sampling[[b]] = theta_hat_sampling[[b]] - theta_model[[B]]
    
    
    # accuracy estimators
    vdesign_sampling[[b]] = Vdesign(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_imp,
      domain = "geslacht_syn",
      target = "onderwijs",
      C = c,
      design="SRS"
    )
    
    CSRM_vdesign_sampling[[b]] = sqrt(vdesign_sampling[[b]])/theta_hat_sampling[[b]]
    
    gmse_results = GMSE(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_imp,
      domain = "geslacht_syn",
      target = "onderwijs",
      C = c,
      imp_model = imp_model
    )
    
    gmse_sampling[[b]] = gmse_results$GMSE_stoc
    CSRM_gmse_sampling[[b]] = sqrt(gmse_sampling[[b]])/theta_hat_sampling[[b]]
    
    PEM_results = PEM_srs(
      sample_data = samp,
      register_data = sim_register_b,
      C = c,
      target = "onderwijs",
      domain = "geslacht_syn",
      matrix = "non-specific"
    )
    
    
    pem_mse_sampling[[b]] = PEM_results$MSE
    pem_bias_sampling[[b]] = PEM_results$bias
    pem_variance_sampling[[b]] = PEM_results$variance
    P_sampling[[b]] = cbind(PEM_results$p00, PEM_results$p11)
    
    CSRM_pem_mse_sampling[[b]] = sqrt(pem_mse_sampling[[b]])/theta_hat_sampling[[b]]
    CSRM_pem_bias_sampling[[b]] = pem_bias_sampling[[b]]/theta_hat_sampling[[b]]
    CSRM_pem_variance_sampling[[b]] = sqrt(pem_variance_sampling[[b]])/theta_hat_sampling[[b]]
    
  }
  
  # benchmark across b samples
  theta_hat_model[[B]] = apply(simplify2array(theta_hat_sampling), c(1,2), mean)
  difference_sampling = simplify2array(difference_sampling)
  bias_model[[B]] = apply(difference_sampling, c(1, 2), mean)
  CSRM_bias_model[[B]] = bias_model[[B]]/theta_model[[B]]
  variance_model[[B]] = apply(difference_sampling, c(1,2), var)
  CSRM_variance_model[[B]] = sqrt(variance_model[[B]])/theta_model[[B]]
  MSE_model[[B]] = bias_model[[B]]**2 + variance_model[[B]]
  CSRM_MSE_model[[B]] = sqrt(MSE_model[[B]])/theta_model[[B]]
  
  # accuracy estimators across b samples
  vdesign_model[[B]] = apply(simplify2array(vdesign_sampling), c(1, 2), mean)
  CSRM_vdesign_model[[B]] = apply(simplify2array(CSRM_vdesign_sampling), c(1, 2), mean)
  
  gmse_model[[B]] = apply(simplify2array(gmse_sampling), c(1,2), mean)
  CSRM_gmse_model[[B]] = apply(simplify2array(CSRM_gmse_sampling), c(1,2), mean)
  
  
  pem_mse_model[[B]] = apply(simplify2array(pem_mse_sampling), c(1,2), mean)
  CSRM_pem_mse_model[[B]] = apply(simplify2array(CSRM_pem_mse_sampling), c(1,2), mean)
  pem_bias_model[[B]] = apply(simplify2array(pem_bias_sampling), c(1,2), mean)
  pem_variance_model[[B]] = apply(simplify2array(pem_variance_sampling), c(1,2), mean)
  CSRM_pem_bias_model[[B]] = apply(simplify2array(CSRM_pem_bias_sampling), c(1,2), mean)
  CSRM_pem_variance_model[[B]] = apply(simplify2array(CSRM_pem_variance_sampling), c(1,2), mean)
  
  P_model[[B]] = apply(simplify2array(P_sampling), c(1,2), mean)
  
}

results = results_dfs(results_model_names, results_sampling_names)

# write results to disk
save(P_model, P_sampling, file = "sup1_modelbias_large_small_samp_P.RData")
write.csv(results$results_sampling_df, "sup1_modelbias_large_small_samp_sampling_distribution.csv")
write.csv(results$results_model_df, "sup1_modelbias_large_small_samp_model_distribution.csv")
write.csv(results$results_joint_df, "sup1_modelbias_large_small_samp_joint_distribution.csv")

################################################################################################

                              # Multinomial condition #

################################################################################################

set.seed(4789)
B_pop=1
b_samp = 1000

pop <- read.csv2(file = 'samplonie.csv', header = TRUE,
                 stringsAsFactors = TRUE)

pop <- pop[pop$leeftijd >= 15, ] #select people over 15


#turn income into categorical variable
pop$inkomen3 <- cut(pop$inkomen, 
                    breaks = c(0, 350, 1400, 5000), 
                    include.lowest = TRUE, 
                    labels = c("Low", "Middle", "High"))

# multinomial target variable with unbalanced classes
c = 3

pop$onderwijs <- ifelse(pop$onderwijs %in% c('Geen','Basisonderwijs'), 1,
                        ifelse(pop$onderwijs %in% c('Vmbo','Havo/vwo','Mbo'),2,3))
pop$onderwijs <- factor(pop$onderwijs, labels = c('Low','Middle','High'))

# large N 
pop <- pop[rep(1:nrow(pop), times = 150), ]
pop$id = 1:nrow(pop)

N = nrow(pop) # N = 111750
n = round(0.05*N)  # n = 5588
# non-informative sampling (SRSWOR)
pop$pi = n/N 


# create the superpopulation
# log-scale the numerical variable
pop$leeftijd = scale(log(pop$leeftijd+0.0001)) 

covs = c("geslacht", "leeftijd", "inkomen3") # geslacht will be the domain variable
formula = as.formula(paste("onderwijs ~", paste(covs, collapse = " + "))) 
superpopulation = multinom(formula, data = pop, trace = FALSE)


# define the items to store for each finite population; so storing the model distribution
bias_model = variance_model = theta_model = theta_hat_model = MSE_model = CSRM_MSE_model = 
  vdesign_model = CSRM_vdesign_model = gmse_model = CSRM_gmse_model =
  pem_mse_model = CSRM_pem_mse_model = pem_bias_model = pem_variance_model = CSRM_pem_bias_model = CSRM_pem_variance_model = 
  P_model = CSRM_bias_model = CSRM_variance_model =
  vector("list", B_pop)

for (B in 1:B_pop) {
  
  print(paste(B, "Multinomial"))
  
  # generate y in register; model variability
  sim_register = pop
  p.true = predict(superpopulation, newdata = sim_register, type="probs")
  sim_register$onderwijs = factor(apply(p.true, 1, function(m) which(rmultinom(1, 1, m) == 1)),
                                  levels = 1:ncol(p.true),
                                  labels = colnames(p.true))
  
  theta_model[[B]] = unclass(as.matrix(
    table(sim_register$geslacht, sim_register$onderwijs)))
  
  
  # store results from each sample drawn from the finite population; sampling distribution
  difference_sampling = theta_hat_sampling = vdesign_sampling = CSRM_vdesign_sampling =
    gmse_sampling = CSRM_gmse_sampling =
    pem_mse_sampling = pem_bias_sampling = pem_variance_sampling = CSRM_pem_mse_sampling = 
    CSRM_pem_bias_sampling = 
    CSRM_pem_variance_sampling = P_sampling =
    vector("list", b_samp)
  
  
  for (b in 1:b_samp) { # sampling variability
    
    print(paste("b =", b))
    
    # draw an SRS without replacement
    S = sample(N, n)
    samp = sim_register[S,]
    
    # train an imputation model
    imp_model = multinom(formula, data = samp, trace = FALSE) # parameter variability
    
    # imputation
    predicted_probs = predict(imp_model, newdata = sim_register, type = "probs")
    sim_register$onderwijs_imp = factor(apply(predicted_probs, 1, function(m) which(rmultinom(1, 1, m) == 1)),
                                        levels = 1:ncol(predicted_probs),
                                        labels = colnames(predicted_probs)) # imputation variability
    
    colnames(predicted_probs) = paste0("p", 1:ncol(predicted_probs))
    
    sim_register_b = cbind(sim_register, as.data.frame(predicted_probs))
    
    print("Created imputed values")
    
    # for benchmarks
    theta_hat_sampling[[b]] = unclass(as.matrix(
      table(sim_register$geslacht, sim_register$onderwijs_imp)))
    difference_sampling[[b]] = theta_hat_sampling[[b]] - theta_model[[B]]
    
    
    print("Computing Vdesign")
    
    # accuracy estimators
    vdesign_sampling[[b]] = Vdesign(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_1,
      domain = "geslacht",
      target = "onderwijs",
      C = c,
      design="SRS"
    )
    
    CSRM_vdesign_sampling[[b]] = sqrt(vdesign_sampling[[b]])/theta_hat_sampling[[b]]
    
    print("Computing GMSE")
    
    gmse_results = GMSE(
      sample_data = samp,
      register_data = sim_register_b,
      covariates = covs_1,
      domain = "geslacht",
      target = "onderwijs",
      C = c,
      imp_model = imp_model
    )
    
    gmse_sampling[[b]] = gmse_results$GMSE_stoc
    CSRM_gmse_sampling[[b]] = sqrt(gmse_sampling[[b]])/theta_hat_sampling[[b]]
    
    
    print("Computing PEM")
    
    PEM_results = PEM_srs(
      sample_data = samp,
      register_data = sim_register_b,
      C = c,
      target = "onderwijs",
      domain = "geslacht",
      matrix = "domain-specific"
    )
    
    
    pem_mse_sampling[[b]] = PEM_results$MSE
    pem_bias_sampling[[b]] = PEM_results$bias
    pem_variance_sampling[[b]] = PEM_results$variance
    P_sampling[[b]] = PEM_results$P
    
    CSRM_pem_mse_sampling[[b]] = sqrt(pem_mse_sampling[[b]])/theta_hat_sampling[[b]]
    CSRM_pem_bias_sampling[[b]] = pem_bias_sampling[[b]]/theta_hat_sampling[[b]]
    CSRM_pem_variance_sampling[[b]] = sqrt(pem_variance_sampling[[b]])/theta_hat_sampling[[b]]
    
  }
  
  # benchmark across b samples
  theta_hat_model[[B]] = apply(simplify2array(theta_hat_sampling), c(1,2), mean)
  difference_sampling = simplify2array(difference_sampling)
  bias_model[[B]] = apply(difference_sampling, c(1, 2), mean)
  CSRM_bias_model[[B]] = bias_model[[B]]/theta_model[[B]]
  variance_model[[B]] = apply(difference_sampling, c(1,2), var)
  CSRM_variance_model[[B]] = sqrt(variance_model[[B]])/theta_model[[B]]
  MSE_model[[B]] = bias_model[[B]]**2 + variance_model[[B]]
  CSRM_MSE_model[[B]] = sqrt(MSE_model[[B]])/theta_model[[B]]
  
  # accuracy estimators across b samples
  vdesign_model[[B]] = apply(simplify2array(vdesign_sampling), c(1, 2), mean)
  CSRM_vdesign_model[[B]] = apply(simplify2array(CSRM_vdesign_sampling), c(1, 2), mean)
  
  gmse_model[[B]] = apply(simplify2array(gmse_sampling), c(1,2), mean)
  CSRM_gmse_model[[B]] = apply(simplify2array(CSRM_gmse_sampling), c(1,2), mean)
  
  
  pem_mse_model[[B]] = apply(simplify2array(pem_mse_sampling), c(1,2), mean)
  CSRM_pem_mse_model[[B]] = apply(simplify2array(CSRM_pem_mse_sampling), c(1,2), mean)
  pem_bias_model[[B]] = apply(simplify2array(pem_bias_sampling), c(1,2), mean)
  pem_variance_model[[B]] = apply(simplify2array(pem_variance_sampling), c(1,2), mean)
  CSRM_pem_bias_model[[B]] = apply(simplify2array(CSRM_pem_bias_sampling), c(1,2), mean)
  CSRM_pem_variance_model[[B]] = apply(simplify2array(CSRM_pem_variance_sampling), c(1,2), mean)
  
  transposed = lapply(seq_along(P_sampling[[1]]), function(j) {
    lapply(P_sampling, `[[`, j)
  })
  P_model[[B]] = lapply(transposed, function(mats) {
    Reduce("+", mats) / length(mats)
  })
  
}


list_of_results_model = mget(results_model_names)

results_model = list()

for (i in seq_along(list_of_results_model)) {
  result = list_of_results_model[[i]]
  result_name = results_model_names[i]
  
  result_df = bind_rows(lapply(seq_along(result), function(j) {
    df = as.data.frame(result[[j]])
    df$Iteration = j
    df$Estimate = result_name
    df$Domain = rownames(df)
    rownames(df) = NULL
    return(df)
  }))
  
  results_model[[i]] = result_df
}

results_model_df = bind_rows(results_model)
results_joint_df = results_model_df %>% group_by(Estimate, Domain) %>%
  summarise(across(all_of(levels(pop$onderwijs)), mean))

list_results_sampling = mget(results_sampling_names)
results_sampling = list()

for (i in seq_along(list_results_sampling)) {
  result = list_results_sampling[[i]]
  result_name = results_sampling_names[i]
  
  result_df = bind_rows(lapply(seq_along(result), function(j) {
    df = as.data.frame(result[[j]])
    df$Iteration = j
    df$Estimate = result_name
    df$Domain = rownames(df)
    rownames(df) = NULL
    return(df)
  }))
  
  results_sampling[[i]] = result_df
}

results_sampling_df = bind_rows(results_sampling)  

save(P_model, P_sampling, file = "sup1_c3_P.RData")
write.csv(results_sampling_df, "sup1_c3_sampling_distribution.csv")
write.csv(results_model_df, "sup1_c3_model_distribution.csv")
write.csv(results_joint_df, "sup1_c3_joint_distribution.csv")



################################################################################################

# Bind the datasets based on the distribution

model_paths = c(
  "sup1_baseline_model_distribution.csv",
  "sup1_small_sample_model_distribution.csv",
  "sup1_small_population_model_distribution.csv",
  "sup1_nonignorable_model_distribution.csv",
  "sup1_overparameterisation_model_distribution.csv",
  "sup1_designbias_model_distribution.csv",
  "sup1_modelbias_small_model_distribution.csv",
  "sup1_modelbias_large_model_distribution.csv",
  "sup1_modelbias_large_small_samp_model_distribution.csv"
)

sampling_paths = gsub("_model_", "_sampling_", model_paths)
joint_paths = gsub("_model_", "_joint_", model_paths)

conditions = c(
  "Baseline",
  "Small sample",
  "Small population",
  "Non-ignorable",
  "Overparameterization",
  "Design bias",
  "Model bias (small)",
  "Model bias (large)",
  "Model bias (large) X small sample")


model_distribution_sup1 = do.call(rbind, lapply(seq_along(model_paths), function(i) {
  df = read.csv(model_paths[i])
  df$Condition <- conditions[i]
  return(df)
}))


sampling_distribution_sup1 = do.call(rbind, lapply(seq_along(sampling_paths), function(i) {
  df = read.csv(sampling_paths[i])
  df$Condition = conditions[i]
  return(df)
}))


joint_distribution_sup1 = do.call(rbind, lapply(seq_along(joint_paths), function(i) {
  df = read.csv(joint_paths[i])
  df$Condition = conditions[i]
  return(df)
}))


# Rename the variables

joint_distribution_sup1$Estimate = as.factor(joint_distribution_sup1$Estimate)
levels(joint_distribution_sup1$Estimate) = c(
  "Bias_joint", "RB_joint", "RRMSE_GMSE",
  "RRMSE_joint", "RB_PEM", "RRMSE_PEM",
  "CV_PEM", "CV_joint", "CV_Vdesign",
  "GMSE", "MSE_joint", "PEM_bias",
  "PEM_MSE", "PEM_variance", "Est_total",
  "True_total", "Variance_joint", "Vdesign"
)

sampling_distribution_sup1$Estimate = as.factor(sampling_distribution_sup1$Estimate)
levels(sampling_distribution_sup1$Estimate) = c(
  "RRMSE_GMSE", "RB_PEM", "RRMSE_PEM",
  "CV_PEM", "CV_Vdesign", "GMSE",
  "PEM_bias", "PEM_MSE", "PEM_variance",
  "Est_total", "Vdesign"
)

model_distribution_sup1$Estimate = as.factor(model_distribution_sup1$Estimate)
levels(model_distribution_sup1$Estimate) = c(
  "Bias_joint", "RB_joint", "RRMSE_GMSE",
  "RRMSE_joint", "RB_PEM", "RRMSE_PEM",
  "CV_PEM", "CV_joint", "CV_Vdesign",
  "GMSE", "MSE_joint", "PEM_bias",
  "PEM_MSE", "PEM_variance", "Est_total",
  "True_total", "Variance_joint", "Vdesign"
)

results_sampling_df$Estimate = as.factor(results_sampling_df$Estimate)
levels(results_sampling_df$Estimate) = c(
  "RRMSE_GMSE", "RB_PEM", "RRMSE_PEM",
  "CV_PEM", "CV_Vdesign", "GMSE",
  "PEM_bias", "PEM_MSE", "PEM_variance",
  "Est_total", "Vdesign"
)
sampling_distribution_sup1_c3 = results_sampling_df

results_model_df$Estimate = as.factor(results_model_df$Estimate)
levels(results_model_df$Estimate) = c(
  "Bias_joint", "RB_joint", "RRMSE_GMSE",
  "RRMSE_joint", "RB_PEM", "RRMSE_PEM",
  "CV_PEM", "CV_joint", "CV_Vdesign",
  "GMSE", "MSE_joint", "PEM_bias",
  "PEM_MSE", "PEM_variance", "Est_total",
  "True_total", "Variance_joint", "Vdesign"
)

sampling_distribution_sup1_c3_summary = results_model_df

# save the files for analysis
save(sampling_distribution_sup1,
     model_distribution_sup1,
     joint_distribution_sup1,
     sampling_distribution_sup1_c3,
     sampling_distribution_sup1_c3_summary,
     file = "sup1_data.RData")
