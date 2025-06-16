#install.packages("tidyverse)
library(tidyverse)


##############################################################################################

                                    # Superpopulation 1

load("sup1_data.RData")

# Estimated and true totals; binomial target variable

detach("package:MASS", unload = TRUE) #conflict with dplyr::select


totals_joint = joint_distribution_sup1 %>% filter(Estimate %in% c("Est_total", "True_total")) %>%
  rename("Value" = "joint") %>%
  pivot_wider(
    names_from = Estimate,
    values_from = Value,
    id_cols = c(Domain, Condition)
  ) %>% select(Condition, Domain, True_total, Est_total)

totals_sampling = model_distribution_sup1 %>% 
  filter(Estimate %in% c("Est_total", "True_total") & Iteration == 100) %>%
  rename("Value" = "Wo") %>%
  pivot_wider(
    names_from = Estimate,
    values_from = Value,
    id_cols = c(Domain, Condition)
  ) %>% select(Condition, Domain, True_total, Est_total)

# Component specific relative measures (benchmarks and their approximations); binomial target variable

estimators_joint = joint_distribution_sup1 %>% 
  filter(Estimate %in% c("RRMSE_joint", "CV_joint", "RB_joint",
                         "RRMSE_GMSE", "CV_Vdesign", "RRMSE_PEM",
                         "CV_PEM", "RB_PEM")) %>%
  rename("Value" = "joint") %>%
  mutate(Value = Value*100) %>% #convert to percentage
  pivot_wider(
    names_from = Estimate,
    values_from = Value,
    id_cols = c(Domain, Condition)
  ) %>% select(Condition, Domain, RRMSE_joint, RB_joint, CV_joint, CV_Vdesign, RRMSE_GMSE, RRMSE_PEM, RB_PEM, CV_PEM)


estimators_sampling = model_distribution_sup1 %>% 
  filter(Estimate %in% c("RRMSE_joint", "CV_joint", "RB_joint",
                         "RRMSE_GMSE", "CV_Vdesign", "RRMSE_PEM",
                         "CV_PEM", "RB_PEM") & Iteration == 100) %>%
  rename("Value" = "Wo") %>%
  mutate(Value = Value*100) %>% #convert to percentage
  pivot_wider(
    names_from = Estimate,
    values_from = Value,
    id_cols = c(Domain, Condition)
  ) %>% select(Condition, Domain, RRMSE_joint, RB_joint, CV_joint, CV_Vdesign, RRMSE_GMSE, RRMSE_PEM, RB_PEM, CV_PEM)


# Multinomial target variable

sup1_c3 = read.csv("sup1_c3_joint_distribution.csv")

sup1_c3$Estimate = as.factor(sup1_c3$Estimate)
levels(sup1_c3$Estimate) = c(
  "Bias_joint", "RB_joint", "RRMSE_GMSE",
  "RRMSE_joint", "RB_PEM", "RRMSE_PEM",
  "CV_PEM", "CV_joint", "CV_Vdesign",
  "GMSE", "MSE_joint", "PEM_bias",
  "PEM_MSE", "PEM_variance", "Est_total",
  "True_total", "Variance_joint", "Vdesign"
)

totals_sampling = sup1_c3 %>% filter(Estimate %in% c("Est_total", "True_total"))

estimates_sampling = sup1_c3 %>%
  filter(Estimate %in% c("RRMSE_joint", "CV_joint", "RB_joint",
                         "RRMSE_GMSE", "CV_Vdesign", "RRMSE_PEM",
                         "CV_PEM", "RB_PEM")) %>%
  mutate(across(c(Low, Middle, High), ~ .x * 100))
  

##############################################################################################

# Superpopulation 2

load("sup2_data.RData")

# Estimated and true totals; binomial target variable

detach("package:MASS", unload = TRUE) #conflict with dplyr::select


totals_joint = joint_distribution_sup2 %>% filter(Estimate %in% c("Est_total", "True_total")) %>%
  rename("Value" = "joint") %>%
  pivot_wider(
    names_from = Estimate,
    values_from = Value,
    id_cols = c(Domain, Condition)
  ) %>% select(Condition, Domain, True_total, Est_total)

totals_sampling = model_distribution_sup2 %>% 
  filter(Estimate %in% c("Est_total", "True_total") & Iteration == 100) %>%
  rename("Value" = "Wo") %>%
  pivot_wider(
    names_from = Estimate,
    values_from = Value,
    id_cols = c(Domain, Condition)
  ) %>% select(Condition, Domain, True_total, Est_total)

# Component specific relative measures (benchmarks and their approximations); binomial target variable

estimators_joint = joint_distribution_sup2 %>% 
  filter(Estimate %in% c("RRMSE_joint", "CV_joint", "RB_joint",
                         "RRMSE_GMSE", "CV_Vdesign", "RRMSE_PEM",
                         "CV_PEM", "RB_PEM")) %>%
  rename("Value" = "joint") %>%
  mutate(Value = Value*100) %>% #convert to percentage
  pivot_wider(
    names_from = Estimate,
    values_from = Value,
    id_cols = c(Domain, Condition)
  ) %>% select(Condition, Domain, RRMSE_joint, RB_joint, CV_joint, CV_Vdesign, RRMSE_GMSE, RRMSE_PEM, RB_PEM, CV_PEM)


estimators_sampling = model_distribution_sup2 %>% 
  filter(Estimate %in% c("RRMSE_joint", "CV_joint", "RB_joint",
                         "RRMSE_GMSE", "CV_Vdesign", "RRMSE_PEM",
                         "CV_PEM", "RB_PEM") & Iteration == 100) %>%
  rename("Value" = "Wo") %>%
  mutate(Value = Value*100) %>% #convert to percentage
  pivot_wider(
    names_from = Estimate,
    values_from = Value,
    id_cols = c(Domain, Condition)
  ) %>% select(Condition, Domain, RRMSE_joint, RB_joint, CV_joint, CV_Vdesign, RRMSE_GMSE, RRMSE_PEM, RB_PEM, CV_PEM)


# Multinomial target variable

sup2_c3 = read.csv("sup2_c3_joint_distribution.csv")

sup2_c3$Estimate = as.factor(sup2_c3$Estimate)
levels(sup2_c3$Estimate) = c(
  "Bias_joint", "RB_joint", "RRMSE_GMSE",
  "RRMSE_joint", "RB_PEM", "RRMSE_PEM",
  "CV_PEM", "CV_joint", "CV_Vdesign",
  "GMSE", "MSE_joint", "PEM_bias",
  "PEM_MSE", "PEM_variance", "Est_total",
  "True_total", "Variance_joint", "Vdesign"
)

totals_sampling = sup2_c3 %>% filter(Estimate %in% c("Est_total", "True_total"))

estimates_sampling = sup2_c3 %>%
  filter(Estimate %in% c("RRMSE_joint", "CV_joint", "RB_joint",
                         "RRMSE_GMSE", "CV_Vdesign", "RRMSE_PEM",
                         "CV_PEM", "RB_PEM")) %>%
  mutate(across(c(Low, Middle, High), ~ .x * 100))

  
  