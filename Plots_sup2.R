#install.packages(tidyverse)
#install.packages(patchwork)
library(tidyverse)
library(patchwork)
detach("package:MASS", unload = TRUE) #conflict with dplyr::select

load("sup2_data.RData")

################################################################################################


# The script below produces the plots for change in benchmark estimators, baseline condition,
# multinomial condition, variance and bias condition for superpopulation 2.


################################################################################################

# prepare plot data

custom_colors_male <- c(
  "RRMSE_GMSE" = "#1874CD",
  "CV_Vdesign" = "#7AC5CD",
  "RRMSE_PEM" = "deepskyblue3"
)

custom_colors_female <- c(
  "RRMSE_GMSE" = "deeppink4",
  "RRMSE_PEM" = "hotpink2",
  "CV_Vdesign" = "#CD6090"
)

sampling_distribution_sup2$Domain = ifelse(sampling_distribution_sup2$Domain == "A", "Man",  # label A and B in domain_random Male and Female for convenience
                                           ifelse(sampling_distribution_sup2$Domain == "B", "Vrouw",
                                                  sampling_distribution_sup2$Domain))
sampling_distribution_sup2$Domain = as.factor(sampling_distribution_sup1$Domain)
levels(sampling_distribution_sup2$Domain) = c("Male", "Female")

model_distribution_sup2$Domain = ifelse(model_distribution_sup2$Domain == "A", "Man", 
                                        ifelse(model_distribution_sup2$Domain == "B", "Vrouw",
                                               model_distribution_sup2$Domain))
model_distribution_sup2$Domain = as.factor(model_distribution_sup2$Domain)
levels(model_distribution_sup2$Domain) = c("Male", "Female")

joint_distribution_sup2$Domain = ifelse(joint_distribution_sup2$Domain == "A", "Man", 
                                        ifelse(joint_distribution_sup2$Domain == "B", "Vrouw",
                                               joint_distribution_sup2$Domain))
joint_distribution_sup2$Domain = as.factor(joint_distribution_sup2$Domain)
levels(joint_distribution_sup2$Domain) = c("Male", "Female")

#############################################################################################

# Benchmark estimators

CV_lineplot = joint_distribution_sup2 %>% 
  filter(Estimate %in% c("RRMSE_joint", "CV_joint", "RB_joint")) %>%
  mutate(Estimate_perc = abs(joint*100))

CV_lineplot_samp = model_distribution_sup2 %>% 
  filter(Iteration == 100 & Estimate %in% c("RRMSE_joint", "CV_joint", "RB_joint")) %>%
  mutate(Estimate_perc = abs(joint*100))

CV_lineplot$Condition = factor(CV_lineplot$Condition, levels=c("Baseline",
                                                               "Small sample",
                                                               "Small population",
                                                               "Non-ignorable",
                                                               "Overparameterization",
                                                               "Design bias",
                                                               "Model bias (small)",
                                                               "Model bias (large)",
                                                               "Model bias (large) X small sample"))

CV_lineplot_samp$Condition = factor(CV_lineplot_samp$Condition, levels = c("Baseline",
                                                                           "Small sample",
                                                                           "Small population",
                                                                           "Non-ignorable",
                                                                           "Overparameterization",
                                                                           "Design bias",
                                                                           "Model bias (small)",
                                                                           "Model bias (large)",
                                                                           "Model bias (large) X small sample"))


# Joint distribution
ggplot() +
  
  geom_line(data = CV_lineplot, aes(
    x = Condition, y = Estimate_perc,
    group = interaction(Estimate, Domain),
    color = Domain,
    linetype = Estimate
  ),
  linewidth = 1, alpha = 0.9) +
  
  geom_point(data = CV_lineplot, aes(
    x = Condition, y = Estimate_perc,
    group = interaction(Estimate, Domain),
    color = Domain
  ),
  size = 2) +
  
  scale_color_manual(values = c("Male" = "#27408B", "Female" = "#8B3A62")) +
  scale_linetype_manual(values = c(
    "RRMSE_joint" = "solid",
    "RB_joint" = "dashed",
    "CV_joint" = "dotted"
  )) +
  labs(
    x = NULL,
    y = "CSRM (%)",
  ) +
  facet_grid(~ Domain) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 15, hjust = 0),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14)
  ) + guides(color = "none")



# Sampling distribution
ggplot() +
  geom_line(data = CV_lineplot_samp, aes(
    x = Condition, y = Estimate_perc,
    group = interaction(Estimate, Domain),
    color = Domain,
    linetype = Estimate
  ),
  size = 1, alpha = 0.85) +
  
  geom_point(data = CV_lineplot_samp, aes(
    x = Condition, y = Estimate_perc,
    group = interaction(Estimate, Domain),
    color = Domain
  ),
  size = 2) +
  scale_color_manual(values = c("Male" = "#27408B", "Female" = "#8B3A62")) +
  scale_linetype_manual(values = c(
    "MSE" = "solid",
    "Bias" = "dashed",
    "Variance" = "dotted"
  )) +
  labs(
    x = NULL,
    y = "CSRM (%)",
  ) +
  facet_grid(~ Domain) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(size = 15, hjust = 0),
    legend.text = element_text(size = 10)
  ) + guides(color = "none")



#############################################################################################

### Baseline condition

# sampling distribution


boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Baseline") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Baseline")) %>% mutate(Estimate_perc = Wo*100)



ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +
  scale_y_continuous(limits = c(1, 2)) +  
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 16),     
    axis.text = element_text(size = 14),       
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Baseline") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Baseline")) %>% mutate(Estimate_perc = Wo*100)


ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +
  scale_y_continuous(limits = c(0.95, 1.9)) +  
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 16),   
    axis.text = element_text(size = 14),         
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


# joint distribution

boxplot = model_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Baseline") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Baseline")) %>% mutate(Estimate_perc = joint*100)



ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +
  scale_y_continuous(limits = c(1, 2)) +  
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 16), 
    axis.text = element_text(size = 14),         
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = model_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Baseline") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Baseline")) %>% mutate(Estimate_perc = joint*100)


ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +
  scale_y_continuous(limits = c(0.95, 1.9)) +  
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 16),       
    axis.text = element_text(size = 14),      
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


# PEM (bias)

# Sampling distribution

boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Baseline") &
           Estimate == "RB_PEM") %>% mutate(Estimate_perc = abs(Wo*100)) # plot absolute bias for convenience

boxplot$Estimate = droplevels(boxplot$Estimate)


benchmark_bias = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RB_joint" & Domain == "Male" & 
           Condition %in%  c("Baseline")) %>% mutate(Estimate_perc = abs(Wo*100))

ggplot(boxplot, aes(x = Condition, y = Estimate_perc)) +
  geom_boxplot(
    fill = "deepskyblue3",
    position = position_dodge(width = 0.8),
    width = 0.6
  ) +
  geom_hline(
    data = benchmark_bias,
    aes(yintercept = Estimate_perc),
    linetype = "solid",
    color = "orange",
    linewidth = 1.25
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(y = "RB (%)", x = "PEM (bias)") +
  theme(
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 16),
    legend.position = "none"
  )


boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Baseline") &
           Estimate == "RB_PEM") %>% mutate(Estimate_perc = abs(Wo*100))

boxplot$Estimate = droplevels(boxplot$Estimate)


benchmark_bias = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RB_joint" & Domain == "Female" & 
           Condition %in%  c("Baseline")) %>% mutate(Estimate_perc = abs(Wo*100))

ggplot(boxplot, aes(x = Condition, y = Estimate_perc)) +
  geom_boxplot(
    fill = "hotpink2",
    position = position_dodge(width = 0.8),
    width = 0.6
  ) +
  geom_hline(
    data = benchmark_bias,
    aes(yintercept = Estimate_perc),
    linetype = "solid",
    color = "orange",
    linewidth = 1.25
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(y = "RB (%)", x = "PEM (bias)") +
  theme(
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )

# Joint distribution

boxplot = model_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Baseline") &
           Estimate == "RB_PEM") %>% mutate(Estimate_perc = abs(Wo*100))

boxplot$Estimate = droplevels(boxplot$Estimate)


benchmark_bias = joint_distribution_sup2 %>%
  filter(Estimate == "RB_joint" & Domain == "Male" & 
           Condition %in%  c("Baseline")) %>% mutate(Estimate_perc = abs(joint*100))


ggplot(boxplot, aes(x = Condition, y = Estimate_perc)) +
  geom_boxplot(
    fill = "deepskyblue3",
    position = position_dodge(width = 0.8),
    width = 0.6
  ) +
  geom_hline(
    data = benchmark_bias,
    aes(yintercept = Estimate_perc),
    linetype = "solid",
    color = "orange",
    linewidth = 1.25
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(y = "RB (%)", x = "PEM (bias)") +
  theme(
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 16),
    legend.position = "none"
  )

boxplot = model_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Baseline") &
           Estimate == "RB_PEM") %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)


benchmark_bias = joint_distribution_sup2 %>%
  filter(Estimate == "RB_joint" & Domain == "Female" & 
           Condition %in%  c("Baseline")) %>% mutate(Estimate_perc = joint*100)

ggplot(boxplot, aes(x = Condition, y = Estimate_perc)) +
  geom_boxplot(
    fill = "hotpink2",
    position = position_dodge(width = 0.8),
    width = 0.6
  ) +
  geom_hline(
    data = benchmark_bias,
    aes(yintercept = Estimate_perc),
    linetype = "solid",
    color = "orange",
    linewidth = 1.25
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(y = "RB (%)", x = "PEM (bias)") +
  theme(
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )

#############################################################################################

### Multinomial condition

# prepare data for plotting

sampling_distribution_sup2_c3$Domain = as.factor(sampling_distribution_sup2_c3$Domain)
levels(sampling_distribution_sup2_c3$Domain) = c("Male", "Female")

c3_sup2 = sampling_distribution_sup2_c3 %>%
  select(Low, Middle, High, Estimate, Domain) %>%
  filter(Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>%
  mutate(Low = Low*100,
         Middle = Middle*100,
         High = High*100) 


c3_sup2$Estimate = as.factor(c3_sup2$Estimate)
c3_sup2$Estimate = droplevels(c3_sup2$Estimate)
c3_sup2$Estimate = factor(c3_sup2$Estimate, levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


c3_sup2_long <- c3_sup2 %>%
  pivot_longer(
    cols = c(Low, Middle, High),
    names_to = "Category",
    values_to = "Estimate_perc"
  ) %>%
  mutate(
    fill_group = interaction(Domain, Estimate, sep = " - ")
  )

c3_sup2_long$Category = as.factor(c3_sup2_long$Category)
c3_sup2_long$Category = factor(c3_sup2_long$Category,
                               levels = c("Low",
                                          "Middle",
                                          "High"))


fill_colors = c(
  "Male - RRMSE_GMSE" = "#1874CD",
  "Male - CV_Vdesign" =  "#7AC5CD",
  "Male - RRMSE_PEM" = "deepskyblue3",
  "Female - RRMSE_GMSE" = "deeppink4",
  "Female - CV_Vdesign" = "hotpink2",
  "Female - RRMSE_PEM" = "#CD6090"
)



label_data <- data.frame(
  Category = rep(c("Low", "Middle", "High"), each = 3),
  y = rep(c(3.5, 1.5, 2.2), each = 3),
  label = rep(c("Vdesign  GMSE  PEM (MSE)"), times = 3),
  Domain = rep(c("Male", "Female"), each = 9)  
)

label_data[label_data$Category == "Low" & label_data$Domain == "Female", "y"] = 4
label_data$Domain = as.factor(label_data$Domain)
label_data$Domain = factor(label_data$Domain,
                           levels = c("Male", "Female"))

ggplot(c3_sup2_long, aes(x = Category, y = Estimate_perc, fill = fill_group)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  facet_wrap(~ Domain) +
  scale_fill_manual(values = custom_fill_colors) +
  
  geom_text(
    data = label_data,
    aes(x = Category, y = y, label = label),
    size = 4,
    hjust = 0.45,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = data.frame(
      x = c(0.4, 0.4, 1.4, 1.4, 2.4, 2.4),
      xend = c(1.6, 1.6, 2.6, 2.6, 3.6, 3.6),
      y = c(3.87, 4.77, 1.87, 1.83, 2.52, 2.51),
      yend = c(3.87, 4.77, 1.87, 1.83, 2.52, 2.51),
      Domain = factor(c("Male", "Female", "Male", "Female", "Male", "Female"), levels = c("Male", "Female"))      
    ),
    aes(x = x, xend = xend, y = y, yend = yend),
    color = "orange",
    linewidth = 1,
    inherit.aes = FALSE
  )+
  
  labs(
    x = NULL,
    y = "CSRM (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 14)
  )

# Bias

c3_sup2 = sampling_distribution_sup2_c3 %>%
  select(Low, Middle, High, Estimate, Domain) %>%
  filter(Estimate %in% c("RB_PEM")) %>%
  mutate(Low = abs(Low*100),
         Middle = abs(Middle*100),
         High = abs(High*100)) 


c3_sup2$Estimate = as.factor(c3_sup2$Estimate)
c3_sup2$Estimate = droplevels(c3_sup2$Estimate)


c3_sup2_long <- c3_sup2 %>%
  pivot_longer(
    cols = c(Low, Middle, High),
    names_to = "Category",
    values_to = "Estimate_perc"
  ) %>%
  mutate(
    fill_group = interaction(Domain, Estimate, sep = " - ")
  )

c3_sup2_long$Category = as.factor(c3_sup2_long$Category)
c3_sup2_long$Category = factor(c3_sup2_long$Category,
                               levels = c("Low",
                                          "Middle",
                                          "High"))

ggplot(c3_sup2_long, aes(x = Category, y = Estimate_perc, fill = Domain)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  facet_wrap(~ Domain)  +
  geom_segment(
    data = data.frame(
      x = c(0.4, 0.4, 1.4, 1.4, 2.4, 2.4),
      xend = c(1.6, 1.6, 2.6, 2.6, 3.6, 3.6),
      y = c(0.10387249, 0.06448135, 0.02296947, 0.06497989, 0.09944083, 0.05208205),
      yend = c(0.10387249, 0.06448135, 0.02296947, 0.06497989, 0.09944083, 0.05208205),
      Domain = factor(c("Male", "Female", "Male", "Female", "Male", "Female"), levels = c("Male", "Female"))      # facet value to match
    ),
    aes(x = x, xend = xend, y = y, yend = yend),
    color = "orange",
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c(Male = "deepskyblue3","Female" = "#CD6090")) +
  
  labs(
    x = NULL,
    y = "RB (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 14)
  )


#############################################################################################

### Variance conditions

#sampling distribution; male

boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Small sample") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Small sample")) %>% mutate(Estimate_perc = Wo*100)




small_sample_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(2, 4.5)) +
  theme(
    axis.title.y = element_text(size = 16),   
    axis.text = element_text(size = 14),     
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )

boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Small population") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Small population")) %>% mutate(Estimate_perc = Wo*100)


small_population_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(1.5, 3)) +
  theme(
    axis.title.y = element_text(size = 16),    
    axis.text = element_text(size = 14),     
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )

boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Non-ignorable") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Non-ignorable")) %>% mutate(Estimate_perc = Wo*100)




non_ignorable_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") + 
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(1, 3.1)) +
  theme(
    axis.title.y = element_text(size = 16),    
    axis.text = element_text(size = 14),    
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Overparameterization") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Overparameterization")) %>% mutate(Estimate_perc = Wo*100)




overparameterization_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(1.2, 2)) +
  theme(
    axis.title.y = element_text(size = 16),    
    axis.text = element_text(size = 14),     
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


combined_plot = (small_sample_plot | small_population_plot) /
  (non_ignorable_plot | overparameterization_plot)


combined_plot


# sampling distribution; female

boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Small sample") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Small sample")) %>% mutate(Estimate_perc = Wo*100)


small_sample_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(2, 4.5)) +
  theme(
    axis.title.y = element_text(size = 16),   
    axis.text = element_text(size = 14),     
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Small population") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Small population")) %>% mutate(Estimate_perc = Wo*100)



small_population_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") + 
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(1.5, 3)) +
  theme(
    axis.title.y = element_text(size = 16),  
    axis.text = element_text(size = 14),          
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )



boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Non-ignorable") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Non-ignorable")) %>% mutate(Estimate_perc = Wo*100)



non_ignorable_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") + 
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(1, 3.1)) +
  theme(
    axis.title.y = element_text(size = 16),   
    axis.text = element_text(size = 14),   
    strip.text = element_text(size = 16), 
    legend.position = "none"
  )



boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Overparameterization") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Overparameterization")) %>% mutate(Estimate_perc = Wo*100)

overparameterization_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(1.2, 2)) +
  theme(
    axis.title.y = element_text(size = 16),   
    axis.text = element_text(size = 14),   
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


combined_plot = (small_sample_plot | small_population_plot) /
  (non_ignorable_plot | overparameterization_plot)


combined_plot


### joint distribution; male

boxplot = model_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Small sample") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Small sample")) %>% mutate(Estimate_perc = joint*100)



small_sample_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") + 
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(2, 4.5)) +
  theme(
    axis.title.y = element_text(size = 16),   
    axis.text = element_text(size = 14),   
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )

boxplot = model_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Small population") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Small population")) %>% mutate(Estimate_perc = joint*100)


small_population_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(1.5, 3)) +
  theme(
    axis.title.y = element_text(size = 16),   
    axis.text = element_text(size = 14),  
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = model_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Non-ignorable") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Non-ignorable")) %>% mutate(Estimate_perc = joint*100)



non_ignorable_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(1, 3.1)) +
  theme(
    axis.title.y = element_text(size = 16),  
    axis.text = element_text(size = 14),    
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = model_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Overparameterization") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Overparameterization")) %>% mutate(Estimate_perc = joint*100)



overparameterization_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(1.2, 2)) +
  theme(
    axis.title.y = element_text(size = 16),  
    axis.text = element_text(size = 14),     
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )



combined_plot = (small_sample_plot | small_population_plot) /
  (non_ignorable_plot | overparameterization_plot)

combined_plot


### joint distribution; female

boxplot = model_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Small sample") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Small sample")) %>% mutate(Estimate_perc = joint*100)




small_sample_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") + 
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(2, 4.5)) +
  theme(
    axis.title.y = element_text(size = 16),    
    axis.text = element_text(size = 14),      
    strip.text = element_text(size = 16), 
    legend.position = "none"
  )

boxplot = model_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Small population") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Small population")) %>% mutate(Estimate_perc = joint*100)



small_population_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(1.5, 3)) +
  theme(
    axis.title.y = element_text(size = 16),  
    axis.text = element_text(size = 14),         
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = model_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Non-ignorable") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Non-ignorable")) %>% mutate(Estimate_perc = joint*100)


non_ignorable_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(1, 3.1)) +
  theme(
    axis.title.y = element_text(size = 16),   
    axis.text = element_text(size = 14),   
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = model_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Overparameterization") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Overparameterization")) %>% mutate(Estimate_perc = joint*100)



overparameterization_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(1.2, 2)) +
  theme(
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 14),      
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


combined_plot = (small_sample_plot | small_population_plot) /
  (non_ignorable_plot | overparameterization_plot)


combined_plot


####### Bias conditions

# sampling distribution; male

boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Design bias") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Design bias")) %>% mutate(Estimate_perc = Wo*100)


benchmark_var = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "CV_joint" & Domain == "Male" & 
           Condition %in%  c("Design bias")) %>% mutate(Estimate_perc = Wo*100)

design_bias_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 2)) +
  theme(
    axis.title.y = element_text(size = 16),   
    axis.text = element_text(size = 14),   
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )

boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Model bias (small)") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Model bias (small)")) %>% mutate(Estimate_perc = Wo*100)


benchmark_var = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "CV_joint" & Domain == "Male" & 
           Condition %in%  c("Model bias (small)")) %>% mutate(Estimate_perc = Wo*100)




mbs_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") + 
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 2)) +
  theme(
    axis.title.y = element_text(size = 16), 
    axis.text = element_text(size = 14),     
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Model bias (large)") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Model bias (large)")) %>% mutate(Estimate_perc = Wo*100)


benchmark_var = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "CV_joint" & Domain == "Male" & 
           Condition %in%  c("Model bias (large)")) %>% mutate(Estimate_perc = Wo*100)


mbl_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 7)) +
  theme(
    axis.title.y = element_text(size = 16), 
    axis.text = element_text(size = 14),     
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Model bias (large) X small sample") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Model bias (large) X small sample")) %>% mutate(Estimate_perc = Wo*100)


benchmark_var = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "CV_joint" & Domain == "Male" & 
           Condition %in%  c("Model bias (large) X small sample")) %>% mutate(Estimate_perc = Wo*100)



mblss_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 7)) +
  theme(
    axis.title.y = element_text(size = 16),  
    axis.text = element_text(size = 14),  
    strip.text = element_text(size = 14),  
    legend.position = "none"
  )



combined_plot = (design_bias_plot | mbs_plot) /
  (mbl_plot | mblss_plot)


combined_plot


### sampling distribution; female

boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Design bias") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Design bias")) %>% mutate(Estimate_perc = Wo*100)


benchmark_var = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "CV_joint" & Domain == "Female" & 
           Condition %in%  c("Design bias")) %>% mutate(Estimate_perc = Wo*100)


design_bias_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") + 
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 2)) +
  theme(
    axis.title.y = element_text(size = 16),  
    axis.text = element_text(size = 14),   
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )

boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Model bias (small)") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Model bias (small)")) %>% mutate(Estimate_perc = Wo*100)


benchmark_var = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "CV_joint" & Domain == "Female" & 
           Condition %in%  c("Model bias (small)")) %>% mutate(Estimate_perc = Wo*100)

mbs_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") + 
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 2)) +
  theme(
    axis.title.y = element_text(size = 16),  
    axis.text = element_text(size = 14),  
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Model bias (large)") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Model bias (large)")) %>% mutate(Estimate_perc = Wo*100)


benchmark_var = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "CV_joint" & Domain == "Female" & 
           Condition %in%  c("Model bias (large)")) %>% mutate(Estimate_perc = Wo*100)


mbl_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 10)) +
  theme(
    axis.title.y = element_text(size = 16),  
    axis.text = element_text(size = 14),          
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = sampling_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Model bias (large) X small sample") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Model bias (large) X small sample")) %>% mutate(Estimate_perc = Wo*100)


benchmark_var = model_distribution_sup2 %>%
  filter(Iteration == 100 & Estimate == "CV_joint" & Domain == "Female" & 
           Condition %in%  c("Model bias (large) X small sample")) %>% mutate(Estimate_perc = Wo*100)

mblss_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 10)) +
  theme(
    axis.title.y = element_text(size = 16),    
    axis.text = element_text(size = 14),   
    strip.text = element_text(size = 14),  
    legend.position = "none"
  )



combined_plot = (design_bias_plot | mbs_plot) /
  (mbl_plot | mblss_plot)


combined_plot


### joint distribution

boxplot = model_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Design bias") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Design bias")) %>% mutate(Estimate_perc = joint*100)


benchmark_var = joint_distribution_sup2 %>%
  filter(Estimate == "CV_joint" & Domain == "Male" & 
           Condition %in%  c("Design bias")) %>% mutate(Estimate_perc = joint*100)


design_bias_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") + 
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 2)) +
  theme(
    axis.title.y = element_text(size = 16),   
    axis.text = element_text(size = 14), 
    strip.text = element_text(size = 16), 
    legend.position = "none"
  )

boxplot = model_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Model bias (small)") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Model bias (small)")) %>% mutate(Estimate_perc = joint*100)


benchmark_var = joint_distribution_sup2 %>%
  filter(Estimate == "CV_joint" & Domain == "Male" & 
           Condition %in%  c("Model bias (small)")) %>% mutate(Estimate_perc = joint*100)




mbs_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 2)) +
  theme(
    axis.title.y = element_text(size = 16),   
    axis.text = element_text(size = 14),     
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = model_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Model bias (large)") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Model bias (large)")) %>% mutate(Estimate_perc = joint*100)


benchmark_var = joint_distribution_sup2 %>%
  filter(Estimate == "CV_joint" & Domain == "Male" & 
           Condition %in%  c("Model bias (large)")) %>% mutate(Estimate_perc = joint*100)


mbl_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 7)) +
  theme(
    axis.title.y = element_text(size = 16),  
    axis.text = element_text(size = 14),      
    strip.text = element_text(size = 16), 
    legend.position = "none"
  )


boxplot = model_distribution_sup2 %>%
  filter(Domain == "Male" &                                           # Male
           Condition %in% c("Model bias (large) X small sample") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Male" & 
           Condition %in%  c("Model bias (large) X small sample")) %>% mutate(Estimate_perc = joint*100)


benchmark_var = joint_distribution_sup2 %>%
  filter(Estimate == "CV_joint" & Domain == "Male" & 
           Condition %in%  c("Model bias (large) X small sample")) %>% mutate(Estimate_perc = joint*100)


mblss_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_male) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 7)) +
  theme(
    axis.title.y = element_text(size = 16),  
    axis.text = element_text(size = 14),   
    strip.text = element_text(size = 14),  
    legend.position = "none"
  )



combined_plot = (design_bias_plot | mbs_plot) /
  (mbl_plot | mblss_plot)

combined_plot


### model distribution; female

boxplot = model_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Design bias") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Design bias")) %>% mutate(Estimate_perc = joint*100)


benchmark_var = joint_distribution_sup2 %>%
  filter(Estimate == "CV_joint" & Domain == "Female" & 
           Condition %in%  c("Design bias")) %>% mutate(Estimate_perc = joint*100)


design_bias_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") + 
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 2)) +
  theme(
    axis.title.y = element_text(size = 16), 
    axis.text = element_text(size = 14),          
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )

boxplot = model_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Model bias (small)") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Model bias (small)")) %>% mutate(Estimate_perc = joint*100)


benchmark_var = joint_distribution_sup2 %>%
  filter(Estimate == "CV_joint" & Domain == "Female" & 
           Condition %in%  c("Model bias (small)")) %>% mutate(Estimate_perc = joint*100)




mbs_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 2)) +
  theme(
    axis.title.y = element_text(size = 16),  
    axis.text = element_text(size = 14),       
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = model_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Model bias (large)") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Model bias (large)")) %>% mutate(Estimate_perc = joint*100)


benchmark_var = joint_distribution_sup2 %>%
  filter(Estimate == "CV_joint" & Domain == "Female" & 
           Condition %in%  c("Model bias (large)")) %>% mutate(Estimate_perc = joint*100)


mbl_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") + 
  labs(y = "CSRM (%)", x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 10)) +
  theme(
    axis.title.y = element_text(size = 16),    
    axis.text = element_text(size = 14),       
    strip.text = element_text(size = 16),  
    legend.position = "none"
  )


boxplot = model_distribution_sup2 %>%
  filter(Domain == "Female" &                                           # Female
           Condition %in% c("Model bias (large) X small sample") &
           Estimate %in% c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM")) %>% mutate(Estimate_perc = Wo*100)

boxplot$Estimate = droplevels(boxplot$Estimate)
boxplot$Estimate = factor(boxplot$Estimate,
                          levels = c("CV_Vdesign", "RRMSE_GMSE", "RRMSE_PEM"))


benchmark_MSE = joint_distribution_sup2 %>%
  filter(Estimate == "RRMSE_joint" & Domain == "Female" & 
           Condition %in%  c("Model bias (large) X small sample")) %>% mutate(Estimate_perc = joint*100)


benchmark_var = joint_distribution_sup2 %>%
  filter(Estimate == "CV_joint" & Domain == "Female" & 
           Condition %in%  c("Model bias (large) X small sample")) %>% mutate(Estimate_perc = joint*100)


mblss_plot = ggplot(boxplot, aes(x = Estimate, y = Estimate_perc, fill = Estimate)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = custom_colors_female) +
  geom_hline(data = benchmark_MSE,
             aes(yintercept = Estimate_perc),
             linetype = "solid", alpha = 1, color = "orange",
             linewidth = 1.25) +
  geom_hline(data = benchmark_var,
             aes(yintercept = Estimate_perc),
             linetype = "dotted", alpha = 1, color = "orange",
             linewidth = 1.25) +
  facet_wrap(~ Condition, scales = "free") +  
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 10)) +
  theme(
    axis.title.y = element_text(size = 16),   
    axis.text = element_text(size = 14),  
    strip.text = element_text(size = 14), 
    legend.position = "none"
  )



combined_plot = (design_bias_plot | mbs_plot) /
  (mbl_plot | mblss_plot)

combined_plot