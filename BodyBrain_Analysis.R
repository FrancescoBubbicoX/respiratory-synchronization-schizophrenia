######################################################
#                                                    #
#            BODY-BRAIN-BEHAVIOUR ANALYSIS           #
#                                                    #
######################################################

# Written by Francesco Bubbico 
# Last updated: June 2026

library(readxl)
library(readr)
library(dplyr)
library(nlme)
library(emmeans)
library(ggplot2)
library(tidyr)
library(ggrain)
library(tibble)

# Load custom analysis functions
source("C:/Users/onesh/OneDrive/Desktop/Open data/Code/BodyBrain_Analysis_Functions.R")

# Load respiratory synchronization dataset
data_respSync <- read_excel(
  "C:/Users/onesh/OneDrive/Desktop/Open data/Data/Final datasets/RespSync_FinalDataset.xlsx"
)

# Prepare variables
data_respSync <- data_respSync %>%
  rename(
    subjectID = ID,
    condition = Condition
  ) %>%
  mutate(
    group = dplyr::recode(
      group,
      "Controls" = "controls",
      "Patients" = "patients"
    ),
    condition = dplyr::recode(
      condition,
      "Resp1_LP"  = "lp",
      "Resp1_HP"  = "hp",
      "Resp2_Pat" = "pat",
      "Resp2_Hea" = "hea"
    ),
    subjectID = factor(subjectID),
    condition = factor(condition),
    group = factor(group)
  )

# Create condition-specific datasets
data_respSync_hplp <- data_respSync %>%
  filter(condition %in% c("hp", "lp"))

data_respSync_heapat <- data_respSync %>%
  filter(condition %in% c("hea", "pat"))

# Set contrast 
options(contrasts = c("contr.sum", "contr.poly"))

###############################
#  Step 1:                    #   
#  Global effects analysis    #      
###############################

#######################################################################
#                   --- 1.1: Predictability ---                       #
#######################################################################

# Define condition setup
cond_setup_pred <- make_condition_setup(
  cond_levels = c("hp", "lp"),
  cond_labels = c("High pred.", "Low pred.")
)

# Define TE input file
te_file_pred <- "C:/Users/onesh/OneDrive/Desktop/Open data/Data/Final datasets/TE_electrode_hp_lp.csv"

# Prepare TE-behaviour dataset
data_pred <- prepare_te_electrode_data(
  te_file = te_file_pred,
  data_respSync = data_respSync_hplp,
  cond_setup = cond_setup_pred
)

# Behavioural and survey outcome variables
performance_vars_global <- c(
  "WCC_mean_maxCorr",
  "WCC_mean_maxLag_sec",
  "PLV",
  "mean_relative_phase_deg",
  "Survey1.RESP",
  "Survey2.RESP",
  "Survey3.RESP",
  "Survey4.RESP",
  "Survey5.RESP"
)

# Test whether global TE predicts behavioural and survey outcomes
global_all_results_pred <- run_global_te_behaviour(
  data = data_pred,
  cond_setup = cond_setup_pred,
  performance_vars_global = performance_vars_global
)

# Print significant TE-behaviour associations
cat("\n===== SIGNIFICANT GLOBAL TE-BEHAVIOUR ASSOCIATIONS: PREDICTABILITY =====\n")

global_all_results_pred %>%
  dplyr::filter(p_fdr < 0.05) %>%
  dplyr::arrange(question, performance_var, p_fdr) %>%
  dplyr::mutate(
    beta  = round(beta, 3),
    se    = round(se, 3),
    t     = round(t, 3),
    p     = round(p, 4),
    p_fdr = round(p_fdr, 4)
  ) %>%
  print(n = Inf)

# Extract significant global TE predictors
sig_predictors_pred <- global_all_results_pred %>%
  dplyr::filter(p_fdr < 0.05) %>%
  dplyr::distinct(interaction, metricVar)

cat("\n===== SELECTED GLOBAL TE PREDICTORS: PREDICTABILITY =====\n")
print(sig_predictors_pred, n = Inf)

# Keep only behaviourally relevant TE predictors
data_pred_selected <- data_pred %>%
  dplyr::semi_join(
    sig_predictors_pred,
    by = c("interaction", "metricVar")
  )

# Test Condition × Group effects for selected TE predictors
global_out_selected_pred <- run_global_lmms(
  data = data_pred_selected,
  cond_setup = cond_setup_pred,
  output_dir = NULL,
  prefix = "GLOBAL_Predictability_selected"
)

# Extract outputs
results_global_selected_pred <- global_out_selected_pred$all_results
sig_fdr_selected_pred        <- global_out_selected_pred$sig_fdr
posthoc_selected_pred        <- global_out_selected_pred$posthoc_fdr

# Print significant effects
cat("\n===== SIGNIFICANT GLOBAL EFFECTS: PREDICTABILITY =====\n")

sig_fdr_selected_pred %>%
  dplyr::mutate(
    F_value       = round(F_value, 3),
    p             = round(p, 4),
    p_fdr_by_type = round(p_fdr_by_type, 4)
  ) %>%
  print(n = Inf)

# Print post-hoc results
cat("\n===== POST-HOC RESULTS: PREDICTABILITY =====\n")

posthoc_selected_pred %>%
  dplyr::mutate(
    estimate = round(estimate, 4),
    SE       = round(SE, 4),
    t.ratio  = round(t.ratio, 3),
    p.value  = round(p.value, 4)
  ) %>%
  print(n = Inf)

# Output directory for body-brain-behaviour plots
plot_dir_bbb <- "C:/Users/onesh/OneDrive/Desktop/Open data/Plots/bodybrain_behavioural"
if (!dir.exists(plot_dir_bbb)) {
  dir.create(plot_dir_bbb, recursive = TRUE)
}


### Global TE-behaviour plots ###

# Respiration → Alpha1: TE association with WCC peak correlation
p_resp_alpha1 <- plot_global_te_behaviour(
  data_global = data_pred,
  interaction_pick = "resp_alpha1",
  metric_pick = "te_lin",
  performance_var = "WCC_mean_maxCorr",
  cond_setup = cond_setup_pred,
  effect_type = "group"
)

print(p_resp_alpha1$plot)

ggplot2::ggsave(
  file.path(plot_dir_bbb, "global_resp_alpha1_WCC_peakCorr_group_pred.png"),
  p_resp_alpha1$plot,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

# Post-hoc simple slopes by group
m_resp_alpha1 <- p_resp_alpha1$model

posthoc_resp_alpha1_group <- emmeans::emtrends(
  m_resp_alpha1,
  ~ group,
  var = "TE_z"
)

summary(posthoc_resp_alpha1_group)
pairs(posthoc_resp_alpha1_group)

# Raincloud plot of global TE values
p_global_resp_alpha1 <- plot_raincloud_2x2(
  data_cluster = data_pred %>%
    dplyr::filter(
      electrode == "GLOBAL",
      interaction == "resp_alpha1",
      metricVar == "te_lin"
    ) %>%
    dplyr::rename(TE = metricValue),
  plot_title = "Respiration → Alpha1 Global (Linear)",
  ylab = "Global TE",
  xlab = "Predictability condition",
  condition_labels = c(
    hp = "High pred.",
    lp = "Low pred."
  ),
  group_labels = c(
    controls = "Controls",
    patients = "Patients"
  )
)

print(p_global_resp_alpha1)

ggplot2::ggsave(
  file.path(plot_dir_bbb, "global_resp_alpha1_raincloud_pred.png"),
  p_global_resp_alpha1,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)


# Respiration → Gamma: TE association with WCC peak lag
p_resp_gamma <- plot_global_te_behaviour(
  data_global = data_pred,
  interaction_pick = "resp_gamma",
  metric_pick = "te_knn",
  performance_var = "WCC_mean_maxLag_sec",
  cond_setup = cond_setup_pred,
  effect_type = "group"
)

print(p_resp_gamma$plot)

ggplot2::ggsave(
  file.path(plot_dir_bbb, "global_resp_gamma_WCC_peakLag_group_pred.png"),
  p_resp_gamma$plot,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

# Post-hoc simple slopes by group
m_resp_gamma <- p_resp_gamma$model

posthoc_resp_gamma_group <- emmeans::emtrends(
  m_resp_gamma,
  ~ group,
  var = "TE_z"
)

summary(posthoc_resp_gamma_group)
pairs(posthoc_resp_gamma_group)

# Raincloud plot of global TE values
p_global_resp_gamma <- plot_raincloud_2x2(
  data_cluster = data_pred %>%
    dplyr::filter(
      electrode == "GLOBAL",
      interaction == "resp_gamma",
      metricVar == "te_knn"
    ) %>%
    dplyr::rename(TE = metricValue),
  plot_title = "Respiration → Gamma Global (Nonlinear)",
  ylab = "Global TE",
  xlab = "Predictability condition",
  condition_labels = c(
    hp = "High pred.",
    lp = "Low pred."
  ),
  group_labels = c(
    controls = "Controls",
    patients = "Patients"
  )
)

print(p_global_resp_gamma)

ggplot2::ggsave(
  file.path(plot_dir_bbb, "global_resp_gamma_raincloud_pred.png"),
  p_global_resp_gamma,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)


# Alpha1 → Respiration: TE association with WCC peak correlation
p_alpha1_resp <- plot_global_te_behaviour(
  data_global = data_pred,
  interaction_pick = "alpha1_resp",
  metric_pick = "te_knn",
  performance_var = "WCC_mean_maxCorr",
  cond_setup = cond_setup_pred,
  effect_type = "condition"
)

print(p_alpha1_resp$plot)

ggplot2::ggsave(
  file.path(plot_dir_bbb, "global_alpha1_resp_WCC_peakCorr_condition_pred.png"),
  p_alpha1_resp$plot,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

# Post-hoc simple slopes by condition
m_alpha1_resp <- p_alpha1_resp$model

posthoc_alpha1_resp_condition <- emmeans::emtrends(
  m_alpha1_resp,
  ~ condition,
  var = "TE_z"
)

summary(posthoc_alpha1_resp_condition)
pairs(posthoc_alpha1_resp_condition)


# Alpha → Respiration: TE association with WCC peak correlation
p_alpha_resp <- plot_global_te_behaviour(
  data_global = data_pred,
  interaction_pick = "alpha_resp",
  metric_pick = "te_knn",
  performance_var = "WCC_mean_maxCorr",
  cond_setup = cond_setup_pred,
  effect_type = "condition"
)

print(p_alpha_resp$plot)

ggplot2::ggsave(
  file.path(plot_dir_bbb, "global_alpha_resp_WCC_peakCorr_condition_pred.png"),
  p_alpha_resp$plot,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

# Post-hoc simple slopes by condition
m_alpha_resp <- p_alpha_resp$model

posthoc_alpha_resp_condition <- emmeans::emtrends(
  m_alpha_resp,
  ~ condition,
  var = "TE_z"
)

summary(posthoc_alpha_resp_condition)
pairs(posthoc_alpha_resp_condition)


#######################################################################
#                         --- 1.2: Bias ---                           #
#######################################################################

# Define condition setup
cond_setup_bias <- make_condition_setup(
  cond_levels = c("hea", "pat"),
  cond_labels = c("Healthy bias", "Patient bias")
)

# Define TE input file
te_file_bias <- "C:/Users/onesh/OneDrive/Desktop/Open data/Data/Final datasets/TE_electrode_hea_pat.csv"

# Prepare TE-behaviour dataset
data_bias <- prepare_te_electrode_data(
  te_file = te_file_bias,
  data_respSync = data_respSync_heapat,
  cond_setup = cond_setup_bias
)

# Behavioural and survey outcome variables
performance_vars_global <- c(
  "WCC_mean_maxCorr",
  "WCC_mean_maxLag_sec",
  "PLV",
  "mean_relative_phase_deg",
  "Survey1.RESP",
  "Survey2.RESP",
  "Survey3.RESP",
  "Survey4.RESP",
  "Survey5.RESP"
)

# Test whether global TE predicts behavioural and survey outcomes
global_all_results_bias <- run_global_te_behaviour(
  data = data_bias,
  cond_setup = cond_setup_bias,
  performance_vars_global = performance_vars_global
)

# Print significant TE-behaviour associations
cat("\n===== SIGNIFICANT GLOBAL TE-BEHAVIOUR ASSOCIATIONS: BIAS =====\n")

global_all_results_bias %>%
  dplyr::filter(p_fdr < 0.05) %>%
  dplyr::arrange(question, performance_var, p_fdr) %>%
  dplyr::mutate(
    beta  = round(beta, 3),
    se    = round(se, 3),
    t     = round(t, 3),
    p     = round(p, 4),
    p_fdr = round(p_fdr, 4)
  ) %>%
  print(n = Inf)

# Extract significant global TE predictors
sig_predictors_bias <- global_all_results_bias %>%
  dplyr::filter(p_fdr < 0.05) %>%
  dplyr::distinct(interaction, metricVar)

cat("\n===== SELECTED GLOBAL TE PREDICTORS: BIAS =====\n")
print(sig_predictors_bias, n = Inf)

# Keep only behaviourally relevant TE predictors
data_bias_selected <- data_bias %>%
  dplyr::semi_join(
    sig_predictors_bias,
    by = c("interaction", "metricVar")
  )

# Test Condition × Group effects for selected TE predictors
global_out_selected_bias <- run_global_lmms(
  data = data_bias_selected,
  cond_setup = cond_setup_bias,
  output_dir = NULL,
  prefix = "GLOBAL_Bias_selected"
)

# Extract outputs
results_global_selected_bias <- global_out_selected_bias$all_results
sig_fdr_selected_bias        <- global_out_selected_bias$sig_fdr
posthoc_selected_bias        <- global_out_selected_bias$posthoc_fdr

# Print significant effects
cat("\n===== SIGNIFICANT GLOBAL EFFECTS: BIAS =====\n")

sig_fdr_selected_bias %>%
  dplyr::mutate(
    F_value       = round(F_value, 3),
    p             = round(p, 4),
    p_fdr_by_type = round(p_fdr_by_type, 4)
  ) %>%
  print(n = Inf)

# Print post-hoc results
cat("\n===== POST-HOC RESULTS: BIAS =====\n")

posthoc_selected_bias %>%
  dplyr::mutate(
    estimate = round(estimate, 4),
    SE       = round(SE, 4),
    t.ratio  = round(t.ratio, 3),
    p.value  = round(p.value, 4)
  ) %>%
  print(n = Inf)


### Global TE-behaviour plots ###

# Respiration → Theta: TE association with IOS
p_resp_theta_survey5 <- plot_global_te_behaviour(
  data_global = data_bias,
  interaction_pick = "resp_theta",
  metric_pick = "te_lin",
  performance_var = "Survey5.RESP",
  cond_setup = cond_setup_bias,
  effect_type = "condition"
)

print(p_resp_theta_survey5$plot)

ggplot2::ggsave(
  file.path(plot_dir_bbb, "global_resp_theta_Survey5_condition_bias.png"),
  p_resp_theta_survey5$plot,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

# Post-hoc simple slopes by bias condition
m_resp_theta_survey5 <- p_resp_theta_survey5$model

posthoc_resp_theta_condition <- emmeans::emtrends(
  m_resp_theta_survey5,
  ~ condition,
  var = "TE_z"
)

summary(posthoc_resp_theta_condition)
pairs(posthoc_resp_theta_condition)


###############################
#  Step 2:                    #   
#  Cluster analysis           #      
###############################

# Plot folders
plot_dir <- "C:/Users/onesh/OneDrive/Desktop/Open data/Plots/bodybrain_behavioural"
plot_dir_bbb <- "C:/Users/onesh/OneDrive/Desktop/Open data/Plots/bodybrain_behavioural"


#######################################################################
#                   --- 2.1: Predictability ---                       #
#######################################################################

# Load and prepare cluster-level TE data
data_cluster_pred <- read_csv(
  "C:/Users/onesh/OneDrive/Desktop/Open data/Data/Final datasets/TE_cluster_hp_lp.csv",
  show_col_types = FALSE
) %>%
  mutate(
    subjectID = ifelse(as.character(subjectID) == "031", "31", as.character(subjectID)),
    subjectID = factor(subjectID),
    group     = relevel(factor(group), ref = "controls"),
    condition = relevel(factor(condition), ref = "hp"),
    interaction = factor(interaction),
    metricVar   = factor(metricVar),
    clusterID   = factor(clusterID),
    TE = metric_mean
  ) %>%
  left_join(
    data_respSync_hplp,
    by = c("subjectID", "condition", "group")
  )

# Plot labels
group_map <- c(
  controls = "Controls",
  patients = "Patients"
)

cond_levels <- c("hp", "lp")
cond_labels <- c("High pred.", "Low pred.")

# Significant clusters from permutation analysis
cluster_info_sig <- data_cluster_pred %>%
  filter(isSigCluster == TRUE) %>%
  distinct(clusterID, clusterP, clusterMass, electrodes, nElectrodes) %>%
  arrange(clusterP)

cluster_info_sig


# Alpha1 → RR cluster
cluster_pick <- "te_lin_alpha1_rr_clu_1"

data_cluster <- data_cluster_pred %>%
  filter(clusterID == cluster_pick, !is.na(TE))

ttl <- make_plot_title(data_cluster, cluster_pick)
file_stem <- paste(paste(cond_levels, collapse = "_"), cluster_pick, sep = "_")

p_rain <- plot_raincloud_2x2(
  data_cluster = data_cluster,
  clusterID = cluster_pick,
  plot_title = ttl,
  condition_order = cond_levels,
  condition_labels = stats::setNames(cond_labels, cond_levels),
  group_labels = group_map
)

print(p_rain)

ggplot2::ggsave(
  file.path(plot_dir, paste0(file_stem, "_raincloud.png")),
  p_rain,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

cluster_alpha1_rr <- run_cluster_lmm_analysis(
  data = data_cluster_pred,
  cluster_pick = cluster_pick
)


# RR → Beta2 cluster
cluster_pick <- "te_lin_rr_beta2_clu_1"

data_cluster <- data_cluster_pred %>%
  filter(clusterID == cluster_pick, !is.na(TE))

ttl <- make_plot_title(data_cluster, cluster_pick)
file_stem <- paste(paste(cond_levels, collapse = "_"), cluster_pick, sep = "_")

p_rain <- plot_raincloud_2x2(
  data_cluster = data_cluster,
  clusterID = cluster_pick,
  plot_title = ttl,
  condition_order = cond_levels,
  condition_labels = stats::setNames(cond_labels, cond_levels),
  group_labels = group_map
)

print(p_rain)

ggplot2::ggsave(
  file.path(plot_dir, paste0(file_stem, "_raincloud.png")),
  p_rain,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

cluster_rr_beta2 <- run_cluster_lmm_analysis(
  data = data_cluster_pred,
  cluster_pick = cluster_pick
)


# Respiration → Gamma cluster
cluster_pick <- "te_lin_resp_gamma_clu_1"

data_cluster <- data_cluster_pred %>%
  filter(clusterID == cluster_pick, !is.na(TE))

ttl <- make_plot_title(data_cluster, cluster_pick)
file_stem <- paste(paste(cond_levels, collapse = "_"), cluster_pick, sep = "_")

p_rain <- plot_raincloud_2x2(
  data_cluster = data_cluster,
  clusterID = cluster_pick,
  plot_title = ttl,
  condition_order = cond_levels,
  condition_labels = stats::setNames(cond_labels, cond_levels),
  group_labels = group_map
)

print(p_rain)

ggplot2::ggsave(
  file.path(plot_dir, paste0(file_stem, "_raincloud.png")),
  p_rain,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

cluster_resp_gamma <- run_cluster_lmm_analysis(
  data = data_cluster_pred,
  cluster_pick = cluster_pick
)


# Theta → Respiration cluster
cluster_pick <- "te_lin_theta_resp_clu_1"

data_cluster <- data_cluster_pred %>%
  filter(clusterID == cluster_pick, !is.na(TE))

ttl <- make_plot_title(data_cluster, cluster_pick)
file_stem <- paste(paste(cond_levels, collapse = "_"), cluster_pick, sep = "_")

p_rain <- plot_raincloud_2x2(
  data_cluster = data_cluster,
  clusterID = cluster_pick,
  plot_title = ttl,
  condition_order = cond_levels,
  condition_labels = stats::setNames(cond_labels, cond_levels),
  group_labels = group_map
)

print(p_rain)

ggplot2::ggsave(
  file.path(plot_dir, paste0(file_stem, "_raincloud.png")),
  p_rain,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

cluster_theta_resp <- run_cluster_lmm_analysis(
  data = data_cluster_pred,
  cluster_pick = cluster_pick
)


### Predictability: cluster-level TE-behaviour delta analysis ###

# Behavioural and survey outcome variables
performance_vars_delta <- c(
  "WCC_mean_maxCorr",
  "WCC_mean_maxLag_sec",
  "PLV",
  "mean_relative_phase_deg",
  "Survey1.RESP",
  "Survey2.RESP",
  "Survey3.RESP",
  "Survey4.RESP",
  "Survey5.RESP"
)

# Prepare significant-cluster dataset
data_sig_pred <- data_cluster_pred %>%
  dplyr::filter(
    isSigCluster == TRUE,
    condition %in% cond_setup_pred$levels
  ) %>%
  dplyr::mutate(
    TE = as.numeric(metric_mean),
    subjectID = factor(subjectID),
    group = relevel(factor(group), ref = "controls"),
    condition = relevel(
      factor(condition, levels = cond_setup_pred$levels),
      ref = cond_setup_pred$ref
    ),
    clusterID = factor(clusterID)
  ) %>%
  droplevels()

# Run delta TE-behaviour models
delta_results_pred <- data.frame()

for (performance_var in performance_vars_delta) {
  
  tmp_data <- data_sig_pred %>%
    dplyr::mutate(performance = as.numeric(.data[[performance_var]]))
  
  for (clu in unique(tmp_data$clusterID)) {
    
    dat_clu <- tmp_data %>%
      dplyr::filter(clusterID == clu) %>%
      droplevels()
    
    dup_check <- dat_clu %>%
      dplyr::count(subjectID, condition)
    
    if (any(dup_check$n > 1)) {
      warning(paste("Duplicate rows detected in", clu))
    }
    
    delta_dat <- dat_clu %>%
      dplyr::select(
        subjectID,
        group,
        clusterID,
        interaction,
        metricVar,
        electrodes,
        condition,
        TE,
        performance
      ) %>%
      tidyr::pivot_wider(
        names_from = condition,
        values_from = c(TE, performance)
      ) %>%
      dplyr::mutate(
        delta_TE =
          .data[[paste0("TE_", cond_setup_pred$cond2)]] -
          .data[[paste0("TE_", cond_setup_pred$cond1)]],
        
        delta_performance =
          .data[[paste0("performance_", cond_setup_pred$cond2)]] -
          .data[[paste0("performance_", cond_setup_pred$cond1)]],
        
        delta_TE_z = as.numeric(scale(delta_TE)),
        delta_performance_z = as.numeric(scale(delta_performance)),
        group = relevel(factor(group), ref = "controls")
      ) %>%
      droplevels()
    
    m_general <- tryCatch(
      lm(delta_performance_z ~ delta_TE_z + group, data = delta_dat),
      error = function(e) NULL
    )
    
    m_group <- tryCatch(
      lm(delta_performance_z ~ delta_TE_z * group, data = delta_dat),
      error = function(e) NULL
    )
    
    delta_results_pred <- rbind(
      delta_results_pred,
      
      extract_lm_term(
        m_general,
        "^delta_TE_z$",
        "general_delta_association",
        unique(delta_dat$clusterID),
        unique(delta_dat$interaction),
        unique(delta_dat$metricVar),
        unique(delta_dat$electrodes),
        performance_var
      ),
      
      extract_lm_term(
        m_group,
        "delta_TE_z:group",
        "group_interaction_delta",
        unique(delta_dat$clusterID),
        unique(delta_dat$interaction),
        unique(delta_dat$metricVar),
        unique(delta_dat$electrodes),
        performance_var
      )
    )
  }
}

# FDR correction
delta_results_pred <- delta_results_pred %>%
  dplyr::group_by(question, performance_var) %>%
  dplyr::mutate(p_fdr = p.adjust(p, method = "fdr")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(question, performance_var, p)

# Print FDR-significant delta results
delta_results_pred %>%
  dplyr::filter(p_fdr < 0.05) %>%
  dplyr::arrange(question, performance_var, p_fdr) %>%
  dplyr::mutate(
    beta  = round(beta, 3),
    se    = round(se, 3),
    t     = round(t, 3),
    p     = round(p, 4),
    p_fdr = round(p_fdr, 4)
  ) %>%
  print(n = Inf)


# Selected cluster delta plot

delta_rr_beta2_wcc <- plot_delta_group_interaction(
  data_sig = data_sig_pred,
  cluster_pick = "te_lin_rr_beta2_clu_1",
  performance_var = "WCC_mean_maxCorr",
  cond_setup = cond_setup_pred
)

print(delta_rr_beta2_wcc$plot)

# Simple slopes by group
posthoc_delta_rr_beta2_group <- emmeans::emtrends(
  delta_rr_beta2_wcc$model,
  ~ group,
  var = "delta_TE_z"
)

summary(posthoc_delta_rr_beta2_group)
pairs(posthoc_delta_rr_beta2_group)

ggplot2::ggsave(
  file.path(plot_dir_bbb, "delta_rr_beta2_WCC_peak_correlation_pred.png"),
  delta_rr_beta2_wcc$plot,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)


#######################################################################
#                         --- 2.2: Bias ---                           #
#######################################################################

# Load and prepare cluster-level TE data
data_cluster_bias <- read_csv(
  "C:/Users/onesh/OneDrive/Desktop/Open data/Data/Final datasets/TE_cluster_hea_pat.csv",
  show_col_types = FALSE
) %>%
  mutate(
    subjectID = ifelse(as.character(subjectID) == "031", "31", as.character(subjectID)),
    subjectID = factor(subjectID),
    group     = relevel(factor(group), ref = "controls"),
    condition = relevel(factor(condition), ref = "hea"),
    interaction = factor(interaction),
    metricVar   = factor(metricVar),
    clusterID   = factor(clusterID),
    TE = metric_mean
  ) %>%
  left_join(
    data_respSync_heapat,
    by = c("subjectID", "condition", "group")
  )

# Plot labels
cond_levels <- c("hea", "pat")
cond_labels <- c("Healthy bias", "Patient bias")

# Significant clusters from permutation analysis
cluster_info_sig_bias <- data_cluster_bias %>%
  dplyr::filter(isSigCluster == TRUE) %>%
  dplyr::distinct(clusterID, clusterP, clusterMass, electrodes, nElectrodes) %>%
  dplyr::arrange(clusterP)

cluster_info_sig_bias


# Beta1 → Respiration cluster
cluster_pick <- "te_lin_beta1_resp_clu_1"

data_cluster <- data_cluster_bias %>%
  dplyr::filter(clusterID == cluster_pick, !is.na(TE))

ttl <- make_plot_title(data_cluster, cluster_pick)
file_stem <- paste(paste(cond_levels, collapse = "_"), cluster_pick, sep = "_")

p_rain <- plot_raincloud_2x2(
  data_cluster = data_cluster,
  clusterID = cluster_pick,
  plot_title = ttl,
  condition_order = cond_levels,
  condition_labels = stats::setNames(cond_labels, cond_levels),
  group_labels = group_map
)

print(p_rain)

ggplot2::ggsave(
  file.path(plot_dir_bbb, paste0(file_stem, "_raincloud.png")),
  p_rain,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

cluster_beta1_resp <- run_cluster_lmm_analysis(
  data = data_cluster_bias,
  cluster_pick = cluster_pick
)


# RR → Alpha cluster
cluster_pick <- "te_knn_rr_alpha_clu_2"

data_cluster <- data_cluster_bias %>%
  dplyr::filter(clusterID == cluster_pick, !is.na(TE))

ttl <- make_plot_title(data_cluster, cluster_pick)
file_stem <- paste(paste(cond_levels, collapse = "_"), cluster_pick, sep = "_")

p_rain <- plot_raincloud_2x2(
  data_cluster = data_cluster,
  clusterID = cluster_pick,
  plot_title = ttl,
  condition_order = cond_levels,
  condition_labels = stats::setNames(cond_labels, cond_levels),
  group_labels = group_map
)

print(p_rain)

ggplot2::ggsave(
  file.path(plot_dir_bbb, paste0(file_stem, "_raincloud.png")),
  p_rain,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

cluster_rr_alpha <- run_cluster_lmm_analysis(
  data = data_cluster_bias,
  cluster_pick = cluster_pick
)

# RR → Beta2 cluster
cluster_pick <- "te_knn_rr_beta2_clu_2"

data_cluster <- data_cluster_bias %>%
  dplyr::filter(clusterID == cluster_pick, !is.na(TE))

ttl <- make_plot_title(data_cluster, cluster_pick)
file_stem <- paste(paste(cond_levels, collapse = "_"), cluster_pick, sep = "_")

p_rain <- plot_raincloud_2x2(
  data_cluster = data_cluster,
  clusterID = cluster_pick,
  plot_title = ttl,
  condition_order = cond_levels,
  condition_labels = stats::setNames(cond_labels, cond_levels),
  group_labels = group_map
)

print(p_rain)

ggplot2::ggsave(
  file.path(plot_dir_bbb, paste0(file_stem, "_raincloud.png")),
  p_rain,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

cluster_rr_beta2_bias <- run_cluster_lmm_analysis(
  data = data_cluster_bias,
  cluster_pick = cluster_pick
)

# RR → Beta1 cluster
cluster_pick <- "te_lin_rr_beta1_clu_1"

data_cluster <- data_cluster_bias %>%
  dplyr::filter(clusterID == cluster_pick, !is.na(TE))

ttl <- make_plot_title(data_cluster, cluster_pick)
file_stem <- paste(paste(cond_levels, collapse = "_"), cluster_pick, sep = "_")

p_rain <- plot_raincloud_2x2(
  data_cluster = data_cluster,
  clusterID = cluster_pick,
  plot_title = ttl,
  condition_order = cond_levels,
  condition_labels = stats::setNames(cond_labels, cond_levels),
  group_labels = group_map
)

print(p_rain)

ggplot2::ggsave(
  file.path(plot_dir_bbb, paste0(file_stem, "_raincloud.png")),
  p_rain,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

cluster_rr_beta1 <- run_cluster_lmm_analysis(
  data = data_cluster_bias,
  cluster_pick = cluster_pick
)

# Alpha1 → Respiration cluster
cluster_pick <- "te_lin_alpha1_resp_clu_1"

data_cluster <- data_cluster_bias %>%
  dplyr::filter(clusterID == cluster_pick, !is.na(TE))

ttl <- make_plot_title(data_cluster, cluster_pick)
file_stem <- paste(paste(cond_levels, collapse = "_"), cluster_pick, sep = "_")

p_rain <- plot_raincloud_2x2(
  data_cluster = data_cluster,
  clusterID = cluster_pick,
  plot_title = ttl,
  condition_order = cond_levels,
  condition_labels = stats::setNames(cond_labels, cond_levels),
  group_labels = group_map
)

print(p_rain)

ggplot2::ggsave(
  file.path(plot_dir_bbb, paste0(file_stem, "_raincloud.png")),
  p_rain,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

cluster_alpha1_resp_bias <- run_cluster_lmm_analysis(
  data = data_cluster_bias,
  cluster_pick = cluster_pick
)

# Alpha2 → RR cluster
cluster_pick <- "te_knn_alpha2_rr_clu_1"

data_cluster <- data_cluster_bias %>%
  dplyr::filter(clusterID == cluster_pick, !is.na(TE))

ttl <- make_plot_title(data_cluster, cluster_pick)
file_stem <- paste(paste(cond_levels, collapse = "_"), cluster_pick, sep = "_")

p_rain <- plot_raincloud_2x2(
  data_cluster = data_cluster,
  clusterID = cluster_pick,
  plot_title = ttl,
  condition_order = cond_levels,
  condition_labels = stats::setNames(cond_labels, cond_levels),
  group_labels = group_map
)

print(p_rain)

ggplot2::ggsave(
  file.path(plot_dir_bbb, paste0(file_stem, "_raincloud.png")),
  p_rain,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

cluster_alpha2_rr <- run_cluster_lmm_analysis(
  data = data_cluster_bias,
  cluster_pick = cluster_pick
)

# Respiration → Alpha2 cluster
cluster_pick <- "te_lin_resp_alpha2_clu_1"

data_cluster <- data_cluster_bias %>%
  dplyr::filter(clusterID == cluster_pick, !is.na(TE))

ttl <- make_plot_title(data_cluster, cluster_pick)
file_stem <- paste(paste(cond_levels, collapse = "_"), cluster_pick, sep = "_")

p_rain <- plot_raincloud_2x2(
  data_cluster = data_cluster,
  clusterID = cluster_pick,
  plot_title = ttl,
  condition_order = cond_levels,
  condition_labels = stats::setNames(cond_labels, cond_levels),
  group_labels = group_map
)

print(p_rain)

ggplot2::ggsave(
  file.path(plot_dir_bbb, paste0(file_stem, "_raincloud.png")),
  p_rain,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

cluster_resp_alpha2 <- run_cluster_lmm_analysis(
  data = data_cluster_bias,
  cluster_pick = cluster_pick
)


### Bias: cluster-level TE-behaviour delta analysis ###

data_sig_bias <- data_cluster_bias %>%
  dplyr::filter(
    isSigCluster == TRUE,
    condition %in% cond_setup_bias$levels
  ) %>%
  dplyr::mutate(
    TE = as.numeric(metric_mean),
    subjectID = factor(subjectID),
    group = relevel(factor(group), ref = "controls"),
    condition = relevel(
      factor(condition, levels = cond_setup_bias$levels),
      ref = cond_setup_bias$ref
    ),
    clusterID = factor(clusterID)
  ) %>%
  droplevels()

delta_results_bias <- data.frame()

for (performance_var in performance_vars_delta) {
  
  tmp_data <- data_sig_bias %>%
    dplyr::mutate(performance = as.numeric(.data[[performance_var]]))
  
  for (clu in unique(tmp_data$clusterID)) {
    
    dat_clu <- tmp_data %>%
      dplyr::filter(clusterID == clu) %>%
      droplevels()
    
    dup_check <- dat_clu %>%
      dplyr::count(subjectID, condition)
    
    if (any(dup_check$n > 1)) {
      warning(paste("Duplicate rows detected in", clu))
    }
    
    delta_dat <- dat_clu %>%
      dplyr::select(
        subjectID, group, clusterID, interaction, metricVar,
        electrodes, condition, TE, performance
      ) %>%
      tidyr::pivot_wider(
        names_from = condition,
        values_from = c(TE, performance)
      ) %>%
      dplyr::mutate(
        delta_TE =
          .data[[paste0("TE_", cond_setup_bias$cond2)]] -
          .data[[paste0("TE_", cond_setup_bias$cond1)]],
        delta_performance =
          .data[[paste0("performance_", cond_setup_bias$cond2)]] -
          .data[[paste0("performance_", cond_setup_bias$cond1)]],
        delta_TE_z = as.numeric(scale(delta_TE)),
        delta_performance_z = as.numeric(scale(delta_performance)),
        group = relevel(factor(group), ref = "controls")
      ) %>%
      droplevels()
    
    m_general <- tryCatch(
      lm(delta_performance_z ~ delta_TE_z + group, data = delta_dat),
      error = function(e) NULL
    )
    
    m_group <- tryCatch(
      lm(delta_performance_z ~ delta_TE_z * group, data = delta_dat),
      error = function(e) NULL
    )
    
    delta_results_bias <- rbind(
      delta_results_bias,
      extract_lm_term(
        m_general,
        "^delta_TE_z$",
        "general_delta_association",
        unique(delta_dat$clusterID),
        unique(delta_dat$interaction),
        unique(delta_dat$metricVar),
        unique(delta_dat$electrodes),
        performance_var
      ),
      extract_lm_term(
        m_group,
        "delta_TE_z:group",
        "group_interaction_delta",
        unique(delta_dat$clusterID),
        unique(delta_dat$interaction),
        unique(delta_dat$metricVar),
        unique(delta_dat$electrodes),
        performance_var
      )
    )
  }
}

delta_results_bias <- delta_results_bias %>%
  dplyr::group_by(question, performance_var) %>%
  dplyr::mutate(p_fdr = p.adjust(p, method = "fdr")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(question, performance_var, p)

delta_results_bias %>%
  dplyr::filter(p_fdr < 0.05) %>%
  dplyr::arrange(question, performance_var, p_fdr) %>%
  dplyr::mutate(
    beta  = round(beta, 3),
    se    = round(se, 3),
    t     = round(t, 3),
    p     = round(p, 4),
    p_fdr = round(p_fdr, 4)
  ) %>%
  print(n = Inf)


delta_beta1_resp_ios <- plot_delta_group_interaction(
  data_sig = data_sig_bias,
  cluster_pick = "te_lin_beta1_resp_clu_1",
  performance_var = "Survey5.RESP",
  cond_setup = cond_setup_bias
)

print(delta_beta1_resp_ios$plot)

posthoc_delta_beta1_resp_group <- emmeans::emtrends(
  delta_beta1_resp_ios$model,
  ~ group,
  var = "delta_TE_z"
)

summary(posthoc_delta_beta1_resp_group)
pairs(posthoc_delta_beta1_resp_group)

ggplot2::ggsave(
  file.path(plot_dir_bbb, "delta_beta1_resp_Survey5_bias.png"),
  delta_beta1_resp_ios$plot,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)
