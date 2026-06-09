######################################################
#                                                    #
#      RESPIRATORY SYNCHRONIZATION TASK ANALYSIS     #
#                                                    #
######################################################

# Written by Francesco Bubbico 
# Last updated: June 2026

library(readxl)
library(writexl)
library(dplyr)
library(nlme)
library(emmeans)

############################
#  Step 1:                 #   
#  Create final dataset    #      
############################

# Define paths
path_info_survey <- "C:/Users/onesh/OneDrive/Desktop/Open data/Data/Final datasets/ParticipantInfoAndSurvey.xlsx"
path_resp_task   <- "C:/Users/onesh/OneDrive/Desktop/Open data/Data/Final datasets/RespiratorySync_dataset.csv"

output_path_xlsx <- "C:/Users/onesh/OneDrive/Desktop/Open data/Data/Final datasets/RespSync_FinalDataset.xlsx"

# Load participant information and survey data
info_and_survey <- read_excel(path_info_survey) %>%
  mutate(
    ID        = as.numeric(ID),
    Condition = as.character(Condition),
    group     = as.character(group)
  )

# Load respiratory synchronization task data
respTask <- read.csv(path_resp_task) %>%
  mutate(
    ID        = as.numeric(ID),
    Condition = as.character(Condition),
    group     = as.character(group)
  )


# Harmonise condition labels
respTask <- respTask %>%
  mutate(
    Condition = dplyr::recode(
      Condition,
      "RespSync_LP_data"  = "Resp1_LP",
      "RespSync_HP_data"  = "Resp1_HP",
      "RespSync_HEA_data" = "Resp2_Hea",
      "RespSync_PAT_data" = "Resp2_Pat"
    )
  )

# Merge respiratory task data with participant/survey data
dataset_respTask <- respTask %>%
  full_join(
    info_and_survey,
    by = c("ID", "Condition", "group")
  )

# Clean group labels
dataset_respTask <- dataset_respTask %>%
  mutate(
    group = dplyr::recode(
      group,
      "PAZIENTI"  = "Patients",
      "CONTROLLI" = "Controls"
    )
  )

# Save final dataset
write_xlsx(dataset_respTask, path = output_path_xlsx)


############################
#  Step 2:                 #   
#  Sample characteristics  #      
############################

participant_level <- dataset_respTask %>%
  distinct(ID, group, age, gender)

# Age: descriptive statistics
age_summary <- participant_level %>%
  group_by(group) %>%
  summarise(
    n    = sum(!is.na(age)),
    mean = round(mean(age, na.rm = TRUE), 2),
    sd   = round(sd(age, na.rm = TRUE), 2),
    min  = min(age, na.rm = TRUE),
    max  = max(age, na.rm = TRUE),
    .groups = "drop"
  )
age_summary

# Age: independent-samples t-test
age_ttest <- t.test(age ~ group, data = participant_level)
age_ttest

# Gender: frequency table and Fisher's exact test
gender_table  <- table(participant_level$group, participant_level$gender)
gender_table
gender_fisher <- fisher.test(gender_table)
gender_fisher


############################
#  Step 3:                 #   
#  Statistical analysis    #      
############################

# Prepare variables
dataset_respTask <- dataset_respTask %>%
  mutate(
    ID        = factor(ID),
    Condition = factor(Condition),
    group     = factor(group, levels = c("Controls", "Patients")),
    gender    = factor(gender),
    age       = as.numeric(age)
  )

options(contrasts = c("contr.sum", "contr.poly"))


# --- 3.1: Predictability models ---

respTask_predictability <- dataset_respTask %>%
  filter(
    Condition %in% c("Resp1_HP", "Resp1_LP"),
    group %in% c("Controls", "Patients")
  ) %>%
  mutate(
    Condition = droplevels(Condition),
    Condition = relevel(Condition, ref = "Resp1_HP"),
    group     = relevel(group, ref = "Controls")
  )

# PLV
m_pred_plv <- lme(
  PLV ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_predictability,
  na.action = na.exclude
)

summary(m_pred_plv)
anova(m_pred_plv, type = "marginal")
emmeans(m_pred_plv, ~ Condition * group)
pairs(emmeans(m_pred_plv, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_pred_plv, ~ group | Condition), adjust = "holm")

# Relative phase
m_pred_phase <- lme(
  mean_relative_phase_deg ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_predictability,
  na.action = na.exclude
)

summary(m_pred_phase)
anova(m_pred_phase, type = "marginal")
emmeans(m_pred_phase, ~ Condition * group)
pairs(emmeans(m_pred_phase, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_pred_phase, ~ group | Condition), adjust = "holm")

# WCC maximum correlation
m_pred_wcc_corr <- lme(
  WCC_mean_maxCorr ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_predictability,
  na.action = na.exclude
)

summary(m_pred_wcc_corr)
anova(m_pred_wcc_corr, type = "marginal")
emmeans(m_pred_wcc_corr, ~ Condition * group)
pairs(emmeans(m_pred_wcc_corr, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_pred_wcc_corr, ~ group | Condition), adjust = "holm")

# WCC lag
m_pred_wcc_lag <- lme(
  WCC_mean_maxLag_sec ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_predictability,
  na.action = na.exclude
)

summary(m_pred_wcc_lag)
anova(m_pred_wcc_lag, type = "marginal")
emmeans(m_pred_wcc_lag, ~ Condition * group)
pairs(emmeans(m_pred_wcc_lag, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_pred_wcc_lag, ~ group | Condition), adjust = "holm")

# Survey 1
m_pred_survey1 <- lme(
  Survey1.RESP ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_predictability,
  na.action = na.exclude
)

summary(m_pred_survey1)
anova(m_pred_survey1, type = "marginal")
emmeans(m_pred_survey1, ~ Condition * group)
pairs(emmeans(m_pred_survey1, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_pred_survey1, ~ group | Condition), adjust = "holm")

# Survey 2
m_pred_survey2 <- lme(
  Survey2.RESP ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_predictability,
  na.action = na.exclude
)

summary(m_pred_survey2)
anova(m_pred_survey2, type = "marginal")
emmeans(m_pred_survey2, ~ Condition * group)
pairs(emmeans(m_pred_survey2, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_pred_survey2, ~ group | Condition), adjust = "holm")

# Survey 3
m_pred_survey3 <- lme(
  Survey3.RESP ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_predictability,
  na.action = na.exclude
)

summary(m_pred_survey3)
anova(m_pred_survey3, type = "marginal")
emmeans(m_pred_survey3, ~ Condition * group)
pairs(emmeans(m_pred_survey3, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_pred_survey3, ~ group | Condition), adjust = "holm")

# Survey 4
m_pred_survey4 <- lme(
  Survey4.RESP ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_predictability,
  na.action = na.exclude
)

summary(m_pred_survey4)
anova(m_pred_survey4, type = "marginal")
emmeans(m_pred_survey4, ~ Condition * group)
pairs(emmeans(m_pred_survey4, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_pred_survey4, ~ group | Condition), adjust = "holm")

# Survey 5
m_pred_survey5 <- lme(
  Survey5.RESP ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_predictability,
  na.action = na.exclude
)

summary(m_pred_survey5)
anova(m_pred_survey5, type = "marginal")
emmeans(m_pred_survey5, ~ Condition * group)
pairs(emmeans(m_pred_survey5, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_pred_survey5, ~ group | Condition), adjust = "holm")


############################################################
# Step 3.2: Bias models
############################################################

respTask_bias <- dataset_respTask %>%
  filter(
    Condition %in% c("Resp2_Hea", "Resp2_Pat"),
    group %in% c("Controls", "Patients")
  ) %>%
  mutate(
    Condition = droplevels(Condition),
    Condition = relevel(Condition, ref = "Resp2_Hea"),
    group     = relevel(group, ref = "Controls")
  )

# PLV
m_bias_plv <- lme(
  PLV ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_bias,
  na.action = na.exclude
)

summary(m_bias_plv)
anova(m_bias_plv, type = "marginal")
emmeans(m_bias_plv, ~ Condition * group)
pairs(emmeans(m_bias_plv, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_bias_plv, ~ group | Condition), adjust = "holm")

# Relative phase
m_bias_phase <- lme(
  mean_relative_phase_deg ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_bias,
  na.action = na.exclude
)

summary(m_bias_phase)
anova(m_bias_phase, type = "marginal")
emmeans(m_bias_phase, ~ Condition * group)
pairs(emmeans(m_bias_phase, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_bias_phase, ~ group | Condition), adjust = "holm")

# WCC maximum correlation
m_bias_wcc_corr <- lme(
  WCC_mean_maxCorr ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_bias,
  na.action = na.exclude
)

summary(m_bias_wcc_corr)
anova(m_bias_wcc_corr, type = "marginal")
emmeans(m_bias_wcc_corr, ~ Condition * group)
pairs(emmeans(m_bias_wcc_corr, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_bias_wcc_corr, ~ group | Condition), adjust = "holm")

# WCC lag
m_bias_wcc_lag <- lme(
  WCC_mean_maxLag_sec ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_bias,
  na.action = na.exclude
)

summary(m_bias_wcc_lag)
anova(m_bias_wcc_lag, type = "marginal")
emmeans(m_bias_wcc_lag, ~ Condition * group)
pairs(emmeans(m_bias_wcc_lag, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_bias_wcc_lag, ~ group | Condition), adjust = "holm")

# Survey 1
m_bias_survey1 <- lme(
  Survey1.RESP ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_bias,
  na.action = na.exclude
)

summary(m_bias_survey1)
anova(m_bias_survey1, type = "marginal")
emmeans(m_bias_survey1, ~ Condition * group)
pairs(emmeans(m_bias_survey1, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_bias_survey1, ~ group | Condition), adjust = "holm")

# Survey 2
m_bias_survey2 <- lme(
  Survey2.RESP ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_bias,
  na.action = na.exclude
)

summary(m_bias_survey2)
anova(m_bias_survey2, type = "marginal")
emmeans(m_bias_survey2, ~ Condition * group)
pairs(emmeans(m_bias_survey2, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_bias_survey2, ~ group | Condition), adjust = "holm")

# Survey 3
m_bias_survey3 <- lme(
  Survey3.RESP ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_bias,
  na.action = na.exclude
)

summary(m_bias_survey3)
anova(m_bias_survey3, type = "marginal")
emmeans(m_bias_survey3, ~ Condition * group)
pairs(emmeans(m_bias_survey3, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_bias_survey3, ~ group | Condition), adjust = "holm")

# Survey 4
m_bias_survey4 <- lme(
  Survey4.RESP ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_bias,
  na.action = na.exclude
)

summary(m_bias_survey4)
anova(m_bias_survey4, type = "marginal")
emmeans(m_bias_survey4, ~ Condition * group)
pairs(emmeans(m_bias_survey4, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_bias_survey4, ~ group | Condition), adjust = "holm")

# Survey 5
m_bias_survey5 <- lme(
  Survey5.RESP ~ Condition * group,
  random = ~ 1 | ID,
  data = respTask_bias,
  na.action = na.exclude
)

summary(m_bias_survey5)
anova(m_bias_survey5, type = "marginal")
emmeans(m_bias_survey5, ~ Condition * group)
pairs(emmeans(m_bias_survey5, ~ Condition | group), adjust = "holm")
pairs(emmeans(m_bias_survey5, ~ group | Condition), adjust = "holm")


