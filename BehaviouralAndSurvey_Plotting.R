######################################################
#                                                    #
#           BEHAVIOURAL AND SURVEY PLOTTING          #
#                                                    #
######################################################

# Written by Francesco Bubbico 
# Last updated: June 2026

library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrain)
library(patchwork)

# Load final dataset
dataset_respTask <- read_excel("C:/Users/onesh/OneDrive/Desktop/Open data/Data/Final datasets/RespSync_FinalDataset.xlsx")

# Create subsets
respTask_predictability <- dataset_respTask %>%
  filter(Condition %in% c("Resp1_HP", "Resp1_LP"))

respTask_bias <- dataset_respTask %>%
  filter(Condition %in% c("Resp2_Hea", "Resp2_Pat"))

# Output directory
plot_dir <- "C:/Users/onesh/OneDrive/Desktop/Open data/Plots/behavioural"
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# Plotting function
plot_behaviour_raincloud <- function(data,
                                     outcome,
                                     ylab,
                                     title,
                                     condition_labels,
                                     group_colors,
                                     group_labels,
                                     base_size = 15,
                                     show_legend = FALSE) {
  
  df <- data %>%
    dplyr::select(ID, Condition, group, value = all_of(outcome)) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::mutate(
      ID = factor(ID),
      Condition = factor(Condition),
      group = factor(group)
    )
  
  p <- ggplot(
    df,
    aes(
      x = Condition,
      y = value,
      fill = group,
      color = group
    )
  ) +
    
    ggrain::geom_rain(
      id.long.var = NULL,
      rain.side = "f2x2",
      alpha = 0.45,
      seed = 42,
      
      boxplot.args = list(
        width = 0.08,
        alpha = 0.40,
        outlier.shape = NA,
        linewidth = 0.55
      ),
      
      point.args = list(
        size = 1.15,
        alpha = 0.40
      )
    ) +
    
    scale_x_discrete(
      labels = condition_labels,
      expand = expansion(mult = c(0.14, 0.14))
    ) +
    
    scale_color_manual(
      values = group_colors,
      labels = group_labels,
      name = NULL
    ) +
    
    scale_fill_manual(
      values = group_colors,
      labels = group_labels,
      name = NULL
    ) +
    
    labs(
      title = title,
      x = NULL,
      y = ylab
    ) +
    
    theme_classic(base_size = base_size) +
    
    theme(
      panel.grid = element_blank(),
      
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = base_size
      ),
      
      axis.title.y = element_text(
        face = "bold",
        size = base_size
      ),
      
      axis.text.x = element_text(
        color = "black",
        size = base_size - 1
      ),
      
      axis.text.y = element_text(
        color = "black",
        size = base_size - 1
      ),
      
      axis.line = element_line(linewidth = 0.7),
      axis.ticks = element_line(linewidth = 0.7),
      
      legend.position = ifelse(show_legend, "bottom", "none"),
      legend.justification = "center",
      legend.box.just = "center",
      legend.text = element_text(size = base_size - 1),
      
      plot.margin = margin(8, 12, 8, 12)
    )
  
  return(p)
}


############################################################
# PREDICTABILITY BEHAVIOURAL FIGURE
############################################################

condition_labels_pred <- c(
  "Resp1_HP" = "High\npred.",
  "Resp1_LP" = "Low\npred."
)

group_colors <- c(
  "Controls" = "#0072B2",
  "Patients" = "#D55E00"
)

group_labels <- c(
  "Controls" = "Controls",
  "Patients" = "Patients"
)

# Four panel
p_relphase <- plot_behaviour_raincloud(
  data = respTask_predictability,
  outcome = "mean_relative_phase_deg",
  ylab = "Relative phase (degrees)",
  title = "Relative phase",
  condition_labels = condition_labels_pred,
  group_colors = group_colors,
  group_labels = group_labels
)

p_plv <- plot_behaviour_raincloud(
  data = respTask_predictability,
  outcome = "PLV",
  ylab = "PLV",
  title = "PLV",
  condition_labels = condition_labels_pred,
  group_colors = group_colors,
  group_labels = group_labels
)

p_wcc_lag <- plot_behaviour_raincloud(
  data = respTask_predictability,
  outcome = "WCC_mean_maxLag_sec",
  ylab = "WCC peak lag (s)",
  title = "WCC peak lag",
  condition_labels = condition_labels_pred,
  group_colors = group_colors,
  group_labels = group_labels,
  show_legend = TRUE
)

p_wcc_corr <- plot_behaviour_raincloud(
  data = respTask_predictability,
  outcome = "WCC_mean_maxCorr",
  ylab = "WCC peak correlation",
  title = "WCC peak correlation",
  condition_labels = condition_labels_pred,
  group_colors = group_colors,
  group_labels = group_labels
)



fig_behaviour_predictability <- 
  ( p_plv | p_relphase) /
  (p_wcc_corr | p_wcc_lag) +
  patchwork::plot_annotation(tag_levels = "A")

fig_behaviour_predictability


ggsave(
  file.path(plot_dir, "Behaviour_Predictability_Raincloud.png"),
  fig_behaviour_predictability,
  width = 11,
  height = 8,
  dpi = 600,
  bg = "white"
)


############################################################
# BIAS BEHAVIOURAL FIGURE
############################################################

condition_labels_bias <- c(
  "Resp2_Hea" = "Healthy\nbias",
  "Resp2_Pat" = "Patient\nbias"
)

group_colors <- c(
  "Controls" = "#0072B2",
  "Patients" = "#D55E00"
)

group_labels <- c(
  "Controls" = "Controls",
  "Patients" = "Patients"
)

# Four panels
p_bias_relphase <- plot_behaviour_raincloud(
  data = respTask_bias,
  outcome = "mean_relative_phase_deg",
  ylab = "Relative phase (degrees)",
  title = "Relative phase",
  condition_labels = condition_labels_bias,
  group_colors = group_colors,
  group_labels = group_labels
)

p_bias_plv <- plot_behaviour_raincloud(
  data = respTask_bias,
  outcome = "PLV",
  ylab = "PLV",
  title = "PLV",
  condition_labels = condition_labels_bias,
  group_colors = group_colors,
  group_labels = group_labels
)

p_bias_wcc_corr <- plot_behaviour_raincloud(
  data = respTask_bias,
  outcome = "WCC_mean_maxCorr",
  ylab = "WCC peak correlation",
  title = "WCC peak correlation",
  condition_labels = condition_labels_bias,
  group_colors = group_colors,
  group_labels = group_labels
)

p_bias_wcc_lag <- plot_behaviour_raincloud(
  data = respTask_bias,
  outcome = "WCC_mean_maxLag_sec",
  ylab = "WCC peak lag (s)",
  title = "WCC peak lag",
  condition_labels = condition_labels_bias,
  group_colors = group_colors,
  group_labels = group_labels,
  show_legend = TRUE
)


fig_behaviour_bias <- 
  (p_bias_plv | p_bias_relphase) /
  (p_bias_wcc_corr | p_bias_wcc_lag) +
  patchwork::plot_annotation(tag_levels = "A")

fig_behaviour_bias


ggsave(
  file.path(plot_dir, "Behaviour_Bias_Raincloud.png"),
  fig_behaviour_bias,
  width = 11,
  height = 8,
  dpi = 600,
  bg = "white"
)



############################################################
# SURVEY PLOTS
############################################################

plot_survey_bars <- function(data,
                             conditions,
                             condition_labels,
                             plot_title,
                             out_file,
                             plot_dir,
                             base_size = 14) {
  
  survey_labels <- c(
    "Survey1.RESP" = "Survey 1",
    "Survey2.RESP" = "Survey 2",
    "Survey3.RESP" = "Survey 3",
    "Survey4.RESP" = "Survey 4",
    "Survey5.RESP" = "Survey 5"
  )
  
  
  fill_colors <- c(
    "Controls_Light" = "#9ECAE1",
    "Controls_Dark"  = "#0072B2",
    "Patients_Light" = "#FDBB84",
    "Patients_Dark"  = "#D55E00"
  )

  df_long <- data %>%
    pivot_longer(
      cols = matches("^Survey[1-5]\\.RESP$"),
      names_to = "Survey",
      values_to = "Response"
    ) %>%
    filter(Condition %in% conditions) %>%
    mutate(
      Condition = factor(Condition, levels = conditions),
      
      group = case_when(
        group %in% c("CONTROLLI", "Controls", "controls") ~ "Controls",
        group %in% c("PAZIENTI", "Patients", "patients") ~ "Patients",
        TRUE ~ as.character(group)
      ),
      
      group = factor(group, levels = c("Controls", "Patients")),
      
      Survey = factor(
        Survey,
        levels = names(survey_labels)
      ),
      
      fill_group = case_when(
        group == "Controls" & Condition == conditions[1] ~ "Controls_Dark",
        group == "Patients" & Condition == conditions[1] ~ "Patients_Dark",
        group == "Controls" & Condition == conditions[2] ~ "Controls_Light",
        group == "Patients" & Condition == conditions[2] ~ "Patients_Light"
      )
    )

  df_summary <- df_long %>%
    group_by(Survey, Condition, group, fill_group) %>%
    summarise(
      mean_response = mean(Response, na.rm = TRUE),
      se_response = sd(Response, na.rm = TRUE) /
        sqrt(sum(!is.na(Response))),
      .groups = "drop"
    )
  
  pd <- position_dodge(width = 0.70)
  
  p <- ggplot(
    df_summary,
    aes(
      x = Condition,
      y = mean_response,
      fill = fill_group
    )
  ) +
    
    geom_col(
      aes(group = group),
      width = 0.60,
      color = "black",
      linewidth = 0.25,
      position = pd
    ) +
    
    geom_errorbar(
      aes(
        ymin = mean_response - se_response,
        ymax = mean_response + se_response,
        group = group
      ),
      width = 0.16,
      linewidth = 0.55,
      position = pd
    ) +
    
    facet_wrap(
      ~ Survey,
      ncol = 5,
      labeller = labeller(
        Survey = as_labeller(survey_labels)
      )
    ) +
    
    scale_x_discrete(
      labels = condition_labels
    ) +
    
    scale_fill_manual(
      values = fill_colors,
      breaks = c(
        "Controls_Dark",
        "Patients_Dark"
      ),
      labels = c(
        "Controls",
        "Patients"
      ),
      name = NULL
    ) +
    
    scale_y_continuous(
      limits = c(0, 7),
      breaks = 1:7,
      expand = expansion(mult = c(0, 0.05))
    ) +
    
    labs(
      title = plot_title,
      x = NULL,
      y = "Mean response (1–7)"
    ) +
    
    theme_classic(base_size = base_size) +
    
    theme(
      panel.grid = element_blank(),
      
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      
      strip.background = element_rect(
        fill = "white",
        color = "black"
      ),
      
      strip.text = element_text(
        face = "bold"
      ),
      
      axis.text.x = element_text(
        color = "black",
        size = base_size - 2
      ),
      
      axis.text.y = element_text(
        color = "black"
      ),
      
      axis.title.y = element_text(
        face = "bold"
      ),
      
      legend.position = "bottom",
      legend.title = element_blank(),
      
      legend.text = element_text(
        size = base_size - 2
      )
    )
  
  ggsave(
    file.path(plot_dir, out_file),
    p,
    width = 14,
    height = 5,
    dpi = 600,
    bg = "white"
  )
  
  return(p)
}


# PREDICTABILITY SURVEY FIGURE

p_survey_predictability <- plot_survey_bars(
  data = dataset_respTask,
  conditions = c("Resp1_HP", "Resp1_LP"),
  condition_labels = c(
    "Resp1_HP" = "High\npred.",
    "Resp1_LP" = "Low\npred."
  ),
  plot_title = "Subjective ratings: predictability",
  out_file = "Survey_Predictability.png",
  plot_dir = plot_dir
)


# BIAS SURVEY FIGURE

p_survey_bias <- plot_survey_bars(
  data = dataset_respTask,
  conditions = c("Resp2_Hea", "Resp2_Pat"),
  condition_labels = c(
    "Resp2_Hea" = "Healthy\nbias",
    "Resp2_Pat" = "Patient\nbias"
  ),
  plot_title = "Subjective ratings: social bias",
  out_file = "Survey_Bias.png",
  plot_dir = plot_dir
)