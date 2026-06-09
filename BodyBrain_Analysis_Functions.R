#########################################################################################
##################      BRAIN BODY INTERACTION ANALYSIS FUNCTIONS      ##################      
#########################################################################################


############################################################
# 1. GENERIC FUNCTIONS
############################################################

############################################################
# LABELS AND COLORS
############################################################

performance_labels <- c(
  WCC_mean_maxCorr = "WCC peak correlation",
  WCC_mean_maxLag_sec = "WCC peak lag",
  PLV = "PLV",
  RelativePhase_MeanCircular_deg = "Relative respiratory phase",
  Survey1.RESP = "Survey 1",
  Survey2.RESP = "Survey 2",
  Survey3.RESP = "Survey 3",
  Survey4.RESP = "Survey 4",
  Survey5.RESP = "Survey 5",
  SurveyMean = "Survey Mean"
)

group_colors <- c(
  controls = "#0072B2",
  patients = "#D55E00"
)

group_labels <- c(
  controls = "Controls",
  patients = "Patients"
)

############################################################
# GENERIC CONDITION SETUP
############################################################

make_condition_setup <- function(cond_levels, cond_labels = NULL) {
  
  if (is.null(cond_labels)) {
    cond_labels <- cond_levels
  }
  
  stopifnot(length(cond_levels) == 2)
  stopifnot(length(cond_labels) == 2)
  
  names(cond_labels) <- cond_levels
  
  list(
    levels = cond_levels,
    labels = cond_labels,
    ref = cond_levels[1],
    cond1 = cond_levels[1],
    cond2 = cond_levels[2],
    contrast_label = paste0(cond_labels[2], " − ", cond_labels[1]),
    colors = setNames(c("#009E73", "#CC79A7"), cond_levels)
  )
}

############################################################
# PRETTY LABEL FUNCTIONS
############################################################

pretty_metric <- function(metricVar) {
  dplyr::case_when(
    metricVar == "te_lin" ~ "Linear TE",
    metricVar == "te_knn" ~ "Nonlinear TE",
    TRUE ~ as.character(metricVar)
  )
}

pretty_node <- function(x) {
  dplyr::case_when(
    x == "rr" ~ "RR",
    x == "resp" ~ "Respiration",
    x == "delta" ~ "Delta",
    x == "theta" ~ "Theta",
    x == "alpha" ~ "Alpha",
    x == "alpha1" ~ "Alpha1",
    x == "alpha2" ~ "Alpha2",
    x == "beta" ~ "Beta",
    x == "beta1" ~ "Beta1",
    x == "beta2" ~ "Beta2",
    x == "gamma" ~ "Gamma",
    TRUE ~ x
  )
}

pretty_interaction <- function(interaction) {
  parts <- strsplit(as.character(interaction), "_")[[1]]
  paste(pretty_node(parts[1]), "→", pretty_node(parts[2]))
}


make_bbb_plot_title <- function(plot_dat,
                                scope = NULL) {
  scope_txt <- if (!is.null(scope)) {
    paste0(" ", scope)
  } else {
    ""
  }
  paste0(
    pretty_interaction(unique(plot_dat$interaction)),
    scope_txt,
    " (",
    pretty_metric(unique(plot_dat$metricVar)),
    ")"
  )
}


#############################################################
# 2. FUNCTIONS FOR GLOBAL EFFECTS ANALYSIS           
############################################################

############################################################
# Prepare electrode/global TE dataset
############################################################

prepare_te_electrode_data <- function(te_file,
                                      data_respSync,
                                      cond_setup) {
  
  data <- readr::read_csv(te_file, show_col_types = FALSE)
  
  data <- data %>%
    dplyr::mutate(
      subjectID = ifelse(subjectID == "031", "31", subjectID),
      group       = factor(group),
      condition   = factor(condition, levels = cond_setup$levels),
      interaction = factor(interaction),
      metricVar   = factor(metricVar),
      electrode   = factor(electrode)
    ) %>%
    dplyr::filter(condition %in% cond_setup$levels) %>%
    dplyr::mutate(
      condition = relevel(condition, ref = cond_setup$ref),
      group = relevel(factor(group), ref = "controls")
    )
  
  data_resp_sub <- data_respSync %>%
    dplyr::filter(condition %in% cond_setup$levels) %>%
    dplyr::mutate(
      subjectID = factor(subjectID),
      condition = factor(condition, levels = cond_setup$levels),
      group = factor(group)
    )
  
  data <- dplyr::left_join(
    data,
    data_resp_sub,
    by = c("subjectID", "condition", "group")
  )
  
  data <- data %>%
    dplyr::mutate(
      gender = factor(gender),
      age = as.numeric(age),
      condition = droplevels(condition)
    )
  
  data
}

############################################################
# Add interaction direction type
############################################################
add_interaction_type <- function(data) {
  
  data %>%
    dplyr::mutate(
      interaction_type = dplyr::case_when(
        
        # BODY–BODY
        interaction == "resp_rr" ~ "resp_to_heart",
        interaction == "rr_resp" ~ "heart_to_resp",
        
        # BODY–BRAIN
        grepl("^resp_", interaction) ~ "resp_to_brain",
        grepl("_resp$", interaction) ~ "brain_to_resp",
        grepl("^rr_", interaction)   ~ "heart_to_brain",
        grepl("_rr$", interaction)   ~ "brain_to_heart",
        
        TRUE ~ "other"
      )
    )
}


############################################################
# GLOBAL TE-behaviour associations
############################################################

run_global_te_behaviour <- function(data,
                                        cond_setup,
                                        performance_vars_global) {
  
  data_global_all <- data %>%
    dplyr::filter(
      electrode == "GLOBAL",
      condition %in% cond_setup$levels
    ) %>%
    dplyr::mutate(
      TE = as.numeric(metricValue),
      subjectID = factor(subjectID),
      group = relevel(factor(group), ref = "controls"),
      condition = relevel(
        factor(condition, levels = cond_setup$levels),
        ref = cond_setup$ref
      ),
      interaction = factor(interaction),
      metricVar = factor(metricVar)
    ) %>%
    add_interaction_type() %>%
    droplevels()
  
  global_all_results <- data.frame()
  
  for (performance_var in performance_vars_global) {
    
    tmp_data <- data_global_all %>%
      dplyr::mutate(performance = as.numeric(.data[[performance_var]]))
    
    combos <- tmp_data %>%
      dplyr::distinct(interaction, metricVar, interaction_type)
    
    for (i in seq_len(nrow(combos))) {
      
      b  <- combos$interaction[i]
      mv <- combos$metricVar[i]
      it <- combos$interaction_type[i]
      
      dat_tmp <- tmp_data %>%
        dplyr::filter(interaction == b, metricVar == mv) %>%
        dplyr::mutate(
          TE_z = as.numeric(scale(TE)),
          performance_z = as.numeric(scale(performance))
        ) %>%
        droplevels()
      
      if (dplyr::n_distinct(dat_tmp$subjectID) < 5) next
      
      m_general <- tryCatch(
        nlme::lme(
          performance_z ~ TE_z + group + condition,
          random = ~1 | subjectID,
          data = dat_tmp,
          method = "REML",
          na.action = na.exclude
        ),
        error = function(e) NULL
      )
      
      m_group <- tryCatch(
        nlme::lme(
          performance_z ~ TE_z * group + condition,
          random = ~1 | subjectID,
          data = dat_tmp,
          method = "REML",
          na.action = na.exclude
        ),
        error = function(e) NULL
      )
      
      m_condition <- tryCatch(
        nlme::lme(
          performance_z ~ TE_z * condition + group,
          random = ~1 | subjectID,
          data = dat_tmp,
          method = "REML",
          na.action = na.exclude
        ),
        error = function(e) NULL
      )
      
      extract_lme_term <- function(model, target_pattern, question_label) {
        
        if (is.null(model)) return(NULL)
        
        tab <- summary(model)$tTable
        target_row <- grep(target_pattern, rownames(tab), value = TRUE)
        if (length(target_row) == 0) return(NULL)
        target_row <- target_row[1]
        
        data.frame(
          question = question_label,
          interaction = as.character(b),
          interaction_type = as.character(it),
          metricVar = as.character(mv),
          performance_var = performance_var,
          term = target_row,
          beta = tab[target_row, "Value"],
          se = tab[target_row, "Std.Error"],
          t = tab[target_row, "t-value"],
          p = tab[target_row, "p-value"],
          stringsAsFactors = FALSE
        )
      }
      
      global_all_results <- rbind(
        global_all_results,
        extract_lme_term(m_general, "^TE_z$", "global_general_TE_performance_association"),
        extract_lme_term(m_group, "TE_z:group", "global_group_interaction"),
        extract_lme_term(m_condition, "TE_z:condition", "global_condition_interaction")
      )
    }
  }
  
  global_all_results <- global_all_results %>%
    dplyr::group_by(question, performance_var, metricVar) %>%
    dplyr::mutate(p_fdr = p.adjust(p, method = "fdr")) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(question, performance_var, metricVar, p)
  
  global_all_results
}

############################################################
# GLOBAL LMMs
############################################################

run_global_lmms <- function(data,
                                    cond_setup,
                                    output_dir = NULL,
                                    prefix = "GLOBAL") {
  
  dat_g <- data %>%
    dplyr::filter(electrode == "GLOBAL") %>%
    dplyr::mutate(
      subjectID = factor(subjectID),
      condition = relevel(
        factor(condition, levels = cond_setup$levels),
        ref = cond_setup$ref
      ),
      group = relevel(factor(group), ref = "controls"),
      interaction = factor(interaction),
      metricVar = factor(metricVar),
      age = as.numeric(age),
      gender = factor(gender)
    ) %>%
    add_interaction_type() %>%
    droplevels()
  
  combos <- dat_g %>%
    dplyr::distinct(interaction, metricVar, interaction_type)
  
  results_global <- tibble::tibble()
  posthoc_global <- tibble::tibble()
  
  for (i in seq_len(nrow(combos))) {
    
    b  <- combos$interaction[i]
    mv <- combos$metricVar[i]
    it <- combos$interaction_type[i]
    
    df <- dat_g %>%
      dplyr::filter(interaction == b, metricVar == mv) %>%
      droplevels()
    
    if (
      dplyr::n_distinct(df$condition) < 2 ||
      dplyr::n_distinct(df$group) < 2 ||
      dplyr::n_distinct(df$subjectID) < 5
    ) {
      next
    }
    
    model <- tryCatch(
      nlme::lme(
        metricValue ~ condition * group,# + age + gender,
        random = ~1 | subjectID,
        data = df,
        method = "REML",
        na.action = na.exclude
      ),
      error = function(e) NULL
    )
    
    if (is.null(model)) next
    
    aov_tab <- anova(model, type = "marginal")
    
    target_terms <- c("condition", "group", "condition:group")
    target_terms <- target_terms[target_terms %in% rownames(aov_tab)]
    
    tmp_results <- tibble::tibble(
      electrode = "GLOBAL",
      interaction = as.character(b),
      interaction_type = as.character(it),
      metricVar = as.character(mv),
      effect = target_terms,
      numDF = aov_tab[target_terms, "numDF"],
      denDF = aov_tab[target_terms, "denDF"],
      F_value = aov_tab[target_terms, "F-value"],
      p = aov_tab[target_terms, "p-value"]
    )
    
    results_global <- dplyr::bind_rows(results_global, tmp_results)
    
    sig_terms_unc <- target_terms[
      aov_tab[target_terms, "p-value"] < 0.05
    ]
    
    if ("condition" %in% sig_terms_unc) {
      tmp_cond <- as.data.frame(
        emmeans::contrast(
          emmeans::emmeans(model, ~ condition),
          method = "pairwise",
          adjust = "holm"
        )
      ) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(
          electrode = "GLOBAL",
          interaction = as.character(b),
          interaction_type = as.character(it),
          metricVar = as.character(mv),
          effect = "condition"
        )
      
      posthoc_global <- dplyr::bind_rows(posthoc_global, tmp_cond)
    }
    
    if ("group" %in% sig_terms_unc) {
      tmp_group <- as.data.frame(
        emmeans::contrast(
          emmeans::emmeans(model, ~ group),
          method = "pairwise",
          adjust = "holm"
        )
      ) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(
          electrode = "GLOBAL",
          interaction = as.character(b),
          interaction_type = as.character(it),
          metricVar = as.character(mv),
          effect = "group"
        )
      
      posthoc_global <- dplyr::bind_rows(posthoc_global, tmp_group)
    }
    
    if ("condition:group" %in% sig_terms_unc) {
      tmp_inter <- as.data.frame(
        emmeans::contrast(
          emmeans::emmeans(model, ~ condition | group),
          method = "pairwise",
          adjust = "holm"
        )
      ) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(
          electrode = "GLOBAL",
          interaction = as.character(b),
          interaction_type = as.character(it),
          metricVar = as.character(mv),
          effect = "condition:group"
        )
      
      posthoc_global <- dplyr::bind_rows(posthoc_global, tmp_inter)
    }
  }
  
  results_global <- results_global %>%
  dplyr::group_by(effect, metricVar) %>%
  dplyr::mutate(p_fdr_by_type = p.adjust(p, method = "fdr")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(effect, metricVar, p)
    
  sig_fdr_by_type <- results_global %>%
    dplyr::filter(p_fdr_by_type < 0.05)
  
  sig_fdr_keys <- sig_fdr_by_type %>%
    dplyr::select(interaction, interaction_type, metricVar, effect)
  
  if (nrow(posthoc_global) > 0 && nrow(sig_fdr_keys) > 0) {
    
    posthoc_global_fdr <- posthoc_global %>%
      dplyr::semi_join(
        sig_fdr_keys,
        by = c("interaction", "interaction_type", "metricVar", "effect")
      ) %>%
      dplyr::arrange(effect, metricVar, interaction_type, interaction, p.value)
    
  } else {
    
    posthoc_global_fdr <- tibble::tibble(
      contrast = character(),
      estimate = numeric(),
      SE = numeric(),
      df = numeric(),
      t.ratio = numeric(),
      p.value = numeric(),
      electrode = character(),
      interaction = character(),
      interaction_type = character(),
      metricVar = character(),
      effect = character()
    )
  }
  
  sig_fdr_by_type_clean <- sig_fdr_by_type %>%
    dplyr::mutate(
      F_value = round(F_value, 3),
      p = round(p, 4),
      p_fdr_by_type = round(p_fdr_by_type, 4)
    )
  
  posthoc_global_fdr_clean <- posthoc_global_fdr %>%
    dplyr::mutate(
      estimate = round(estimate, 5),
      SE = round(SE, 5),
      df = round(df, 2),
      t.ratio = round(t.ratio, 3),
      p.value = round(p.value, 4)
    )
  
  if (!is.null(output_dir)) {
    write.csv2(
      results_global,
      file.path(output_dir, paste0(prefix, "_all_omnibus.csv")),
      row.names = FALSE
    )
    
    write.csv2(
      sig_fdr_by_type_clean,
      file.path(output_dir, paste0(prefix, "_sig_FDR.csv")),
      row.names = FALSE
    )
    
    write.csv2(
      posthoc_global_fdr_clean,
      file.path(output_dir, paste0(prefix, "_posthoc_FDR.csv")),
      row.names = FALSE
    )
  }
  
  list(
    all_results = results_global,
    sig_fdr = sig_fdr_by_type_clean,
    posthoc_fdr = posthoc_global_fdr_clean
  )
}



############################################################
# GLOBAL TE–BEHAVIOUR ASSOCIATION PLOT
############################################################

plot_global_te_behaviour <- function(data_global, interaction_pick, metric_pick,
                                     performance_var, cond_setup,
                                     effect_type = c("condition", "group", "general"),
                                     title_type = c("predicts", "association")) {
  
  effect_type <- match.arg(effect_type)
  title_type <- match.arg(title_type)
  
  plot_dat <- data_global %>%
    dplyr::filter(electrode == "GLOBAL",
                  interaction == interaction_pick,
                  metricVar == metric_pick) %>%
    dplyr::mutate(
      condition = relevel(factor(condition, levels = cond_setup$levels), ref = cond_setup$ref),
      group = relevel(factor(group), ref = "controls"),
      performance = as.numeric(.data[[performance_var]]),
      TE_z = as.numeric(scale(metricValue)),
      performance_z = as.numeric(scale(performance))
    ) %>%
    dplyr::filter(!is.na(TE_z), !is.na(performance_z)) %>%
    droplevels()
  
  title_text <- if (title_type == "predicts") {
    paste0(pretty_interaction(interaction_pick), " predicts ", performance_labels[[performance_var]])
  } else {
    paste0(pretty_interaction(interaction_pick), " association with ", performance_labels[[performance_var]])
  }
  
  if (effect_type == "condition") {
    model <- nlme::lme(performance_z ~ TE_z * condition + group,
                       random = ~1 | subjectID, data = plot_dat,
                       method = "REML", na.action = na.exclude)
    
    p <- ggplot2::ggplot(plot_dat, ggplot2::aes(TE_z, performance_z, color = condition)) +
      ggplot2::geom_point(ggplot2::aes(shape = group), size = 2.4, alpha = .75) +
      ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 1.1, alpha = .15) +
      ggplot2::scale_color_manual(values = cond_setup$colors, labels = cond_setup$labels) +
      ggplot2::scale_shape_discrete(labels = group_labels)
  }
  
  if (effect_type == "group") {
    model <- nlme::lme(performance_z ~ TE_z * group + condition,
                       random = ~1 | subjectID, data = plot_dat,
                       method = "REML", na.action = na.exclude)
    
    p <- ggplot2::ggplot(plot_dat, ggplot2::aes(TE_z, performance_z, color = group)) +
      ggplot2::geom_point(ggplot2::aes(shape = condition), size = 2.4, alpha = .75) +
      ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 1.1, alpha = .15) +
      ggplot2::scale_color_manual(values = group_colors, labels = group_labels) +
      ggplot2::scale_shape_discrete(labels = cond_setup$labels)
  }
  
  if (effect_type == "general") {
    model <- nlme::lme(performance_z ~ TE_z + group + condition,
                       random = ~1 | subjectID, data = plot_dat,
                       method = "REML", na.action = na.exclude)
    
    p <- ggplot2::ggplot(plot_dat, ggplot2::aes(TE_z, performance_z)) +
      ggplot2::geom_point(ggplot2::aes(color = group, shape = condition), size = 2.4, alpha = .75) +
      ggplot2::geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1.1, alpha = .15) +
      ggplot2::scale_color_manual(values = group_colors, labels = group_labels) +
      ggplot2::scale_shape_discrete(labels = cond_setup$labels)
  }
  
  p <- p +
    ggplot2::labs(
      title = title_text,
      x = "Global TE (z-score)",
      y = paste0(performance_labels[[performance_var]], " (z-score)")
    ) +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = .5, size = 16),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      legend.position = "right",
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.box.spacing = grid::unit(2, "pt"),
      legend.spacing.y = grid::unit(3, "pt")
    )
  
  print(summary(model))
  list(plot = p, model = model, data = plot_dat)
}



#############################################################
# 2. CLUSTER ANALYSIS AND PLOTS          
############################################################


############################################################
# Raincloud plot function
############################################################

plot_raincloud_2x2 <- function(data_cluster, clusterID = NULL, plot_title = NULL,
                               ylab = "TE", xlab = "Condition",
                               condition_labels = NULL, condition_order = NULL,
                               group_labels = NULL, group_order = NULL,
                               rain_side = NULL, alpha = .50,
                               jitter_width = .03, jitter_height_prop = .02,
                               seed = 42, add_subject_lines = TRUE,
                               palette = NULL, base_size = 15) {
  
  infer_condition_setup <- function(cond_levels) {
    cl <- as.character(cond_levels)
    ord <- if (!is.null(condition_order)) condition_order else if (all(c("hp", "lp") %in% cl)) c("hp", "lp") else if (all(c("hea", "pat") %in% cl)) c("hea", "pat") else cl
    labs <- if (!is.null(condition_labels)) condition_labels else if (all(c("hp", "lp") %in% cl)) c(hp = "High pred.", lp = "Low pred.") else if (all(c("hea", "pat") %in% cl)) c(hea = "Healthy bias", pat = "Patient bias") else setNames(ord, ord)
    list(order = ord, labels = labs)
  }
  
  df <- data_cluster %>%
    dplyr::select(dplyr::any_of(c("subjectID", "group", "condition", "TE", "clusterID",
                                  "interaction", "metricVar"))) %>%
    dplyr::filter(!is.na(TE)) %>%
    droplevels()
  
  if (is.null(clusterID) && "clusterID" %in% names(df)) {
    u <- unique(as.character(df$clusterID))
    if (length(u) == 1) clusterID <- u
  }
  
  if (is.null(plot_title)) {
    metric_now <- unique(as.character(df$metricVar))[1]
    interaction_now <- unique(as.character(df$interaction))[1]
    
    te_type <- if (metric_now == "te_lin") {
      "Linear"
    } else if (metric_now == "te_knn") {
      "Nonlinear"
    } else {
      metric_now
    }
    
    plot_title <- paste0(pretty_interaction(interaction_now), " Cluster (", te_type, ")")
  }
  
  cond_setup <- infer_condition_setup(levels(df$condition))
  df$condition <- factor(df$condition, levels = cond_setup$order)
  df$group <- if (!is.null(group_order)) factor(df$group, levels = group_order) else factor(df$group)
  
  if (!is.null(group_labels)) {
    df$group <- dplyr::recode(as.character(df$group), !!!as.list(group_labels))
    df$group <- factor(df$group)
  }
  
  if (is.null(palette)) {
    palette <- if (all(c("controls", "patients") %in% levels(df$group))) {
      c(controls = "#0072B2", patients = "#D55E00")
    } else {
      setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[seq_along(levels(df$group))],
               levels(df$group))
    }
  } else {
    if (is.null(names(palette))) names(palette) <- levels(df$group)
  }
  
  if (is.null(rain_side)) {
    rain_side <- if (length(levels(df$group)) == 2 && length(levels(df$condition)) == 2) "f2x2" else "l"
  }
  
  y_rng <- diff(range(df$TE, na.rm = TRUE))
  jitter_height <- if (is.finite(y_rng) && y_rng > 0) jitter_height_prop * y_rng else 0
  pos_jit <- ggplot2::position_jitter(width = jitter_width, height = jitter_height, seed = seed)
  
  ggplot2::ggplot(df, ggplot2::aes(condition, TE, fill = group, color = group)) +
    ggrain::geom_rain(
      id.long.var = if (add_subject_lines) "subjectID" else NULL,
      rain.side = rain_side,
      alpha = alpha,
      seed = seed,
      point.args.pos = list(position = pos_jit),
      line.args = list(linewidth = .35, alpha = .25),
      line.args.pos = list(position = pos_jit)
    ) +
    ggplot2::scale_x_discrete(labels = cond_setup$labels) +
    ggplot2::scale_color_manual(values = palette, labels = group_labels) +
    ggplot2::scale_fill_manual(values = palette, labels = group_labels) +
    ggplot2::labs(title = plot_title, x = xlab, y = ylab) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = .5, size = 16),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      legend.position = "right",
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank()
    )
}


# Function for plot title
make_plot_title <- function(data_cluster, cluster_pick = NULL,
                            include_cluster = TRUE,
                            include_number = FALSE) {
  
  metricVar <- unique(as.character(data_cluster$metricVar))
  metricVar <- metricVar[!is.na(metricVar)][1]
  
  interaction <- unique(as.character(data_cluster$interaction))
  interaction <- interaction[!is.na(interaction)][1]
  
  te_label <- dplyr::case_when(
    metricVar == "te_lin" ~ "Linear",
    metricVar == "te_knn" ~ "Nonlinear",
    TRUE ~ metricVar
  )
  
  cluster_txt <- ""
  
  if (include_cluster) {
    cluster_txt <- " Cluster"
    
    if (include_number) {
      if (is.null(cluster_pick) && "clusterID" %in% names(data_cluster)) {
        u <- unique(as.character(data_cluster$clusterID))
        if (length(u) == 1) cluster_pick <- u
      }
      
      clu_num <- suppressWarnings(
        as.integer(sub(".*_clu_(\\d+)$", "\\1", cluster_pick))
      )
      
      if (!is.na(clu_num)) {
        cluster_txt <- paste0(" Cluster ", clu_num)
      }
    }
  }
  
  paste0(
    pretty_interaction(interaction),
    cluster_txt,
    " (",
    te_label,
    ")"
  )
}

############################################################
# LMM ON SELECTED CLUSTER
############################################################

run_cluster_lmm_analysis <- function(data, cluster_pick) {
  
  data_cluster <- data %>%
    dplyr::filter(clusterID == cluster_pick) %>%
    dplyr::filter(!is.na(TE))
  
  ttl <- make_plot_title(data_cluster, cluster_pick)
  
  model <- nlme::lme(
    TE ~ condition * group,
    random = ~ 1 | subjectID,
    data = data_cluster,
    na.action = na.exclude
  )
  
  anova_table <- anova(model)
  
  emm <- emmeans::emmeans(model, ~ condition * group)
  
  post_condition_within_group <- emmeans::contrast(
    emmeans::emmeans(model, ~ condition | group),
    method = "pairwise",
    adjust = "holm"
  )
  
  post_group_within_condition <- emmeans::contrast(
    emmeans::emmeans(model, ~ group | condition),
    method = "pairwise",
    adjust = "holm"
  )
  
  print(ttl)
  print(summary(model))
  print(anova_table)
  print(post_condition_within_group)
  print(post_group_within_condition)
  
  list(
    clusterID = cluster_pick,
    title = ttl,
    data = data_cluster,
    model = model,
    anova = anova_table,
    emmeans = emm,
    post_condition_within_group = post_condition_within_group,
    post_group_within_condition = post_group_within_condition
  )
}

############################################################
# EXTRACT MODEL TERM STATISTICS
############################################################

extract_lm_term <- function(model, target_pattern, question_label,
                            clusterID, interaction, metricVar,
                            electrodes, performance_var) {
  
  if (is.null(model)) return(NULL)
  
  tab <- summary(model)$coefficients
  target_row <- grep(target_pattern, rownames(tab), value = TRUE)
  
  if (length(target_row) == 0) return(NULL)
  
  target_row <- target_row[1]
  
  data.frame(
    question = question_label,
    clusterID = clusterID,
    interaction = interaction,
    metricVar = metricVar,
    electrodes = electrodes,
    performance_var = performance_var,
    term = target_row,
    beta = tab[target_row, "Estimate"],
    se = tab[target_row, "Std. Error"],
    t = tab[target_row, "t value"],
    p = tab[target_row, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}


############################################################
# DELTA GROUP INTERACTION PLOT
############################################################

plot_delta_group_interaction <- function(data_sig, cluster_pick, performance_var, cond_setup,
                                         title_type = c("predicts", "association")) {
  
  title_type <- match.arg(title_type)
  
  plot_dat <- data_sig %>%
    dplyr::filter(clusterID == cluster_pick) %>%
    dplyr::mutate(
      condition = factor(condition, levels = cond_setup$levels),
      performance = as.numeric(.data[[performance_var]])
    ) %>%
    dplyr::select(subjectID, group, clusterID, interaction, metricVar,
                  electrodes, condition, TE, performance) %>%
    tidyr::pivot_wider(names_from = condition, values_from = c(TE, performance)) %>%
    dplyr::mutate(
      delta_TE = .data[[paste0("TE_", cond_setup$cond2)]] -
        .data[[paste0("TE_", cond_setup$cond1)]],
      delta_performance = .data[[paste0("performance_", cond_setup$cond2)]] -
        .data[[paste0("performance_", cond_setup$cond1)]],
      delta_TE_z = as.numeric(scale(delta_TE)),
      delta_performance_z = as.numeric(scale(delta_performance)),
      group = relevel(factor(group), ref = "controls")
    ) %>%
    droplevels()
  
  metric_now <- unique(as.character(plot_dat$metricVar))[1]
  
  te_type <- dplyr::case_when(
    metric_now == "te_lin" ~ "Linear",
    metric_now == "te_knn" ~ "Nonlinear",
    TRUE ~ metric_now
  )
  
  title_text <- if (title_type == "predicts") {
    paste0(
      "Δ ",
      pretty_interaction(unique(plot_dat$interaction)),
      " (", te_type, ") predicts Δ ",
      performance_labels[[performance_var]]
    )
  } else {
    paste0(
      "Δ ",
      pretty_interaction(unique(plot_dat$interaction)),
      " (", te_type, ") association with Δ ",
      performance_labels[[performance_var]]
    )
  }
  
  model <- lm(delta_performance_z ~ delta_TE_z * group, data = plot_dat)
  
  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(delta_TE_z, delta_performance_z, color = group)) +
    ggplot2::geom_point(size = 2.4, alpha = .75) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 1.1, alpha = .15) +
    ggplot2::scale_color_manual(values = group_colors, labels = group_labels) +
    ggplot2::labs(
      title = title_text,
      x = "Δ TE (z-score)",
      y = paste0("Δ ", performance_labels[[performance_var]], " (z-score)"),
      caption = paste0("Δ = ", cond_setup$contrast_label)
    ) +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = .5, size = 16),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12),
      plot.caption = ggplot2::element_text(size = 10, hjust = 1),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      legend.position = "right",
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.box.spacing = grid::unit(2, "pt"),
      legend.spacing.y = grid::unit(3, "pt")
    )
  
  print(summary(model))
  list(plot = p, model = model, data = plot_dat)
}