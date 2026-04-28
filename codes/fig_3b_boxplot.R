library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(agricolae)

load("data/metrics.rdata")

data_processed <- data_original %>%
  separate(sampling_site_order, into = c("site", "year_code"), sep = "_", remove = FALSE) %>%
  mutate(
    year = case_when(
      year_code == "1" ~ "1st",
      year_code == "2" ~ "2nd", 
      year_code == "3" ~ "3rd",
      year_code == "4" ~ "4th",
      TRUE ~ year_code
    ),
    year = factor(year, levels = c("1st", "2nd", "3rd", "4th")),
    year_numeric = as.numeric(year)
  )

data_long <- data_processed %>%
  pivot_longer(
    cols = c(connectance, mean_path_length, modularity, nestedness, vulnerability, robustness),
    names_to = "index",
    values_to = "value"
  ) %>%
  mutate(
    index = factor(index, levels = c("connectance", "mean_path_length", "modularity", "nestedness", "vulnerability", "robustness"))
  )

year_colors <- c(
  "1st" = "#647ADD",
  "2nd" = "#C0C0C0",
  "3rd" = "#FFA500",
  "4th" = "#F5664D"
)

index_names <- c(
  "connectance" = "Connectance",
  "mean_path_length" = "Mean Path Length",
  "modularity" = "Modularity",
  "nestedness" = "Nestedness",
  "vulnerability" = "Vulnerability",
  "robustness" = "Robustness"
)

y_labels <- c(
  "connectance" = "Connectance",
  "mean_path_length" = "Mean Path Length",
  "modularity" = "Modularity",
  "nestedness" = "Nestedness",
  "vulnerability" = "Vulnerability",
  "robustness" = "Robustness"
)

y_params <- list(
  connectance = list(limits = c(0, 0.15)),
  mean_path_length = list(limits = c(1, 2)),
  modularity = list(limits = c(0, 0.18)),
  nestedness = list(limits = c(0, 6.5)),
  robustness = list(limits = c(0, 0.35)),
  vulnerability = list(limits = c(10, 75))
)

for (current_index in levels(data_long$index)) {
  index_data <- data_long %>%
    filter(index == current_index) %>%
    drop_na(value)
  
  index_data$year <- factor(index_data$year, levels = c("1st", "2nd", "3rd", "4th"))
  
  lm_model <- lm(value ~ year, data = index_data)
  anova_res <- anova(lm_model)
  p_value <- anova_res$`Pr(>F)`[1]
  f_value <- anova_res$`F value`[1]
  f_value_text <- sprintf("F = %.2f", f_value)
  
  p_label <- if (is.na(p_value)) {
    "p = NA"
  } else if (p_value < 0.001) {
    "P < 0.001"
  } else if (p_value < 0.01) {
    "P < 0.01"
  } else if (p_value < 0.05) {
    "P < 0.05"
  } else {
    paste0("P = ", sprintf("%.3f", p_value))
  }
  
  anova_text <- paste(f_value_text, '', p_label)
  
  anova_model <- aov(value ~ year, data = index_data)
  tukey_result <- HSD.test(anova_model, "year", group = TRUE, console = FALSE)
  tukey_groups <- tukey_result$groups
  tukey_groups <- tukey_groups[levels(index_data$year), , drop = FALSE]
  
  boxplot_stats <- index_data %>%
    group_by(year) %>%
    dplyr::summarize(
      Q1 = quantile(value, 0.25, na.rm = TRUE),
      Q3 = quantile(value, 0.75, na.rm = TRUE),
      Median = median(value, na.rm = TRUE),
      IQR = IQR(value, na.rm = TRUE),
      Upper_whisker = min(max(value, na.rm = TRUE), Q3 + 1.5 * IQR),
      Lower_whisker = max(min(value, na.rm = TRUE), Q1 - 1.5 * IQR),
      Max_data = max(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  letter_df <- data.frame(
    year = rownames(tukey_groups),
    letter = tukey_groups$groups,
    stringsAsFactors = FALSE
  ) %>%
    mutate(year = factor(year, levels = levels(index_data$year))) %>%
    left_join(boxplot_stats, by = "year")
  
  y_range <- diff(y_params[[current_index]]$limits)
  
  letter_df <- letter_df %>%
    mutate(
      base_y = pmax(Upper_whisker, Q3),
      letter_y = base_y + y_range * 0.085
    )
  
  if (current_index == "nestedness") {
    letter_df <- letter_df %>%
      mutate(letter_y = ifelse(letter_y > 1.5 & letter_y < 2.5, 2.55, letter_y))
  } else if (current_index == "connectance") {
    letter_df <- letter_df %>%
      mutate(letter_y = ifelse(letter_y > 0.3 & letter_y < 0.4, 0.405, letter_y))
  }
  
  y_limit_max <- max(letter_df$letter_y, na.rm = TRUE)
  final_y_limit_max <- max(y_params[[current_index]]$limits[2], y_limit_max) + y_range * 0.09
  y_limits <- c(y_params[[current_index]]$limits[1], final_y_limit_max)
  
  anova_x <- 0.55
  anova_y <- y_params[[current_index]]$limits[1] + y_range * 0.04
  
  p_base <- ggplot(index_data, aes(x = year, y = value, fill = year, color = year)) +
    geom_boxplot(
      width = 0.575,
      alpha = 0,
      color = "black",
      linewidth = 1.1,
      outlier.shape = NA
    ) +
    geom_point(
      position = position_jitter(width = 0.2, height = 0),
      size = 7,
      alpha = 0.7,
      aes(color = year)
    ) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.4,
      linewidth = 0.6,
      color = "black",
      alpha = 0.8
    ) +
    scale_fill_manual(values = year_colors) +
    scale_color_manual(values = year_colors) +
    scale_y_continuous(
      name = y_labels[current_index],
      labels = function(x) {
        ifelse(round(x, 0) == x, as.character(round(x)), sprintf("%.2f", x))
      },
      limits = y_limits,
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_discrete(name = "Year") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 2),
      axis.line.x = element_line(color = "black", linewidth = 1.25),
      axis.line.y = element_line(color = "black", linewidth = 1.25),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 50, color = "black"),
      axis.text.y = element_blank(),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right = element_blank(),
      axis.title.y.right = element_blank(),
      axis.ticks.length = unit(0.3, "cm"),
      legend.position = "none",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  
  output_file <- paste0(
    "figures/fig_3b_boxplot/",
    current_index,
    ".png"
  )
  
  ggsave(
    filename = output_file,
    plot = p_base,
    width = 8,
    height = 6.75,
    dpi = 600
  )
  
  cat("\nANOVA for", index_names[current_index], ":\n")
  cat(anova_text, "\n\n")
  cat("Tukey HSD Results for", index_names[current_index], ":\n")
  print(tukey_groups)
}