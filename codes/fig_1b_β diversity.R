library(ggtern)
library(reshape2)
library(ade4)
library(viridis)
library(ggplot2)

load("data/abundance_by_group.rdata")

datasets_to_analyze <- c(
  "fish", "benthic_invertebrate", "insect", "zooplankton",
  "fungi", "bacteria", "phytoplankton"
)

multigroups <- rbind(
  fish, benthic_invertebrate, insect, zooplankton,
  fungi, bacteria, phytoplankton
)

multigroups[multigroups > 0] <- 1
assign("multigroups", multigroups)

all_datasets <- c(datasets_to_analyze, "multigroups")

analyze_beta_diversity <- function(dataset_name, data) {
  data[data > 0] <- 1
  data <- t(as.matrix(data[, -(1:9)]))
  
  beta_com <- beta.div.comp(data, coef = "J", quant = FALSE, save.abc = FALSE)
  
  repl <- as.matrix(beta_com$repl)
  repl_df <- melt(repl)
  repl_df <- repl_df[repl_df$Var1 != repl_df$Var2, ]
  repl_df <- repl_df[as.numeric(repl_df$Var1) < as.numeric(repl_df$Var2), ]
  rownames(repl_df) <- NULL
  colnames(repl_df)[3] <- "Repl"
  
  rich <- as.matrix(beta_com$rich)
  rich_df <- melt(rich)
  rich_df <- rich_df[rich_df$Var1 != rich_df$Var2, ]
  rich_df <- rich_df[as.numeric(rich_df$Var1) < as.numeric(rich_df$Var2), ]
  rownames(rich_df) <- NULL
  colnames(rich_df)[3] <- "RichDiff"
  
  data_df <- merge(repl_df, rich_df)
  data_df$BDtotal <- data_df$Repl + data_df$RichDiff
  
  data_df <- data_df[
    substr(data_df$Var1, 1, regexpr("_", data_df$Var1) - 1) ==
      substr(data_df$Var2, 1, regexpr("_", data_df$Var2) - 1),
  ]
  rownames(data_df) <- NULL
  
  summary_stats <- data.frame(
    Dataset = dataset_name,
    `BDtotal (total diversity)` = mean(data_df$BDtotal, na.rm = TRUE),
    `Repl (turnover)` = mean(data_df$Repl, na.rm = TRUE),
    `RichDiff (nestedness)` = mean(data_df$RichDiff, na.rm = TRUE),
    `Repl/BDtotal` = mean(data_df$Repl / data_df$BDtotal, na.rm = TRUE),
    `RichDiff/BDtotal` = mean(data_df$RichDiff / data_df$BDtotal, na.rm = TRUE),
    `Number of sample pairs` = nrow(data_df),
    check.names = FALSE
  )
  
  p <- ggtern(
    data = data_df,
    aes(x = RichDiff, y = Repl, z = 1 - BDtotal)
  ) +
    stat_density_tern(
      geom = "polygon",
      aes(fill = ..level.., alpha = ..level..),
      base = "ilr",
      n = 100,
      bins = 20,
      alpha = 0.175
    ) +
    scale_fill_viridis_c(
      breaks = function(x) pretty(x, n = 3),
      guide = guide_colorbar(
        barwidth = 2,
        barheight = 7,
        label.size = 1,
        title.position = "left",
        title.hjust = 0.5
      )
    ) +
    geom_point(size = 3.875, color = "grey", fill = "grey", alpha = 0.95) +
    scale_alpha(range = c(0.2, 0.6), guide = "none") +
    theme_bw() +
    theme_showarrows() +
    labs(
      title = paste0(
        "Beta Diversity Decomposition - ",
        tools::toTitleCase(dataset_name),
        " (Jaccard)"
      ),
      x = "RichDiff",
      y = "Repl",
      z = "Similarity"
    ) +
    theme(
      tern.axis.arrow = element_line(size = 1.2, color = "black", alpha = 0.7),
      tern.axis.arrow.sep = 0.1,
      tern.axis.title.T = element_blank(),
      tern.axis.title.L = element_blank(),
      tern.axis.title.R = element_blank(),
      tern.axis.arrow.text.T = element_text(size = 33, vjust = 0.3),
      tern.axis.arrow.text.L = element_text(size = 33, vjust = 0.3),
      tern.axis.arrow.text.R = element_text(size = 33, vjust = 0.6),
      tern.axis.text = element_text(size = 23.5, color = "black", family = "Arial"),
      plot.title = element_blank(),
      legend.position = c(0.81, 0.81),
      legend.box.just = "right",
      legend.background = element_blank(),
      legend.text = element_text(size = 21.5),
      legend.title = element_text(size = 19),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
    )
  
  ggsave(
    paste0("figures/fig_1b_β diversity/jac_", dataset_name, ".png"),
    p,
    width = 9,
    height = 9,
    dpi = 100
  )
  
  return(summary_stats)
}

all_results <- list()

for (dataset_name in all_datasets) {
  data <- get(dataset_name)
  result <- analyze_beta_diversity(dataset_name, data)
  all_results[[length(all_results) + 1]] <- result
}

final_results <- do.call(rbind, all_results)

print(final_results)