load("data/metrics.rdata")

# ========================
# 1. Trajectory PCoA plot
# ========================

variables <- c(
  "bacteria_richness","benthic_invertebrate_richness",
  "fish_richness","fungi_richness",
  "insect_richness","phytoplankton_richness","zooplankton_richness",
  "bacteria_diversity","benthic_invertebrate_diversity",
  "fish_diversity","fungi_diversity",
  "insect_diversity","phytoplankton_diversity",
  "zooplankton_diversity","bacteria_PC1",
  "benthic_invertebrate_PC1","fish_PC1",
  "fungi_PC1","insect_PC1",
  "phytoplankton_PC1","zooplankton_PC1",
  "connectance","mean_path_length",
  "modularity","nestedness",
  "vulnerability","robustness"
)

data <- data_original[, variables]
data <- as.data.frame(t(data))
data <- cbind(rownames(data),data)
  
colnames(data) <- c(
  "BB_1","BB_2","BB_3","BB_4","BD_1","BD_2","BD_3","BD_4",
  "BN_1","BN_2","BN_3","BN_4","CS_1","CS_2","CS_3","CS_4",
  "CT_1","CT_2","CT_3","CT_4","FD_1","FD_2","FD_3","FD_4",
  "FJ_1","FJ_2","FJ_3","FJ_4","FL_1","FL_2","FL_3","FL_4",
  "JTW_1","JTW_2","JTW_3","JTW_4","LH_1","LH_2","LH_3","LH_4",
  "WS_1","WS_2","WS_3","WS_4","WZ_1","WZ_2","WZ_3","WZ_4",
  "XZY_1","XZY_2","XZY_3","XZY_4","YJC_1","YJC_2","YJC_3","YJC_4",
  "YW_1","YW_2","YW_3","YW_4","YY_1","YY_2","YY_3","YY_4",
  "ZG_1","ZG_2","ZG_3","ZG_4","ZX_1","ZX_2","ZX_3","ZX_4"
)

year_colors <- c(
  '1' = '#496C8D',
  '2' = '#D3D3D3',
  '3' = '#669189',
  '4' = '#D39139'
)
centre_colors <- c(
  '1' = '#496C8D',
  '2' = '#666666',
  '3' = '#669189',
  '4' = '#D39139'
)
group_labels <- c('1st', '2nd', '3rd', '4th')

groups <- data.frame(
  sampling_site_order = c(
    "BB_1","BD_1","BN_1","CS_1","CT_1","FD_1","FJ_1","FL_1","JTW_1","LH_1",
    "WS_1","WZ_1","XZY_1","YJC_1","YW_1","YY_1","ZG_1","ZX_1",
    "BB_2","BD_2","BN_2","CS_2","CT_2","FD_2","FJ_2","FL_2","JTW_2","LH_2",
    "WS_2","WZ_2","XZY_2","YJC_2","YW_2","YY_2","ZG_2","ZX_2",
    "BB_3","BD_3","BN_3","CS_3","CT_3","FD_3","FJ_3","FL_3","JTW_3","LH_3",
    "WS_3","WZ_3","XZY_3","YJC_3","YW_3","YY_3","ZG_3","ZX_3",
    "BB_4","BD_4","BN_4","CS_4","CT_4","FD_4","FJ_4","FL_4","JTW_4","LH_4",
    "WS_4","WZ_4","XZY_4","YJC_4","YW_4","YY_4","ZG_4","ZX_4"
  ),
  year_groups = c(
    rep(2022, 18),
    rep(2023, 18),
    rep(2024, 18),
    rep(2025, 18)
  )
)

entities <- sub('_.*', '', groups$sampling_site_order)
surveys  <- as.numeric(rep(c('1','2','3','4'), times = 18))

sample_names <- colnames(data[, -1])
site  <- sub("_.*", "", sample_names)
year  <- sub(".*_", "", sample_names)

ord <- order(site, as.numeric(year))
data <- data[, c(1, ord + 1)]

rownames(data) <- data[, 1]
data <- scale(t(data[, -1]))
d <- as.matrix(vegdist(data, method = 'euclidean'))

trajectory_pcoa <- trajectoryPCoA(defineTrajectories(d, entities, surveys))
print(trajectoryLengths(defineTrajectories(d, entities, surveys)))

trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
coords <- trajectory_pcoa_point[, 1:2]
coords$site <- entities
coords$group <- surveys

colnames(coords) <- c('dim1', 'dim2', 'site', 'group')

pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]

group_var <- factor(surveys, levels = c('1','2','3','4'), labels = c('1st','2nd','3rd','4th'))
permanova_result <- adonis2(d ~ group_var, permutations = 999)

r2_value <- sprintf("%.3f", permanova_result$R2[1])
p_value  <- permanova_result$`Pr(>F)`[1]

mycol <- year_colors

hulls <- do.call(rbind, lapply(split(coords, coords$group), function(df) {
  df[chull(df$dim1, df$dim2), ]
}))

centroids <- coords %>%
  group_by(group) %>%
  summarise(
    dim1 = mean(dim1),
    dim2 = mean(dim2)
  ) %>%
  arrange(factor(group, levels = c('1','2','3','4'))) %>%
  mutate(label = group_labels)

traj_segments <- data.frame(
  x     = centroids$dim1[-nrow(centroids)],
  y     = centroids$dim2[-nrow(centroids)],
  xend  = centroids$dim1[-1],
  yend  = centroids$dim2[-1]
)

p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
  geom_polygon(
    data = hulls,
    aes(fill = as.factor(group), group = group),
    alpha = 0.2,
    linewidth = 1.5,
    color = NA
  ) +
  geom_point(size = 8, shape = 16, alpha = 0.55) +
  geom_point(
    data = centroids,
    aes(x = dim1, y = dim2),
    shape  = 2,
    color = centre_colors,
    size   = 10,
    stroke = 2,
    alpha = 0.65
  ) +
  geom_segment(
    data = traj_segments, linewidth = 1,
    aes(x = x, y = y, xend = xend, yend = yend),
    color  = 'black',
    alpha = 0.6,
    arrow  = arrow(length = unit(0.4, 'cm'), type = "closed")
  ) +
  xlab(paste0('PC 1 (', pc1, '%)')) +
  ylab(paste0('PC 2 (', pc2, '%)')) +
  scale_color_manual(values = mycol, labels = group_labels, name = 'Year') +
  scale_fill_manual(values = mycol, labels = group_labels, name = 'Year') +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 1, color = "black"),
    panel.grid = element_blank(),
    axis.title = element_text(size = 37),
    axis.text = element_text(size = 27),
    legend.position = 'none',
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed')

ggsave(
  'figures/fig_5a_trajectory/metrics_pcoa.png',
  plot = p, width = 10, height = 8, dpi = 1000
)

# ========================
# 2. PCA loadings plot
# ========================

# Rebuild the transposed, raw variable matrix for PCA
data_raw <- as.data.frame(t(data_original[, variables]))
data_raw <- cbind(variable = rownames(data_raw), data_raw)

sample_names <- colnames(data[, -1])
site <- sub("_.*", "", sample_names)
year <- sub(".*_", "", sample_names)

ord <- order(site, as.numeric(year))
data_raw <- data_raw[, c(1, ord + 1)]
rownames(data_raw) <- data_raw[, 1]
X_scaled <- scale(t(data_raw[, -1]))
X_scaled <- X_scaled[complete.cases(X_scaled), ]

pca_res <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
loadings <- pca_res$rotation[, 1:2]
colnames(loadings) <- c("PC1", "PC2")

loading_df <- data.frame(
  variable = rownames(loadings),
  PC1 = loadings[, "PC1"],
  PC2 = loadings[, "PC2"]
)

loading_df <- loading_df %>%
  mutate(category = case_when(
    grepl("richness", variable)      ~ "Richness",
    grepl("diversity", variable)     ~ "Diversity",
    grepl("PC1", variable)           ~ "PC1",
    variable %in% c("connectance", "mean_path_length",
                    "modularity", "nestedness",
                    "vulnerability", "robustness") ~ "Foodweb",
    TRUE ~ "Other"
  ))

cat_colors <- c(
  "Richness"  = "#496C8D",
  "Diversity" = "#C0C0C0",
  "PC1"       = "#669189",
  "Foodweb"   = "#D39139"
)

var_explained <- round((pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100, 2)
pc1_lab <- "Loadings"
pc2_lab <- "Loadings"

p_load <- ggplot(loading_df, aes(x = PC1, y = PC2, color = category)) +
  geom_segment(aes(xend = 0, yend = 0),
               linewidth = 0.7, alpha = 0.8, show.legend = FALSE) +
  geom_point(aes(color = category),
             shape = 15, size = 6.95, show.legend = TRUE) +
  geom_text_repel(aes(label = gsub("(_richness|_diversity|_PC1)", "", variable)),
                  size = 7, show.legend = FALSE,
                  box.padding = 0.5, point.padding = 0.2,
                  force = 6, segment.color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = cat_colors, name = "Category") +
  labs(x = pc1_lab, y = pc2_lab) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 1, color = "black"),
    panel.grid = element_blank(),
    axis.title = element_text(size = 37),
    axis.text = element_text(size = 27),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.position = c(0.20, 0.9),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  guides(color = guide_legend(
    nrow = 2,
    override.aes = list(shape = 15, size = 4, linetype = "blank")
  ))

ggsave('figures/fig_5a_trajectory/pcoa_loadings.png',
       plot = p_load, width = 10, height = 8, dpi = 1000)

cat("PCA R2 value:\n")
cat("PC1:", var_explained[1], "%\n")
cat("PC2:", var_explained[2], "%\n")

# ===================================
# 3. Segment distance scatter plot
# ===================================

traj_length <- as.data.frame(
  trajectoryLengths(defineTrajectories(d, entities, surveys))[, 1:3]
)
colnames(traj_length) <- c('1st-2nd', '2nd-3rd', '3rd-4th')

entity_list <- unique(entities)
dist_1_4 <- numeric(length(entity_list))
dist_2_4 <- numeric(length(entity_list))

for (i in seq_along(entity_list)) {
  entity <- entity_list[i]
  idx <- which(entities == entity)
  idx <- idx[order(surveys[idx])]
  dist_1_4[i] <- d[idx[1], idx[4]]
  dist_2_4[i] <- d[idx[2], idx[4]]
}

traj_length$'2nd-4th' <- dist_2_4
traj_length$'1st-4th' <- dist_1_4

data_long <- pivot_longer(
  traj_length,
  cols = c('1st-2nd', '2nd-3rd', '3rd-4th', '2nd-4th', '1st-4th'),
  names_to = 'Segment',
  values_to = 'Length'
)

seg_levels <- c('1st-2nd', '2nd-3rd', '3rd-4th', '2nd-4th', '1st-4th')
data_long$Segment <- factor(data_long$Segment, levels = seg_levels)

mean_values <- data_long %>%
  group_by(Segment) %>%
  summarise(Mean = mean(Length), .groups = 'drop') %>%
  mutate(Segment = factor(Segment, levels = seg_levels))

data_long_1 <- data_long[data_long$Segment %in% c("1st-2nd", "2nd-3rd", "3rd-4th"), ]

anova_result <- aov(Length ~ Segment, data = data_long_1)
cat("===== ANOVA Results =====\n")
print(summary(anova_result))
cat("\n===== Tukey HSD =====\n")
print(TukeyHSD(anova_result))

f_value <- summary(anova_result)[[1]]$`F value`[1]
p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1]
annot_text <- paste0('F = ', round(f_value, 2),
                     '  P = ', ifelse(p_value < 0.001, '<0.001', round(p_value, 3)))

segment_colors <- c(
  '1st-2nd' = '#C0C0C0',
  '2nd-3rd' = '#496C8D',
  '3rd-4th' = '#D39139',
  '2nd-4th' = '#669189',
  '1st-4th' = '#996633'
)

p_scatter <- ggplot(data_long, aes(x = Segment, y = Length, fill = Segment)) +
  geom_jitter(
    position = position_jitter(width = 0.105, height = 0),
    color = 'black', stroke = 0, size = 10.5, alpha = 0.4, shape = 21
  ) +
  geom_point(
    data = mean_values,
    aes(x = Segment, y = Mean, fill = Segment),
    color = 'black', shape = 21, size = 13.5, stroke = 1.2
  ) +
  geom_text(
    data = mean_values,
    aes(x = Segment, y = Mean + 0.8, label = round(Mean, 2)),
    size = 9, fontface = 'bold', color = 'black'
  ) +
  xlab('Segment') +
  ylab('Distance') +
  scale_fill_manual(values = segment_colors) +
  scale_color_manual(values = segment_colors) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.text.y = element_text(size = 28, color = 'black'),
    axis.text.x = element_text(size = 27, color = 'black'),
    axis.title = element_text(size = 37),
    axis.ticks.length = unit(0.35, 'cm'),
    legend.position = 'none',
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm')
  )

ggsave(
  'figures/fig_5a_trajectory/traj_distance.png',
  plot = p_scatter, width = 10, height = 8, dpi = 1000
)
