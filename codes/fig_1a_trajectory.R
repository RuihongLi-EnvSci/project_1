library(vegan)
library(ggplot2)
library(dplyr)
library(ecotraj)
library(labdsv)

load("data/abundance_by_group.rdata")

multigroups <- rbind(
  fish, benthic_invertebrate, insect, zooplankton,
  fungi, bacteria, phytoplankton
)

taxa_list <- c(
  'phytoplankton', 'bacteria', 'fungi', 'zooplankton',
  'insect', 'benthic_invertebrate', 'fish', 'multigroups'
)

output_names <- c(
  'pcoa_phytoplankton.png', 'pcoa_bacteria.png', 'pcoa_fungi.png', 'pcoa_zooplankton.png',
  'pcoa_insect.png', 'pcoa_benthic_invertebrate.png', 'pcoa_fish.png', 'pcoa_multigroups.png'
)

year_colors <- c(
  '1' = '#647ADD', '2' = '#C0C0C0', '3' = '#FFA500', '4' = '#F5664D'
)

group_labels <- c('1st', '2nd', '3rd', '4th')

groups <- data.frame(
  sampling_site_order = paste0(
    rep(c("BB","BD","BN","CS","CT","FD","FJ","FL","JTW","LH","WS","WZ","XZY","YJC","YW","YY","ZG","ZX"), 4),
    "_", rep(1:4, each = 18)
  ),
  year_groups = rep(2022:2025, each = 18)
)

entities <- sub('_.*', '', groups$sampling_site_order)
surveys <- as.numeric(rep(1:4, times = 18))

process_taxon <- function(taxon_name, output_name) {
  data <- get(taxon_name)[, -(1:8)]
  sample_names <- colnames(data[, -1])
  site <- sub("_.*", "", sample_names)
  year <- sub(".*_", "", sample_names)
  
  ord <- order(site, as.numeric(year))
  data <- data[, c(1, ord + 1)]
  rownames(data) <- data[,1]
  
  data <- hellinger(as.data.frame(t(data[, -1])))
  d <- as.matrix(vegdist(data, method = 'bray'))
  
  traj <- defineTrajectories(d, entities, surveys)
  trajectory_pcoa <- trajectoryPCoA(traj)
  
  coords <- as.data.frame(trajectory_pcoa$points[,1:2])
  coords$site <- entities
  coords$group <- surveys
  colnames(coords) <- c('dim1','dim2','site','group')
  
  pc_eig <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)
  pc1 <- pc_eig[1]
  pc2 <- pc_eig[2]
  
  group_var <- factor(surveys, levels = 1:4, labels = group_labels)
  permanova <- adonis2(d ~ group_var, permutations = 999)
  
  r2_value <- sprintf("%.3f", permanova$R2[1])
  p_value <- permanova$`Pr(>F)`[1]
  p_label <- ifelse(is.na(p_value), "p = NA",
                    ifelse(p_value < 0.001, " P < 0.001",
                           ifelse(p_value < 0.01, " P < 0.01",
                                  ifelse(p_value < 0.05, " P < 0.05",
                                         paste0(" P = ", sprintf("%.3f", p_value))))))
  
  hulls <- do.call(rbind, lapply(split(coords, coords$group), function(df) df[chull(df$dim1, df$dim2), ]))
  
  centroids <- coords %>%
    group_by(group) %>%
    summarise(dim1 = mean(dim1), dim2 = mean(dim2)) %>%
    arrange(factor(group, levels = 1:4)) %>%
    mutate(label = group_labels)
  
  traj_segments <- data.frame(
    x = centroids$dim1[-nrow(centroids)],
    y = centroids$dim2[-nrow(centroids)],
    xend = centroids$dim1[-1],
    yend = centroids$dim2[-1]
  )
  
  p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
    geom_polygon(data = hulls, aes(fill = as.factor(group), group = group), alpha = 0.2, linewidth = 1.5, color = NA) +
    geom_point(size = 8, shape = 16, alpha = 0.55) +
    geom_point(data = centroids, aes(x = dim1, y = dim2), color = 'black', fill = year_colors[as.character(centroids$group)], shape = 2, size = 10, stroke = 2, alpha = 0.65) +
    geom_segment(data = traj_segments, aes(x = x, y = y, xend = xend, yend = yend), color = 'black', linewidth = 1, arrow = arrow(length = unit(0.5, 'cm'), type = "closed")) +
    xlab(paste0('PC 1 (', pc1, '%)')) + ylab(paste0('PC 2 (', pc2, '%)')) +
    scale_color_manual(values = year_colors, labels = group_labels, name = 'Year') +
    scale_fill_manual(values = year_colors, labels = group_labels, name = 'Year') +
    theme_bw() +
    theme(panel.border = element_blank(), axis.line = element_line(linewidth = 1.225, color = "black"),
          panel.grid = element_blank(), axis.title = element_text(size = 45), axis.text = element_text(size = 32),
          legend.position = 'none', plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed')
  
  ggsave(paste0('figures/fig_1a_trajectory/', output_name), plot = p, width = 8, height = 8, dpi = 600)
  return(p)
}

plots_list <- lapply(seq_along(taxa_list), function(i) process_taxon(taxa_list[i], output_names[i]))
