library(piecewiseSEM)

load("data/metrics.rdata")

model_1 <- piecewiseSEM::psem(
  lm(benthic_invertebrate_richness ~ fish_richness, data_average),
  lm(insect_richness ~ benthic_invertebrate_richness, data_average),
  lm(zooplankton_richness ~ insect_richness, data_average),
  lm(fungi_richness ~ benthic_invertebrate_richness + zooplankton_richness, data_average),
  lm(phytoplankton_synchrony ~ fish_synchrony + benthic_invertebrate_synchrony, data_average),
  lm(zooplankton_synchrony ~ fish_richness + zooplankton_richness + benthic_invertebrate_richness, data_average),
  lm(fungi_synchrony ~ fungi_richness + benthic_invertebrate_synchrony, data_average),
  lm(bacteria_synchrony ~ benthic_invertebrate_synchrony, data_average),
  lm(multigroups_synchrony ~ fish_synchrony + benthic_invertebrate_synchrony + bacteria_synchrony + fish_richness, data_average),
  lm(fish_stability ~ fish_synchrony, data_average),
  lm(benthic_invertebrate_stability ~ benthic_invertebrate_synchrony, data_average),
  lm(insect_stability ~ insect_synchrony, data_average),
  lm(zooplankton_stability ~ zooplankton_synchrony, data_average),
  lm(phytoplankton_stability ~ phytoplankton_synchrony, data_average),
  lm(bacteria_stability ~ bacteria_synchrony, data_average),
  lm(fungi_stability ~ fungi_synchrony, data_average),
  lm(multigroups_stability ~ multigroups_synchrony, data_average)
)

summary(model_1)

model_2 <- piecewiseSEM::psem(
  lm(benthic_invertebrate_diversity ~ fish_diversity, data_average),
  lm(insect_diversity ~ fish_diversity + benthic_invertebrate_diversity, data_average),
  lm(zooplankton_diversity ~ insect_diversity + fish_diversity, data_average),
  lm(fungi_diversity ~ benthic_invertebrate_diversity, data_average),
  lm(insect_synchrony ~ fish_diversity, data_average),
  lm(phytoplankton_synchrony ~ fish_synchrony + benthic_invertebrate_synchrony, data_average),
  lm(zooplankton_synchrony ~ fish_diversity + zooplankton_diversity + insect_diversity + insect_synchrony + benthic_invertebrate_diversity + benthic_invertebrate_synchrony, data_average),
  lm(fungi_synchrony ~ zooplankton_diversity + benthic_invertebrate_synchrony + insect_synchrony, data_average),
  lm(bacteria_synchrony ~ fish_diversity + benthic_invertebrate_synchrony + insect_diversity + insect_synchrony, data_average),
  lm(multigroups_synchrony ~ fish_synchrony + benthic_invertebrate_synchrony + bacteria_synchrony + fish_diversity, data_average),
  lm(fish_stability ~ fish_synchrony, data_average),
  lm(benthic_invertebrate_stability ~ benthic_invertebrate_synchrony, data_average),
  lm(insect_stability ~ insect_synchrony, data_average),
  lm(zooplankton_stability ~ zooplankton_synchrony, data_average),
  lm(phytoplankton_stability ~ phytoplankton_synchrony, data_average),
  lm(bacteria_stability ~ bacteria_synchrony, data_average),
  lm(fungi_stability ~ fungi_synchrony, data_average),
  lm(multigroups_stability ~ multigroups_synchrony, data_average)
)

summary(model_2)

model_3 <- psem(
  lm(mean_path_length ~ fish_richness + zooplankton_richness + benthic_invertebrate_richness, data_original),
  lm(connectance ~ mean_path_length + fish_richness + benthic_invertebrate_richness, data_original),
  lm(nestedness ~ connectance + fish_richness + benthic_invertebrate_richness, data_original),
  lm(modularity ~ mean_path_length + connectance + fish_richness + nestedness, data_original),
  lm(robustness ~ mean_path_length + connectance + fish_richness + zooplankton_richness + benthic_invertebrate_richness, data_original),
  lm(vulnerability ~ nestedness + fish_richness + benthic_invertebrate_richness + zooplankton_richness, data_original)
)

summary(model_3)

model_4 <- psem(
  lm(mean_path_length ~ fish_diversity + zooplankton_diversity + benthic_invertebrate_diversity, data_original),
  lm(connectance ~ mean_path_length, data_original),
  lm(nestedness ~ mean_path_length + connectance + fish_diversity + zooplankton_diversity + benthic_invertebrate_diversity, data_original),
  lm(modularity ~ mean_path_length + connectance + nestedness, data_original),
  lm(robustness ~ mean_path_length + connectance + nestedness + fish_diversity, data_original),
  lm(vulnerability ~ mean_path_length + connectance + nestedness + modularity + robustness + fish_diversity + zooplankton_diversity + benthic_invertebrate_diversity, data_original)
)

summary(model_4)