load("data/metrics.rdata")

run_lmer_models <- function(response_vars, data) {
  model_results <- list()
  summary_results <- list()
  vif_results <- list()
  r2_results <- list()
  
  for (resp in response_vars) {
    cat('\n=====================================\n')
    cat(paste('## environment -', resp), '\n')
    cat('=====================================\n')
    
    formula_str <- paste(
      resp, '~ T + pH + DO + COD + NH3N + TP + TN + EC + TUB + (1|year)'
    )
    model <- lmer(as.formula(formula_str), data = data)
    
    cat('\n--- model summary ---\n')
    model_summary <- summary(model)  # 现在会包含 p 值列
    print(model_summary)
    
    fixed_effects <- as.data.frame(model_summary$coefficients)
    fixed_effects$Variable <- rownames(fixed_effects)
    fixed_effects$Response <- resp
    rownames(fixed_effects) <- NULL
    
    random_effects <- as.data.frame(VarCorr(model))
    random_effects$Response <- resp
    
    cat('\n--- vif (collinearity) ---\n')
    vif_values <- car::vif(model)
    print(vif_values)
    
    vif_df <- data.frame(
      Variable = names(vif_values),
      VIF = as.numeric(vif_values),
      Response = resp,
      row.names = NULL
    )
    
    cat('\n--- r² value ---\n')
    r2_values <- performance::r2(model)
    print(r2_values)
    
    r2_df <- data.frame(
      R2_Type = names(r2_values),
      R2_Value = as.numeric(r2_values),
      Response = resp,
      row.names = NULL
    )
    
    model_results[[resp]] <- model
    summary_results[[resp]] <- list(
      fixed_effects = fixed_effects,
      random_effects = random_effects
    )
    vif_results[[resp]] <- vif_df
    r2_results[[resp]] <- r2_df
    
    cat('\n-------------------------------------\n\n')
  }
  
  all_results <- list(
    models = model_results,
    summary_stats = summary_results,
    vif_stats = do.call(rbind, vif_results),
    r2_stats = do.call(rbind, r2_results),
    all_fixed_effects = do.call(rbind, lapply(summary_results, function(x) x$fixed_effects)),
    all_random_effects = do.call(rbind, lapply(summary_results, function(x) x$random_effects))
  )
  
  return(all_results)
}

response_variables <- c(
  'phytoplankton_richness', 'phytoplankton_diversity',
  'bacteria_richness', 'bacteria_diversity',
  'fungi_richness', 'fungi_diversity',
  'zooplankton_richness', 'zooplankton_diversity',
  'insect_richness', 'insect_diversity',
  'benthic_invertebrate_richness', 'benthic_invertebrate_diversity',
  'fish_richness', 'fish_diversity',
  'multigroups_richness', 'multigroups_diversity',
  'mean_path_length', 'connectance',
  'modularity', 'nestedness',
  'robustness', 'vulnerability'
)

all_models <- run_lmer_models(response_vars = response_variables, data = data_original)
