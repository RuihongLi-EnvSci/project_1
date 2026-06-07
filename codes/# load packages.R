options(timeout = 1800)

pre.packages <- c(
  'openxlsx', 'tidyr', 'vegan', 'dplyr', 'codyn', 'tibble', 'igraph', 
  'bipartite', 'purrr', 'Hmisc', 'piecewiseSEM','tidyverse','rstatix',
  'car','lme4','lmerTest','multilevelTools','extraoperators','JWileymisc',
  'effectsize','influence.ME','GGally','MuMIn','sjstats','labdsv','betapart',
  'reshape2','ggtern','ggplot2', 'ecotraj','stringr', 
  'smacof','ggrepel','adespatial','agricolae','sf') 

installed <- installed.packages()[, "Package"] 

install.packages(pre.packages[!pre.packages %in% installed], dependencies = TRUE) 

lapply(pre.packages, library,character.only = TRUE)

if (!require(pairwiseAdonis)) {
  install.packages("devtools")
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}
library(pairwiseAdonis)
