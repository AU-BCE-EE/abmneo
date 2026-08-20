# Load packages

library(ggplot2)
library(data.table)
devtools::load_all()

sink('../logs/R_versions.txt')
  print(sessionInfo())
sink()
