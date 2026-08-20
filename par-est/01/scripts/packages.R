# Load packages

library(ggplot2)
library(data.table)
library(foreach)
library(doParallel)
library(doRNG)
devtools::load_all()

sink('../logs/R_versions.txt')
  print(sessionInfo())
sink()
