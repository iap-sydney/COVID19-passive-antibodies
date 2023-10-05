# -------------------------------------------------------------------------
#'
#' This file is used to set up everything that's needed across the project
#' It loads libraries, creates functions, sets themes and defaults
#' 
# -------------------------------------------------------------------------


# load packages -----------------------------------------------------------

library(readxl)
library(writexl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(mvtnorm)
library(lme4)
library(patchwork)
library(readr)
library(metafor)
library(lmec)
library(mcr)
library(magic)
library(censReg)
library(lemon)
library(MethComp)
library(ggpubr)
library(data.table)


# set defaults ------------------------------------------------------------

# assume a total plasma volume of 3 liters or a concentration of 45ml/kg:
# (used to compute the fold-convalescent dose)
plasma_L <- 3
plasma_conc <- 45*1e-3 # concentration in l/kg

# other defaults:
conf.level <- 0.95 # confidence level
sd <- 0.4647092 # standard deviation of titers from Khoury et al. (Nature Medicine)
