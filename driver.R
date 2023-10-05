# -------------------------------------------------------------------------
#' 
#' COVID-19 passive Antibody treatment analysis
#' 
#' Provides code to reproduce analysis in Stadler et al. (2023) "Determinants
#' of passive antibody efficacy in SARS-CoV-2 infection: a systematic review
#' and meta-analysis" published in The Lancet Microbe
#' 
# -------------------------------------------------------------------------


# setup -------------------------------------------------------------------

source("setup.R")


# processing --------------------------------------------------------------

# load data and clean data
source("processing/01_Data-preparation.R")


# analysis ----------------------------------------------------------------

source("analysis/01_Data-analysis.R")
source("analysis/02_Visualization-of-efficacy-and-stage.R")
source("analysis/03_Dose-response-curve.R")
source("analysis/04_Dose-response-curve-without-high-RoB.R")
source("analysis/05_Effect-of-treatment-stage.R")
source("analysis/06_Funnel-plots.R")
source("analysis/07_Predicting-efficacy.R")
