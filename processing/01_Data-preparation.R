# -------------------------------------------------------------------------
#' 
#' Preparing the data for further analysis
#' 
# -------------------------------------------------------------------------


# load data ---------------------------------------------------------------

# load efficacy data:
studies_overview <- read_excel("raw-data/Studies passive Abs 2023-02-16.xlsx",sheet="Sheet1",col_names = TRUE)

# rename the columns:
tmp_col <- c(21,28,35,42)
tmp_names <- c("symptomatic","hospitalisation","ventilation","death") 
names(studies_overview)[tmp_col] <- tmp_names
names(studies_overview)[tmp_col+1] <- paste(tmp_names," CI lower",sep="")
names(studies_overview)[tmp_col+2] <- paste(tmp_names," CI upper",sep="")
names(studies_overview)[tmp_col+3] <- paste(tmp_names," events treatment",sep="")
names(studies_overview)[tmp_col+4] <- paste(tmp_names," n treatment",sep="")
names(studies_overview)[tmp_col+5] <- paste(tmp_names," events control",sep="")
names(studies_overview)[tmp_col+6] <- paste(tmp_names," n control",sep="")

# load meta-analysis:
load("raw-data/MetaAnalysisCalibration_Stanford_20230221_h09_m28_s22.RData")
Means_table$Antibodies <- tolower(Means_table$Antibodies)
conv <- Means_table$IC50[Means_table$Antibodies=="convalescent"]

# The code for the meta-analysis of the IC50 of different mAbs against different SARS-CoV-2 variants (using
# data from the Stanford database) that is included here, was previously published and can be accessed here:
# https://github.com/david-s-khoury/COVID19-mAb-prophylaxis
# This code can be used to create the data used here (MetaAnalysisCalibration_Stanford_20230221_h09_m28_s22.RData)
# and to reproduce Figure S10.


# include fold-convalescent dose ------------------------------------------

# add dose in mg for mAb treatments
dose_mg <- studies_overview$Dosage
dose_mg[studies_overview$`Treatment type`!="mAb"] <- NA
dose_mg <- ifelse(dose_mg=="2800mg + 2800mg","5600mg",ifelse(dose_mg=="700mg + 1400mg","2100mg",ifelse(dose_mg=="175mg+700mg+1400mg","2275mg",ifelse(dose_mg=="1000mg/1000mg","2000mg",
                                                                                                                                                     ifelse(dose_mg%in%c("7000mg or 700mg","1.2g, 2.4g, or 8.0g","1.2 or 2.4 g","40mg/kg, 80mg/kg","40mg/kg","80mg/kg","2.4g or 8g"),NA,dose_mg)))))
dose_mg <- as.numeric(gsub('[mg]','',dose_mg))

studies_overview <- cbind(studies_overview,dose_mg)

# compute the fold-convalescent dose:
tmp_mAbs <- studies_overview$Treatment
tmp_mAbs[studies_overview$`Treatment type`!="mAb"] <- NA # only mAb studies
tmp_variant <- studies_overview$`Covid variants`
tmp_variant[studies_overview$`Treatment type`!="mAb"] <- NA
tmp_variant[tmp_variant=="Omicron BA.1/BA.1.1"] <- "Omicron/BA.1.1"
tmp_variant[is.na(tmp_variant)] <- "Wild Type" # if there is no variant information, assume the variant is Wild-Type
tmp_IC50 <- rep(NA,length(tmp_mAbs))
for(i in 1:length(tmp_IC50)){
  if(!is.na(tmp_mAbs[i])){
    tmp_IC50[i] <- Means_table$IC50[Means_table$Antibodies==tmp_mAbs[i] & Means_table$VariantCategory==tmp_variant[i]]
  }
}
dose_conv <- (dose_mg/plasma_L)/(tmp_IC50*1e-3*conv)
# add the fold-convalescent dose for the Eom trial:
dose_conv[studies_overview$Trial=="Eom" & studies_overview$Dosage=="40mg/kg"] <- (40/plasma_conc)/(Means_table$IC50[Means_table$Antibodies=="regdanvimab" & Means_table$VariantCategory=="Wild Type"]*1e-3*conv)
dose_conv[studies_overview$Trial=="Eom" & studies_overview$Dosage=="80mg/kg"] <- (80/plasma_conc)/(Means_table$IC50[Means_table$Antibodies=="regdanvimab" & Means_table$VariantCategory=="Wild Type"]*1e-3*conv)


# fold-convalescent dose for CP studies -----------------------------------

# determine the fold-convalescent dose for CP (convalescent plasma) studies 
# that report hospitalisations after treatment of symptomatic patients

# Korley study:
# data about convalescent neutralizing Ab titers using Broad PRNT extracted from Di Germanio et al. (2021), Figure 2B BROAD PRNT (-> first measurement within 60 days for 15 donors)
DG_data <- read_csv("raw-data/Di Germanio et al (2021), Figure 2B Broad PRNT.csv",col_names = FALSE) # load data extracted from Di Germanio Figure 2B BROAD PRNT
data_tmp <- DG_data$X2
geom_mean <- exp(mean(log(data_tmp))) # geometric mean titer from Di Germanio
Korley_median <- 641/geom_mean # median titer 641 compared to the mean titer from Di Germanio
Korley_conv <- Korley_median*0.25/3 # dilution: 250mL of donor plasma in 3l of recipient plasma

# Libster study:
Lib_S4 <- read_csv("raw-data/Libster et al (2021), Figure S4.csv",col_names = FALSE) # load data extracted from Figure S4
donor_dist <- Lib_S4$X2
# fit a normal distribution to the donor distribution data using a likelihood approach
num_cens <- round(length(donor_dist)*0.72/0.28) # estimated number of censored data using that donors are in the 28th percentile
nllh <- function(p){-sum(log(dnorm(log10(donor_dist),p[1],exp(p[2]))))-num_cens*log(pnorm(log10(1000),p[1],exp(p[2])))}
donor_dist_fit <- nlm(nllh,c(2,3),hessian = TRUE)
# compute CIs for the mean and SD of the normal distribution:
mean_est <- donor_dist_fit$estimate[1]
sd_est <- exp(donor_dist_fit$estimate[2])
mean_CI <- mean_est+c(-1,1)*qnorm((1 + conf.level)/2)*sqrt(solve(donor_dist_fit$hessian)[1,1])
sd_CI <- exp(donor_dist_fit$estimate[2]+c(-1,1)*qnorm((1 + conf.level)/2)*sqrt(solve(donor_dist_fit$hessian)[2,2]))
# mean on a linear scale:
mean_est_lin <- 10^mean_est # -> use the convalescent mean 360 to compute the doses relative to convalescent!
mean_est_lin_CI <- 10^mean_CI
# median donor titer between 1,000 and 3,200:
Lib_med_titer_1 <- 10^(qnorm((pnorm(log10(1000),mean_est,sd_est)+pnorm(log10(3200),mean_est,sd_est))/2,mean_est,sd_est))
Lib_med_titer_1_conv <-(10^(qnorm((pnorm(log10(1000),mean_est,sd_est)+pnorm(log10(3200),mean_est,sd_est))/2,mean_est,sd_est))*0.25/3)/360
# median donor titer above 3,200:
Lib_med_titer_2 <- 10^(qnorm((1+pnorm(log10(3200),mean_est,sd_est))/2,mean_est,sd_est))
Lib_med_titer_2_conv <- (10^(qnorm((1+pnorm(log10(3200),mean_est,sd_est))/2,mean_est,sd_est))*0.25/3)/360

# Millat-Martinez study:
MM_conv <- (386/160)*0.25/3 # median titer 386 compared to a median CP titer of 160  (from Gharbharan 2021) and dilution of 250ml in 3 liters of plasma

# Sullivan study:
# Sullivan estimated that they used the top 60-70% of unselected titers => we use 65% and assume a normal distribution of the log10-titers
N <- 1e6 # number of samples
titers_all <- rnorm(N,0,sd)
titers_selected <- titers_all[titers_all>=quantile(titers_all,0.35)]
fold_conv <- 10^(mean(titers_selected)) # geometric mean fold-convalescent of top 65% of unselected donors
Sullivan_dose_conv <- 1.84*0.25/3 # include dilution of donor plasma

# add the fold-convalescent dose for some CP studies: 
dose_conv[studies_overview$Trial=="Korley"] <- Korley_conv
dose_conv[studies_overview$Trial=="Libster" & studies_overview$Subgroup=="CP titer (1,000-3,200)"] <- Lib_med_titer_1_conv
dose_conv[studies_overview$Trial=="Libster" & studies_overview$Subgroup=="CP titer (>3,200)"] <- Lib_med_titer_2_conv
dose_conv[studies_overview$Trial=="Millat-Martinez"] <- MM_conv
dose_conv[studies_overview$Trial=="Sullivan"] <- Sullivan_dose_conv

# add the fold-convalescent dose to the data
studies_overview <- cbind(studies_overview,dose_conv)

# save processed data: (includes information presented in Tables S1, S2, and S10)
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(studies_overview,glue::glue("output/Processed-data_{Sys.Date()}_{cur_time}.xlsx"))


# cleanup -----------------------------------------------------------------

rm(DG_data,donor_dist_fit,Lib_S4,data_tmp,donor_dist,dose_conv,dose_mg,fold_conv,
   geom_mean,i,Korley_conv,Korley_median,Lib_med_titer_1,Lib_med_titer_1_conv,
   Lib_med_titer_2,Lib_med_titer_2_conv,mean_CI,mean_est,mean_est_lin,mean_est_lin_CI,
   MM_conv,N,num_cens,sd_CI,sd_est,Sullivan_dose_conv,titers_all,tmp_col,tmp_IC50,
   tmp_mAbs,tmp_names,tmp_variant,nllh,titers_selected,cur_time)
