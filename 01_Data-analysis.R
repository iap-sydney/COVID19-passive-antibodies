# -------------------------------------------------------------------------
#' 
#' Analysis of the data for treatment efficacy of passive antibody treatments
#' 
# -------------------------------------------------------------------------


# number of symptomatic patients with outcome events ----------------------

studies_tmp <- studies_overview[is.na(studies_overview$Subgroup) & studies_overview$`Treatment stage`=="symptomatic",]
sympt_outcomes <- data.frame(all=data.frame(colSums(studies_tmp[,names(studies_tmp)%in%c("hospitalisation events treatment","hospitalisation n treatment","hospitalisation events control",
                                                                                         "hospitalisation n control","ventilation events treatment","ventilation n treatment",
                                                                                         "ventilation events control","ventilation n control","death events treatment","death n treatment",
                                                                                         "death events control","death n control")],na.rm = TRUE)),
                             mAb=data.frame(colSums(studies_tmp[studies_tmp$`Treatment type`=="mAb",names(studies_tmp)%in%c("hospitalisation events treatment","hospitalisation n treatment","hospitalisation events control",
                                                                                                                            "hospitalisation n control","ventilation events treatment","ventilation n treatment",
                                                                                                                            "ventilation events control","ventilation n control","death events treatment","death n treatment",
                                                                                                                            "death events control","death n control")],na.rm = TRUE)),
                             plasma=data.frame(colSums(studies_tmp[studies_tmp$`Treatment type`!="mAb",names(studies_tmp)%in%c("hospitalisation events treatment","hospitalisation n treatment","hospitalisation events control",
                                                                                                                               "hospitalisation n control","ventilation events treatment","ventilation n treatment",
                                                                                                                               "ventilation events control","ventilation n control","death events treatment","death n treatment",
                                                                                                                               "death events control","death n control")],na.rm = TRUE)))
names(sympt_outcomes) <- c("all","mAb","plasma")

# save: (Table S12)
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
openxlsx::write.xlsx(sympt_outcomes,glue::glue("output/Symptomatic-patients-numbers_{Sys.Date()}_{cur_time}.xlsx"),
                     rowNames=TRUE)

# cleanup:
rm(studies_tmp,sympt_outcomes,cur_time)


# Efficacy by dose model (therapeutic treatment) --------------------------

# select initial (treatment) stage and outcome stage for therapeutic treatment:
initial_stage <- "symptomatic"
outcome_stage <- "hospitalisation"

# select data with the right initial stage, only extract necessary columns
dose_resp_data <- studies_overview[studies_overview$`Treatment stage`==initial_stage,
                                   names(studies_overview)%in%c("Trial","Subgroup","Treatment","Treatment type","Dosage","dose_mg","dose_conv","Administration",
                                                                "Patient risk","hospitalisation","hospitalisation CI lower","hospitalisation CI upper",
                                                                "hospitalisation events treatment","hospitalisation n treatment",
                                                                "hospitalisation events control","hospitalisation n control",
                                                                "Time from onset of symptoms to drug (median, days)")]
# select dose subgroups, exclude other subgroups: 
dose_resp_data <- rbind(dose_resp_data[grepl("Dose",dose_resp_data$Subgroup),],
                        dose_resp_data[dose_resp_data$Trial%in%c("Dougan (bam+ete, 2.1g)","Gupta","Montgomery","Gottlieb (bam + ete)") & is.na(dose_resp_data$Subgroup),],
                        dose_resp_data[dose_resp_data$`Treatment type`=="CP" & is.na(dose_resp_data$Subgroup) & dose_resp_data$Trial!="Libster",],
                        dose_resp_data[dose_resp_data$Trial=="Libster" & !is.na(dose_resp_data$Subgroup),])
# rename studies to include mAb and dose:
dose_resp_data$Trial[dose_resp_data$`Treatment type`=="mAb"] <- paste(dose_resp_data$Trial[dose_resp_data$`Treatment type`=="mAb"]," (",dose_resp_data$Treatment[dose_resp_data$`Treatment type`=="mAb"],", ",as.numeric(gsub("mg","",dose_resp_data$Dosage[dose_resp_data$`Treatment type`=="mAb"]))/1000,"g)",sep="")
dose_resp_data$Trial[dose_resp_data$Trial=="Eom (regdanvimab, NAg)" & dose_resp_data$Subgroup=="Dose (40mg/kg)"] <- "Eom (regdanvimab, 40mg/kg)"
dose_resp_data$Trial[dose_resp_data$Trial=="Eom (regdanvimab, NAg)" & dose_resp_data$Subgroup=="Dose (80mg/kg)"] <- "Eom (regdanvimab, 80mg/kg)"
dose_resp_data$Trial[dose_resp_data$Trial=="Libster" & dose_resp_data$Subgroup=="CP titer (1,000-3,200)"] <- "Libster (CP 1,000-3,200)"
dose_resp_data$Trial[dose_resp_data$Trial=="Libster" & dose_resp_data$Subgroup=="CP titer (>3,200)"] <- "Libster (CP >,200)"
dose_resp_data$Trial[which(grepl(" \\(bam\\) ",dose_resp_data$Trial))] <- gsub(" \\(bam\\) "," ",dose_resp_data$Trial[which(grepl(" \\(bam\\) ",dose_resp_data$Trial))])
dose_resp_data$Trial[dose_resp_data$Trial=="Dougan (bam+ete, 2.1g) (bamlanivimab + etesevimab, NAg)"] <- "Dougan (bamlanivimab + etesevimab, 2.1g)"
dose_resp_data$Trial[dose_resp_data$Trial=="Gottlieb (bam + ete) (bamlanivimab + etesevimab, NAg)"] <- "Gottlieb (bamlanivimab + etesevimab, 5.6g)"

# exclude studies without fold-convalescent dose or hospitalisation outcomes:
dose_resp_data <- dose_resp_data[!is.na(dose_resp_data$dose_conv),]
dose_resp_data <- dose_resp_data[!is.na(dose_resp_data$hospitalisation),]

outcome_stage_col <- which(names(dose_resp_data)==outcome_stage)

names(dose_resp_data)[outcome_stage_col] <- "outcome"
names(dose_resp_data)[names(dose_resp_data)==paste(outcome_stage," CI lower",sep="")] <- "lower"
names(dose_resp_data)[names(dose_resp_data)==paste(outcome_stage," CI upper",sep="")] <- "upper"
names(dose_resp_data)[names(dose_resp_data)==paste(outcome_stage," events treatment",sep="")] <- "eTreat"
names(dose_resp_data)[names(dose_resp_data)==paste(outcome_stage," n treatment",sep="")] <- "nTreat"
names(dose_resp_data)[names(dose_resp_data)==paste(outcome_stage," events control",sep="")] <- "eCont"
names(dose_resp_data)[names(dose_resp_data)==paste(outcome_stage," n control",sep="")] <- "nCont"

# patient risk:
dose_resp_data$`Patient risk` <- factor(dose_resp_data$`Patient risk`,levels = c("low","mixed","high"))
names(dose_resp_data)[names(dose_resp_data)=="Patient risk"] <- "risk" 

# median time since symptom onset:
names(dose_resp_data)[names(dose_resp_data)=="Time from onset of symptoms to drug (median, days)"] <- "median_time" 
dose_resp_data$median_time[dose_resp_data$Trial=="Montgomery (cilgavimab + tixagevimab, 0.6g)"] <- 5 # median computed from Fig.2 data
dose_resp_data$median_time[dose_resp_data$Trial=="Gupta (sotrovimab, 0.5g)"] <- 2.5 # middle of 0 to 5 days range
dose_resp_data$median_time[dose_resp_data$Trial=="Libster (CP 1,000-3,200)"] <- mean(c(39.6,38.3))/24 # computed from reported mean hours since symptom onset
dose_resp_data$median_time[dose_resp_data$Trial=="Libster (CP >,200)"] <- mean(c(39.6,38.3))/24 # computed from reported mean hours since symptom onset

# Fit a mixed-effects models:
dose_resp_data_me <- data.frame(Trial=rep(dose_resp_data$Trial,2),events=c(dose_resp_data$eTreat,dose_resp_data$eCont),
                                n=c(dose_resp_data$nTreat,dose_resp_data$nCont),treatment=c(rep("treatment",nrow(dose_resp_data)),rep("control",nrow(dose_resp_data))),
                                risk=rep(dose_resp_data$risk,2),dose_conv=rep(dose_resp_data$dose_conv,2),median_time=c(rep(dose_resp_data$median_time,2)))
dose_resp_data_me <- dplyr::mutate(dose_resp_data_me,no_events=n-events)

# change the treatment variable to be coded as 0 for control and 1 for treatment:
dose_resp_data_me$treatment <- ifelse(dose_resp_data_me$treatment=="treatment",1,0)

# mixed-effects model for treatment effect by dose only:
tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment + log10(dose_conv):treatment + (1 | Trial),control = glmerControl(optimizer="bobyqa"),
                   family=binomial("log"), data=dose_resp_data_me)
est <- exp(fixef(tmp.glmer))
tmp.ci <- exp(confint(tmp.glmer))
test.sign <- drop1(tmp.glmer, scope=c("treatment","treatment:log10(dose_conv)"), test="Chisq")
results_dose <- list(est,tmp.ci,test.sign,tmp.glmer)

# save results: 
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
save(results_dose,file=glue::glue("output/Efficacy-by-dose-GLMM_{Sys.Date()}_{cur_time}.RData")) # Table S7


# mixed-effects model for dose and median time since symptom onset
tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment + median_time:treatment + log10(dose_conv):treatment + (1 | Trial),
                   control = glmerControl(optimizer="bobyqa"),family=binomial("log"), data=dose_resp_data_me)
est <- exp(fixef(tmp.glmer))
tmp.ci <- exp(confint(tmp.glmer))
test.sign <- drop1(tmp.glmer, scope=c("treatment","treatment:median_time","treatment:log10(dose_conv)"), test="Chisq")

results_dose_time <- list(est,tmp.ci,test.sign,tmp.glmer)

# save results: 
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
save(results_dose_time,file=glue::glue("output/Efficacy-by-dose-and-time-GLMM_{Sys.Date()}_{cur_time}.RData")) # Table S8


# cleanup:
rm(dose_resp_data,dose_resp_data_me,results_dose,results_dose_time,test.sign,tmp.ci,tmp.glmer,
   cur_time,est,initial_stage,outcome_stage,outcome_stage_col)
