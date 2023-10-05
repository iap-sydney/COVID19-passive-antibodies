# -------------------------------------------------------------------------
#' 
#' Effect of treatment stage and time on the efficacy of treatment
#' 
# -------------------------------------------------------------------------


# data preparation --------------------------------------------------------

# select rows such that no studies & data are included multiple times (by removing subgroup analysis rows)
# for O'Brien (PCR-) study: use data grouped by time to infection (1 week vs 2-4 weeks)
studies_tmp <- rbind(studies_overview[which(is.na(studies_overview$Subgroup) & studies_overview$Trial!="O'Brien (PCR-)"),],
                     studies_overview[which(studies_overview$Trial=="O'Brien (PCR-)" & grepl("Time",studies_overview$Subgroup)),])

# Analysis separately for mAb treatment and plasma (CP/hIVIG) treatment:
studies_tmp_ab <- studies_tmp[studies_tmp$`Treatment type`=="mAb",]
studies_tmp_plas <- studies_tmp[studies_tmp$`Treatment type`%in%c("CP","hIVIG"),]


# effect of the treatment stage for various outcome stages ----------------
# efficacy by outcome stage (lower parts of Tables S4 and S5)

outcome_stages_all <- c("symptomatic","hospitalisation","ventilation","death")
outcome_stages_all_cols <- which(names(studies_tmp)%in%outcome_stages_all)
results_ab <- data.frame(outcome.stage=c(),estimate=c(),pval=c(),var=c()) # save outcome stage, estimated relative risk, and p-value
results_plas <- data.frame(outcome.stage=c(),estimate=c(),pval=c(),var=c())
for(i in 1:length(outcome_stages_all)){
  tmp_outcome_stage <- outcome_stages_all[i]
  # select data for the current outcome stage (exclude studies with no outcomes results)
  data_tmp_ab <- studies_tmp_ab[!is.na(studies_tmp_ab[,outcome_stages_all_cols[i]]),]
  data_tmp_plas <- studies_tmp_plas[!is.na(studies_tmp_plas[,outcome_stages_all_cols[i]]),]
  # mAb data:
  if(nrow(data_tmp_ab)>0 & length(unique(data_tmp_ab$`Treatment stage`))>1){
    data_tmp_ab <- data.frame("Trial"=rep(data_tmp_ab$Trial,2),"events"=c(data_tmp_ab[,outcome_stages_all_cols[i]+3],data_tmp_ab[,outcome_stages_all_cols[i]+5]),
                              "n"=c(data_tmp_ab[,outcome_stages_all_cols[i]+4],data_tmp_ab[,outcome_stages_all_cols[i]+6]),
                              "treatment"=c(rep("treatment",nrow(data_tmp_ab)),rep("control",nrow(data_tmp_ab))),"initial.stage"=rep(data_tmp_ab$`Treatment stage`,2))
    data_tmp_ab$initial.stage <- factor(data_tmp_ab$initial.stage,levels = c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death"))
    data_tmp_ab <- dplyr::mutate(data_tmp_ab,no_events=n-events)
    # mixed effects model
    tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment*initial.stage + (1 | Trial),family=binomial("log"), data=data_tmp_ab)
    est <- exp(fixef(tmp.glmer))
    est <- est[grepl("treatment:initial",names(est))]
    tmp.ci <- exp(confint(tmp.glmer))
    tmp.ci <- tmp.ci[grepl("treatment:initial",rownames(tmp.ci)),]
    test.sign <- drop1(tmp.glmer, test="Chisq")
    tmp.var <- as.data.frame(VarCorr(tmp.glmer))$vcov
    if(is.null(nrow(tmp.ci))){
      results_ab <- rbind(results_ab,data.frame(outcome.stage=tmp_outcome_stage,estimate=est,CI_lower=tmp.ci[1],CI_upper=tmp.ci[2],pval=test.sign$`Pr(Chi)`[2],var=tmp.var))
    }else{
      results_ab <- rbind(results_ab,data.frame(outcome.stage=tmp_outcome_stage,estimate=est,CI_lower=tmp.ci[,1],CI_upper=tmp.ci[,2],pval=test.sign$`Pr(Chi)`[2],var=tmp.var))
    }
    
    # exclude peri-(post-)exposure treatment for outcome death #TODO!!!
    if(tmp_outcome_stage=="death"){
      data_tmp_ab <- data_tmp_ab[data_tmp_ab$initial.stage%in%c("symptomatic","hospitalisation"),]
      # mixed effects model
      tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment*initial.stage + (1 | Trial),family=binomial("log"), data=data_tmp_ab)
      est <- exp(fixef(tmp.glmer))
      est <- est[grepl("treatment:initial",names(est))]
      tmp.ci <- exp(confint(tmp.glmer))
      tmp.ci <- tmp.ci[grepl("treatment:initial",rownames(tmp.ci)),]
      test.sign <- drop1(tmp.glmer, test="Chisq")
      tmp.var <- as.data.frame(VarCorr(tmp.glmer))$vcov
      results_ab <- rbind(results_ab,data.frame(outcome.stage=tmp_outcome_stage,estimate=est,CI_lower=tmp.ci[1],CI_upper=tmp.ci[2],pval=test.sign$`Pr(Chi)`[2],var=tmp.var))
    }
    
    # model checking:
    # allfit_models <- allFit(tmp.glmer, maxfun = 30000)
    # summary(allfit_models)$fixef
  }
  # CP & IVIG data:
  if(nrow(data_tmp_plas)>0 & length(unique(data_tmp_plas$`Treatment stage`))>1){
    data_tmp_plas <- data.frame("Trial"=rep(data_tmp_plas$Trial,2),"events"=c(data_tmp_plas[,outcome_stages_all_cols[i]+3],data_tmp_plas[,outcome_stages_all_cols[i]+5]),
                                "n"=c(data_tmp_plas[,outcome_stages_all_cols[i]+4],data_tmp_plas[,outcome_stages_all_cols[i]+6]),
                                "treatment"=c(rep("treatment",nrow(data_tmp_plas)),rep("control",nrow(data_tmp_plas))),"initial.stage"=rep(data_tmp_plas$`Treatment stage`,2))
    data_tmp_plas$initial.stage <- factor(data_tmp_plas$initial.stage,levels = c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death"))
    data_tmp_plas <- dplyr::mutate(data_tmp_plas,no_events=n-events)
    # mixed effects model
    tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment*initial.stage + (1 | Trial),family=binomial("log"), data=data_tmp_plas)
    est <- exp(fixef(tmp.glmer))
    est <- est[grepl("treatment:initial",names(est))]
    tmp.ci <- exp(confint(tmp.glmer))
    tmp.ci <- tmp.ci[grepl("treatment:initial",rownames(tmp.ci)),]
    test.sign <- drop1(tmp.glmer, test="Chisq")
    tmp.var <- as.data.frame(VarCorr(tmp.glmer))$vcov
    results_plas <- rbind(results_plas,data.frame(outcome.stage=tmp_outcome_stage,estimate=est,CI_lower=tmp.ci[1],CI_upper=tmp.ci[2],pval=test.sign$`Pr(Chi)`[2],var=tmp.var))
    
    # model checking:
    # allfit_models <- allFit(tmp.glmer, maxfun = 30000)
    # summary(allfit_models)$fixef
  }
}
# results_ab
# results_plas

# save results:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(results_ab,glue::glue("output/Effect of timing by outcome stage (mAb)_{Sys.Date()}_{cur_time}.xlsx")) # Table S4 (lower part)
write_xlsx(results_plas,glue::glue("output/Effect of timing by outcome stage (CP, hIG)_{Sys.Date()}_{cur_time}.xlsx")) # Table S5 (lower part)


# effect of the treatment stage for preventing progression ----------------
# progression to the next stage instead of progression to a specific outcome

# make progression data:
initial_stages_all <- unique(studies_tmp$`Treatment stage`)
outcome_stages_all <- c("symptomatic","hospitalisation","ventilation","death")
outcome_stages_all_cols <- which(names(studies_tmp)%in%outcome_stages_all)

# Analysis of Ab data:
progr_any_ab <- data.frame(Trial=c(),events=c(),n=c(),treatment=c(),initial.stage=c(),outcome.stage=c())
for(i in 1:nrow(studies_tmp_ab)){
  trial_tmp <- studies_tmp_ab$Trial[i]
  if(any(!is.na(studies_tmp_ab[i,outcome_stages_all_cols]))){
    outcome_tmp <- outcome_stages_all[which(!is.na(studies_tmp_ab[i,outcome_stages_all_cols]))]
    outcome_tmp_cols <- outcome_stages_all_cols[which(!is.na(studies_tmp_ab[i,outcome_stages_all_cols]))]
    for(j in 1:length(outcome_tmp)){
      progr_tmp <- data.frame(Trial=rep(trial_tmp,2),events=c(studies_tmp_ab[i,outcome_tmp_cols[j]+3],studies_tmp_ab[i,outcome_tmp_cols[j]+5]),
                              n=c(studies_tmp_ab[i,outcome_tmp_cols[j]+4],studies_tmp_ab[i,outcome_tmp_cols[j]+6]),
                              treatment=c("treatment","control"),initial.stage=rep(studies_tmp_ab$`Treatment stage`[i],2),outcome.stage=rep(outcome_tmp[j],2))
      progr_any_ab <- rbind(progr_any_ab,progr_tmp)
    }
  }
}
progr_any_ab <- dplyr::mutate(progr_any_ab,no_events=n-events)

# next stage
progr_next_ab <- progr_any_ab[which((progr_any_ab$initial.stage%in%c("pre-exposure","peri-(post-)exposure") & progr_any_ab$outcome.stage=="symptomatic") |
                                      (progr_any_ab$initial.stage=="symptomatic" & progr_any_ab$outcome.stage=="hospitalisation") |
                                      (progr_any_ab$initial.stage=="hospitalisation" & progr_any_ab$outcome.stage=="death")),]
# with initial stage as a numerical variable:
progr_next_num_ab <- dplyr::mutate(progr_next_ab,initial.stage.num = ifelse(progr_next_ab$initial.stage=="pre-exposure",0,
                                                                            ifelse(progr_next_ab$initial.stage=="peri-(post-)exposure",1,
                                                                                   ifelse(progr_next_ab$initial.stage=="symptomatic",2,3))))

# final model reported in the supplement:
tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment*initial.stage.num + (1 | Trial),family=binomial("log"), data=progr_next_num_ab)
est <- exp(fixef(tmp.glmer))
tmp.ci <- exp(confint(tmp.glmer))
test.sign <- drop1(tmp.glmer, scope = c("treatment","initial.stage.num","treatment:initial.stage.num"), test="Chisq")
rand.var <- as.data.frame(VarCorr(tmp.glmer))$vcov
results_progr_next_num_ab <- list(est,tmp.ci,test.sign,tmp.glmer,rand.var)

# save result:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
save(results_progr_next_num_ab,file=glue::glue("output/Effect of timing by stage progression (mAb)_{Sys.Date()}_{cur_time}.RData")) # Table S4 (upper part)


# Analysis of CP & IVIG data:
progr_any_plas <- data.frame(Trial=c(),events=c(),n=c(),treatment=c(),initial.stage=c(),outcome.stage=c())
for(i in 1:nrow(studies_tmp_plas)){
  trial_tmp <- studies_tmp_plas$Trial[i]
  if(any(!is.na(studies_tmp_plas[i,outcome_stages_all_cols]))){
    outcome_tmp <- outcome_stages_all[which(!is.na(studies_tmp_plas[i,outcome_stages_all_cols]))]
    outcome_tmp_cols <- outcome_stages_all_cols[which(!is.na(studies_tmp_plas[i,outcome_stages_all_cols]))]
    for(j in 1:length(outcome_tmp)){
      progr_tmp <- data.frame(Trial=rep(trial_tmp,2),events=c(studies_tmp_plas[i,outcome_tmp_cols[j]+3],studies_tmp_plas[i,outcome_tmp_cols[j]+5]),
                              n=c(studies_tmp_plas[i,outcome_tmp_cols[j]+4],studies_tmp_plas[i,outcome_tmp_cols[j]+6]),
                              treatment=c("treatment","control"),initial.stage=rep(studies_tmp_plas$`Treatment stage`[i],2),outcome.stage=rep(outcome_tmp[j],2))
      progr_any_plas <- rbind(progr_any_plas,progr_tmp)
    }
  }
}
progr_any_plas <- dplyr::mutate(progr_any_plas,no_events=n-events)

# next stage
progr_next_plas <- progr_any_plas[which((progr_any_plas$initial.stage%in%c("pre-exposure","peri-(post-)exposure") & progr_any_plas$outcome.stage=="symptomatic") |
                                          (progr_any_plas$initial.stage=="symptomatic" & progr_any_plas$outcome.stage=="hospitalisation") |
                                          (progr_any_plas$initial.stage=="hospitalisation" & progr_any_plas$outcome.stage=="death")),]
# with initial stage as a numerical variable:
progr_next_num_plas <- dplyr::mutate(progr_next_plas,initial.stage.num = ifelse(progr_next_plas$initial.stage=="pre-exposure",1,
                                                                                ifelse(progr_next_plas$initial.stage=="peri-(post-)exposure",2,
                                                                                       ifelse(progr_next_plas$initial.stage=="symptomatic",3,4))))

# final model reported in the supplement:
tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment*initial.stage.num + (1 | Trial),family=binomial("log"), data=progr_next_num_plas)
est <- exp(fixef(tmp.glmer))
tmp.ci <- exp(confint(tmp.glmer))
test.sign <- drop1(tmp.glmer, scope=c("treatment","initial.stage.num","treatment:initial.stage.num"), test="Chisq")
rand.var <- as.data.frame(VarCorr(tmp.glmer))$vcov
results_progr_next_num_plas <- list(est,tmp.ci,test.sign,tmp.glmer,rand.var)

# save result:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
save(results_progr_next_num_plas,file=glue::glue("output/Effect of timing by stage progression (CP, hIG)_{Sys.Date()}_{cur_time}.RData")) # Table S5 (upper part)



# visualization of effect of time since symptom onset ---------------------
# efficacy by time of treatment for therapeutic treatment (symptomatic to hospitalisation)

# Extract relevant data, i.e. studies of treatment of symptomatic patients that report hospitalisation outcomes:
studies_tmp <- studies_overview[studies_overview$`Treatment stage`=="symptomatic" & !is.na(studies_overview$hospitalisation),
                                names(studies_overview)%in%c("Trial","Subgroup","Treatment","Treatment type","Dosage","Treatment stage",
                                                             "Time from onset of symptoms to drug (median, days)","Time from onset of symptoms to drug (IQR or range, min, days)",
                                                             "Time from onset of symptoms to drug (IQR or range, max, days)","IQR or range",
                                                             "hospitalisation","hospitalisation CI lower","hospitalisation CI upper",
                                                             "hospitalisation events treatment","hospitalisation n treatment",
                                                             "hospitalisation events control","hospitalisation n control")]
# Remove subgroups except time of symptom onset subgroups:
studies_tmp <- studies_tmp[is.na(studies_tmp$Subgroup) | grepl("Symptom onset",studies_tmp$Subgroup),]
studies_tmp <- studies_tmp[!(studies_tmp$Trial%in%studies_tmp$Trial[grepl("Symptom onset",studies_tmp$Subgroup)] & is.na(studies_tmp$Subgroup)),] # remove the studies whose time since symmpton onset is included
studies_tmp$Trial[grepl("Symptom onset",studies_tmp$Subgroup)] <- paste(studies_tmp$Trial[grepl("Symptom onset",studies_tmp$Subgroup)],gsub("Symptom onset","",studies_tmp$Subgroup[grepl("Symptom onset",studies_tmp$Subgroup)]),sep="")

# add lower bound for the range of days since symptom onset
studies_tmp$`Time from onset of symptoms to drug (IQR or range, min, days)`[studies_tmp$`IQR or range`=="range" & is.na(studies_tmp$`Time from onset of symptoms to drug (IQR or range, min, days)`)] <- 0

# add missing median information:
studies_tmp <- dplyr::mutate(studies_tmp,median_known=!is.na(studies_tmp$`Time from onset of symptoms to drug (median, days)`))
studies_tmp$`Time from onset of symptoms to drug (median, days)`[!studies_tmp$median_known] <- (studies_tmp$`Time from onset of symptoms to drug (IQR or range, min, days)`[!studies_tmp$median_known]+
                                                                                                  studies_tmp$`Time from onset of symptoms to drug (IQR or range, max, days)`[!studies_tmp$median_known])/2
# add manually computed medians for Montgomery study:
studies_tmp$`Time from onset of symptoms to drug (median, days)`[studies_tmp$Trial=="Montgomery (up to 5 days)"] <- 4.5
studies_tmp$`Time from onset of symptoms to drug (median, days)`[studies_tmp$Trial=="Montgomery (more than 5 days)"] <- 7
# add computed mean for Libster study:
studies_tmp$`Time from onset of symptoms to drug (median, days)`[studies_tmp$Trial=="Libster"] <- mean(c(39.6,38.3))/24 # reported mean times since symptom onset for treatment and convalescent group

timing_aes <- data.frame(Trial=sort(unique(studies_tmp$Trial)))
timing_aes <- cbind(timing_aes,line_type=rep("solid",nrow(timing_aes)),shape=rep(0,nrow(timing_aes)))
timing_aes$line_type[timing_aes$Trial%in%studies_tmp$Trial[studies_tmp$`IQR or range`=="IQR"]] <- "dashed"
timing_aes$shape <- c(0,2,5,6,13,16,17,18,1,2,2,4,7)
timing_aes$shape[timing_aes$Trial%in%studies_tmp$Trial[!studies_tmp$median_known] & !grepl("Montgomery",timing_aes$Trial)] <- 3

# Does the median time since symptom onset influence efficacy? -> include a GLMM
timing.data <- data.frame(Trial=rep(studies_tmp$Trial,2),events=c(studies_tmp$`hospitalisation events treatment`,studies_tmp$`hospitalisation events control`),
                          n=c(studies_tmp$`hospitalisation n treatment`,studies_tmp$`hospitalisation n control`),
                          treatment=c(rep("treatment",nrow(studies_tmp)),rep("control",nrow(studies_tmp))),
                          median_time=rep(studies_tmp$`Time from onset of symptoms to drug (median, days)`,2))
# timing.data <- timing.data[timing.data$Trial%in%studies_tmp$Trial[studies_tmp$median_known | grepl("Montgomery",timing.data$Trial)],]
timing.data <- dplyr::mutate(timing.data,no_events=n-events)
timing.data$treatment <- ifelse(timing.data$treatment=="treatment",1,0)

# tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment + treatment:median_time + (1 | Trial),family=binomial("log"), data=timing.data) # no median time covariate in the model
tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment*median_time + (1 | Trial),family=binomial("log"), data=timing.data) # median time as a covariate in the model
est <- exp(fixef(tmp.glmer))
tmp.ci <- exp(confint(tmp.glmer))
test.sign <- drop1(tmp.glmer, scope=c("treatment","median_time","treatment:median_time"), test="Chisq")
rand.var <- as.data.frame(VarCorr(tmp.glmer))$vcov
results_timing <- list(est,tmp.ci,test.sign,tmp.glmer,rand.var)

# save results: 
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
save(results_timing,file=glue::glue("output/Effect of symptom onset time (mAb, CP)_{Sys.Date()}_{cur_time}.RData")) # Table S6 (upper part)


# visualize the efficacy by time since symptom onset with the GLMM estimates
# efficacy estimate by days since symptom onset from the above model:
eff_est <- function(d){100*(1-exp(fixef(tmp.glmer)[2]+fixef(tmp.glmer)[4]*d))}
days <- seq(0,10,length.out=1e3)
estimate <- eff_est(days)
eff_est_data <- data.frame(days,estimate)

my_colours <- c("yellow3","goldenrod1","red","red3","lightpink","orchid3","dodgerblue","turquoise1","aquamarine3","mediumpurple","mediumpurple",
                "springgreen4","purple4")

FigS11A <- ggplot(data=studies_tmp,aes(x=`Time from onset of symptoms to drug (median, days)`,y=hospitalisation,color=Trial,shape=Trial)) +
  # visualize the data:
  # geom_hline(yintercept=0, color="#999999") + # bolder line at 0
  geom_errorbarh(aes(xmin = `Time from onset of symptoms to drug (IQR or range, min, days)`,
                     xmax = `Time from onset of symptoms to drug (IQR or range, max, days)`,
                     linetype = `IQR or range`, height=4), size=0.7) + # horizontal line between stages
  geom_errorbar(aes(ymin = `hospitalisation CI lower`, ymax = `hospitalisation CI upper`, width=0.3), size=0.7) + # vertical line for CIs
  geom_point(size=2) + 
  scale_shape_manual(values=timing_aes$shape) +
  scale_color_manual(values=my_colours) +
  scale_linetype_manual(values=timing_aes$line_type) +
  # add fit to the data:
  geom_line(data=eff_est_data,inherit.aes = FALSE, aes(x=days, y=estimate), size=2) +
  # axis, labels, etc.:
  labs(x="Median time from onset of symptoms [days]", y="Efficacy [%]", title = "Efficacy of therapeutic treatment by time since symptom onset") +
  scale_y_continuous(minor_breaks = c(seq(0, 100, 20)), breaks=c(seq(0, 100, 20)),labels=c(seq(0, 100, 20))) +
  coord_cartesian(ylim=c(0,100),xlim=c(0,10), expand = c(0.01)) + 
  # ylim(0,100) + 
  theme_bw() + theme(panel.grid = element_line(colour="gray95",size = 0.2))
# FigS11A

# save figure:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
pdf(glue::glue("output/FigS11A_{Sys.Date()}_{cur_time}.pdf"),height=4,width=7)
print(FigS11A)
dev.off()


# Same analysis for mAbs only:
CP_studies <- c("Korley","Libster","Millat-Martinez","Sullivan")
timing.data.mAbs <- timing.data[!(timing.data$Trial%in%CP_studies),]

tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment*median_time + (1 | Trial),family=binomial("log"), data=timing.data.mAbs) # median time as a covariate in the model
est <- exp(fixef(tmp.glmer))
tmp.ci <- exp(confint(tmp.glmer))
test.sign <- drop1(tmp.glmer, scope=c("treatment","median_time","treatment:median_time"), test="Chisq")
rand.var <- as.data.frame(VarCorr(tmp.glmer))$vcov

eff_est <- function(d){100*(1-exp(fixef(tmp.glmer)[2]+fixef(tmp.glmer)[4]*d))}
days <- seq(0,10,length.out=1e3)
estimate <- eff_est(days)
eff_est_data_mAbs <- data.frame(days,estimate)
results_timing_mAbs <- list(est,tmp.ci,test.sign,tmp.glmer,rand.var)

# save results: 
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
save(results_timing_mAbs,file=glue::glue("output/Effect of symptom onset time (mAb only)_{Sys.Date()}_{cur_time}.RData")) # Table S6 (lower part)

# visualization:
my_colours <- c("yellow3","goldenrod1","red","red3","lightpink","orchid3","mediumpurple","mediumpurple","purple4")
my_shapes <- c(0,2,5,6,13,NA,2,2,7)
my_lines <- rep("dashed",nrow(timing.data.mAbs)/2)
my_lines[c(2,6,7,8)] <- "solid"

FigS11B <- ggplot(data=studies_tmp[!(studies_tmp$Trial%in%CP_studies),],aes(x=`Time from onset of symptoms to drug (median, days)`,
                                                                          y=hospitalisation,color=Trial,shape=Trial)) +
  # visualize the data:
  # geom_hline(yintercept=0, color="#999999") + # bolder line at 0
  geom_errorbarh(aes(xmin = `Time from onset of symptoms to drug (IQR or range, min, days)`,
                     xmax = `Time from onset of symptoms to drug (IQR or range, max, days)`,
                     linetype = `IQR or range`, height=4), size=0.7) + # horizontal line between stages
  geom_errorbar(aes(ymin = `hospitalisation CI lower`, ymax = `hospitalisation CI upper`, width=0.3), size=0.7) + # vertical line for CIs
  geom_point(size=2) + 
  scale_shape_manual(values=my_shapes) +
  scale_color_manual(values=my_colours) +
  scale_linetype_manual(values=my_lines) +
  # add fit to the data:
  geom_line(data=eff_est_data_mAbs,inherit.aes = FALSE, aes(x=days, y=estimate), size=2) +
  # axis, labels, etc.:
  labs(x="Median time from onset of symptoms [days]", y="Efficacy [%]", title = "Efficacy of therapeutic treatment by time since symptom onset") +
  scale_y_continuous(minor_breaks = c(seq(0, 100, 20)), breaks=c(seq(0, 100, 20)),labels=c(seq(0, 100, 20))) +
  coord_cartesian(ylim=c(0,100),xlim=c(0,10), expand = c(0.01)) + 
  # ylim(0,100) + 
  theme_bw() + theme(panel.grid = element_line(colour="gray95",size = 0.2))
# FigS11B

# save figure:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
pdf(glue::glue("output/FigS11B_{Sys.Date()}_{cur_time}.pdf"),height=4,width=7)
print(FigS11B)
dev.off()

# combine subfigures and save figure:
FigS11 <- FigS11A | FigS11B
# FigS11

cur_time <- format(Sys.time(),"h%H_m%M_s%S")
pdf(glue::glue("output/FigS11_{Sys.Date()}_{cur_time}.pdf"),height=4,width=14)
print(FigS11)
dev.off()

# combine and save efficacy by day data:
eff_est_data <- data.frame(day = eff_est_data$days, fit_all = eff_est_data$estimate, fit_mAbs = eff_est_data_mAbs$estimate)

cur_time <- format(Sys.time(),"h%H_m%M_s%S")
save(eff_est_data,file=glue::glue("output/Estimated efficacy by median treatment time_{Sys.Date()}_{cur_time}.RData"))


# Maximal efficacy in a non-naive population with scaling -----------------

# Same analysis as for Fig. 4B but using the scaled dose-response curve instead:

# scaling factor for the dose-response curve: efficacy for mAb treatment on the day of symptom onset
new_max <- eff_est_data$fit_mAbs[eff_est_data$day==0]/100

# load estimated parameters for the dose-response curve:
load("raw-data/Dose-resp parameters (symp-hosp).RData")

# efficacy function with variable parameters:
eff <- function(d,par){new_max/(1+exp(-exp(par[2])*(log10(d)-log10(exp(par[1])))))}
# efficacy function with best fit parameters:
eff_par <- function(d){new_max/(1+exp(-exp(par_est$estimate[par_est$description=="slope"])*(log10(d)-log10(exp(par_est$estimate[par_est$description=="half-maximal dose"])))))}
# function of the initial Ab dose (d):
max_eff_non_naive <- function(d){max(c(0,(par_est$estimate[par_est$description=="maximal efficacy"]-eff_par(d))/(1-eff_par(d))))}
# function of the initial Ab dose (d) with an x-fold drop in the neutralisation titre against a variant:
# max_eff_non_naive_var <- function(d,x){(par_est$estimate[par_est$description=="maximal efficacy"]-eff_par(d/x))/(1-eff_par(d/x))}
# function to compute dose for given efficacy (same parameters as eff function, with the same scaling):
dose_from_eff <- function(eff,par){exp(par[1])*10^(-log(new_max/eff-1)/exp(par[2]))}
# dose by maximal efficacy in non-naive patients:
dose_from_max_eff <- function(max_eff,par){dose_from_eff((par_est$estimate[par_est$description=="maximal efficacy"]-max_eff)/(1-max_eff),par)}

dose <- 10^seq(-2,4,length.out = 1e3)
model_fit_scaled <- sapply(dose,function(x){100*max_eff_non_naive(x)})

# make 95% CIs using parametric bootstrapping:
N <- 1e5 # number of bootstraps
n_par <- nrow(par_est)
par_bootstr <- rmvnorm(n=N,mean=par_est$estimate,sigma=covar)
# compute efficacy CIs for a given dose:
eff_tmp <- function(d,par){(par[3]+new_max-par_est$estimate[par_est$description=="maximal efficacy"])/(1+exp(-exp(par[2])*(log10(d)-log10(exp(par[1])))))}
max_eff_non_naive_ci <- function(d){quantile(sapply(c(1:nrow(par_bootstr)),function(x){max(c(0,(par_bootstr[x,3]-eff_tmp(d,par_bootstr[x,]))/
                                                                                               (1-eff_tmp(d,par_bootstr[x,]))))}), probs = c((1-conf.level)/2,(1+conf.level)/2))}

model_fit_ci <- sapply(dose,function(x){100*max_eff_non_naive_ci(x)})

# combine all data & save it:
model_data_resc <- data.frame(dose,model_fit_scaled,t(model_fit_ci))#,model_fit_var)
names(model_data_resc) <- c("dose","model_fit","model_fit_lower","model_fit_upper")#),"model_fit_var")

# cur_time <- format(Sys.time(),"h%H_m%M_s%S")
# save(model_data_resc, file =glue::glue("output/Max efficacy non-naive population (rescaled)_{Sys.Date()}_{cur_time}.RData"))

# horizontal line position:
hline <- 30 # horizontal line at 30% efficacy

# neutralisation levels from Stanford database (data extracted 22-02-2023) for "two or more shots, with history" for different variants,
# average for BNT and Mod 2-6M:
var_names <- c("wild type","alpha","beta","gamma","delta","BA.1","BA.1.1","BA.2","BA.2.12.1","BA.2.75.2","XBB.1","XBB.1.5","BA.4/5","BA.4.6","BQ.1.1")
neuts <- c(2976.5,1502,637,1499,14899.5,719,414,1205,732,154,179.5,220,650.5,612,278.5)
neuts_conv <- neuts/conv
neut.data.2vax.hist <- data.frame(variant=var_names,neuts=neuts,neuts_conv=neuts_conv,efficacy=sapply(neuts_conv,function(x){100*max_eff_non_naive(x)}))

# add only some variants to the plot:
add_var <- c("wild type","alpha","BA.1","BA.2","XBB.1","BA.4/5","BQ.1.1")
add_neuts_data <- neut.data.2vax.hist[neut.data.2vax.hist$variant%in%add_var,]

# my_colours <- c("yellow2","yellow3","goldenrod1","darkorange","red","red3","orangered4","lightpink","hotpink","orchid3",
#                          "dodgerblue","turquoise1","palegreen2","aquamarine3","mediumpurple","springgreen4","darkviolet","purple4","black")
my_colours <- c("yellow3","goldenrod1","red","orangered4","hotpink","black","darkviolet")
my_shapes <- c(15,16,17,18,15,16,17,18)

# split the data to plot in below and above the threshold in different colours:
model_data_below <- model_data_resc[model_data_resc$model_fit_lower<=hline,]
model_data_below$model_fit_upper[model_data_below$model_fit_upper>=hline] <- hline
model_data_above <- model_data_resc[model_data_resc$model_fit_upper>=hline,]
model_data_above$model_fit_lower[model_data_above$model_fit_lower<=hline] <- hline

FigS14 <- ggplot() +
  # maximal efficacy below the threshold:
  geom_ribbon(data=model_data_below,aes(x=dose,ymin=model_fit_lower, ymax=model_fit_upper),fill="black", alpha = 0.1)+ # TODO
  geom_line(data=model_data_resc[c(which.min(abs(model_data_resc$model_fit-hline)):nrow(model_data_resc)),],aes(x=dose,y=model_fit),color="black", size=1) +
  # maximal efficacy above the threshold:
  geom_ribbon(data=model_data_above,aes(x=dose,ymin=model_fit_lower, ymax=model_fit_upper),fill="dodgerblue", alpha = 0.1)+ # TODO
  geom_line(data=model_data_resc[model_data_resc$model_fit>=hline,],aes(x=dose,y=model_fit),color="dodgerblue", size=1) +
  # add horizonal and vertical lines with annotations:
  geom_segment(aes(x=1e-3,xend=1e4,y=hline,yend=hline),color="red",linetype="dashed") +
  geom_segment(aes(x=dose_from_max_eff(hline/100,par_est$estimate),xend=dose_from_max_eff(hline/100,par_est$estimate),
                   y=hline,yend=0),color="red",linetype="dashed") +
  annotate(geom="text",x=0.1,y=hline,label=paste(hline,"% Efficacy \n",sep=""),color="red") + 
  annotate(geom="text",x=dose_from_max_eff(hline/100,par_est$estimate),y=hline/2,color="red",angle=90,
           label=paste("Neut. for ",hline,"% \n efficacy: ",round(dose_from_max_eff(hline/100,par_est$estimate),2),sep="")) + 
  # add points for different variants:
  geom_point(data=add_neuts_data, aes(x=neuts_conv,y=efficacy,color=variant,shape=variant),size=3) +
  scale_color_manual(values=my_colours) +
  scale_shape_manual(values=my_shapes) +
  # other annotations:
  geom_segment(aes(x=add_neuts_data$neuts_conv[add_neuts_data$variant=="wild type"],xend=add_neuts_data$neuts_conv[add_neuts_data$variant=="wild type"],
                   y=0,yend=100), color="black", linetype="dashed") +
  # annotate(geom="text",x=50,y=15,color="black",label="< 30% clinical efficacy of passive \n Ab treatment expected") + 
  # annotate(geom="text",x=0.1,y=85,color="black",label="Potential benefit of \n passive Ab treatment \n (if an efficacious Ab \n is available") + 
  # add previous curve for comparison:
  # geom_line(data=model_data,aes(x=dose,y=model_fit),color="springgreen4",linetype="dotted", size=2) + # load the data for the previous curve (Fig.4B) from output: load("output/Max efficacy non-naive population_DATE_TIME.RData")
  # axis, theme, etc.:
  scale_x_log10(breaks=c(1e-2,1e-1,1e0,1e1,1e2,1e3),labels=c("0.01","0.1","1","10","100","1000")) +
  scale_y_continuous(minor_breaks = seq(0, 100, 20),breaks=seq(0, 100, 20),labels=seq(0, 100, 20),expand = c(1e-2,1e-2)) +
  coord_cartesian(ylim=c(-2,101),xlim=c(1e-2,1e4),expand=c(0)) +
  labs(x="Endogenous neutralisation level [fold convalescent]", y="Efficacy [%]", 
       title = "Maximal efficacy of passive Abs in a non-naive population") +
  theme_bw() +
  theme(legend.key.size = unit(0.6, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2))
# FigS14

# save the figure:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
pdf(glue::glue("output/FigS14_{Sys.Date()}_{cur_time}.pdf"),height=4,width=6)
print(FigS14)
dev.off()


# cleanup -----------------------------------------------------------------

rm(add_neuts_data,covar,data_tmp_ab,data_tmp_plas,eff_est_data,eff_est_data_mAbs,FigS11,FigS11A,
   FigS11B,FigS14,model_data_above,model_data_below,model_data_resc,model_fit_ci,neut.data.2vax.hist,par_bootstr,par_est,progr_any_ab,progr_any_plas,progr_next_ab,progr_next_num_ab,
   progr_next_num_plas,progr_next_plas,progr_tmp,results_ab,results_plas,
   results_progr_next_num_ab,results_progr_next_num_plas,results_timing,
   results_timing_mAbs,studies_tmp,studies_tmp_ab,studies_tmp_plas,table_eff,table_eff_fixed_m,test.sign,
   timing_aes,timing.data,timing.data.mAbs,tmp.ci,tmp.glmer,add_var,CP_studies,cur_time,
   days,dose,est,estimate,hline,i,initial_stages_all,j,model_fit_scaled,my_colours,my_lines,
   my_shapes,N,n_par,neuts,neuts_conv,new_max,outcome_stages_all,outcome_stages_all_cols,outcome_tmp,outcome_tmp_cols,
   rand.var,tmp_outcome_stage,tmp.var,trial_tmp,var_names,dose_from_eff,dose_from_max_eff,eff,eff_est,eff_par,eff_tmp,max_eff_non_naive,max_eff_non_naive_ci)
