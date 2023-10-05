# -------------------------------------------------------------------------
#' 
#' Dose-response curve for efficacy of therapeutic treatment
#' 
# -------------------------------------------------------------------------


# data preparation --------------------------------------------------------

# select initial (=treatment) and outcome stage
initial_stage <- "symptomatic"
outcome_stage <- "hospitalisation"

# select data with the right initial stage, only extract necessary columns
dose_resp_data <- studies_overview[studies_overview$`Treatment stage`==initial_stage,
                                   names(studies_overview)%in%c("Trial","Subgroup","Treatment","Treatment type","Dosage","dose_mg","dose_conv","Administration",
                                                                "Patient risk","hospitalisation","hospitalisation CI lower","hospitalisation CI upper",
                                                                "hospitalisation events treatment","hospitalisation n treatment",
                                                                "hospitalisation events control","hospitalisation n control")]
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


# fitting a logistic dose-response curve to the data ----------------------

# efficacy function by dose (d):
# logistic function with 3 parameters: half-maximal dose (EC50, log-transformed), slope (log-transformed), and maximal efficacy
eff <- function(d,par){par[3]/(1+exp(-exp(par[2])*(log10(d)-log10(exp(par[1])))))}
# function to compute dose for given efficacy (same parameters as eff function, with the same scaling):
dose_from_eff <- function(eff,par){exp(par[1])*10^(-log(par[3]/eff-1)/exp(par[2]))}

# negative loglikelihood function:
# parameters: baseline risk for each study and 3 parameters for efficacy by dose (EC50, slope, maximum)
n_studies <- nrow(dose_resp_data) # number of studies
nllh <- function(par){-sum(log(dbinom(dose_resp_data$eCont,dose_resp_data$nCont,par[1:n_studies])))-
    sum(log(dbinom(dose_resp_data$eTreat,dose_resp_data$nTreat,par[1:n_studies]*(1-eff(dose_resp_data$dose_conv,par[n_studies+(1:3)])))))}
# initial parameters: baseline risk = risk for control population, EC50=0, slope=1, maximal efficacy=0.9 (=90%)
par_init <- c(dose_resp_data$eCont/dose_resp_data$nCont,0,log(1),0.9)

# minimize the negative loglikelihood function:
dose_response_fit <- nlm(nllh,par_init,hessian = TRUE)
dose_response_par_fit <- dose_response_fit$estimate
dose_response_par_fit_ci <- cbind(dose_response_par_fit-qnorm((1+0.95)/2)*sqrt(diag(solve(dose_response_fit$hessian))),
                                  dose_response_par_fit+qnorm((1+0.95)/2)*sqrt(diag(solve(dose_response_fit$hessian))))
# best fit parameters for the efficacy function
par_est <- data.frame(description=c("half-maximal dose","slope","maximal efficacy"),
                      estimate=tail(dose_response_par_fit,3),CI_lower=tail(dose_response_par_fit_ci[,1],3),
                      CI_upper=tail(dose_response_par_fit_ci[,2],3))
par_est <- cbind(par_est,estimate_resc=c(exp(par_est$estimate[c(1:2)]),100*par_est$estimate[3]),
                 CI_lower_resc=c(exp(par_est$CI_lower[c(1:2)]),100*par_est$CI_lower[3]),
                 CI_upper_resc=c(exp(par_est$CI_upper[c(1:2)]),100*par_est$CI_upper[3]))


# compute CIs using parametric bootstrapping ------------------------------

# compute 95% CIs for the fit using parametric bootstrap of the parameters for the efficacy curve:
N <- 1e5 # number of bootstraps
n_par <- nrow(par_est)
covar <- solve(dose_response_fit$hessian)[nrow(dose_response_fit$hessian)-c((n_par-1):0),ncol(dose_response_fit$hessian)-c((n_par-1):0)] # covariance matrix for efficacy parameters
par_bootstr <- rmvnorm(n=N,mean=dose_response_par_fit[length(dose_response_par_fit)-c((n_par-1):0)],sigma=covar)

# compute efficacy CIs for a given dose:
eff_ci_from_dose <- function(d){quantile(sapply(c(1:nrow(par_bootstr)),function(x){eff(d,par_bootstr[x,])}), probs = c((1-conf.level)/2,(1+conf.level)/2))}
# compute dose CIs for a given efficacy
dose_ci_from_eff <- function(ef){quantile(sapply(c(1:nrow(par_bootstr)),function(x){dose_from_eff(ef,par_bootstr[x,])}), probs = c((1-conf.level)/2,(1+conf.level)/2),na.rm = TRUE)}

table_eff <- par_est[,names(par_est)%in%c("description","estimate_resc","CI_lower_resc","CI_upper_resc")]
names(table_eff) <- c("description","estimate","CI_lower","CI_upper")
table_eff <- rbind(table_eff,data.frame(description=c("EC90 dose","IC50 dose"),
                                        estimate=c(dose_from_eff(0.9*par_est$estimate_resc[par_est$description=="maximal efficacy"]/100,dose_response_par_fit[n_studies+(1:n_par)]),dose_from_eff(0.5,dose_response_par_fit[n_studies+(1:n_par)])),
                                        CI_lower=c(dose_ci_from_eff(0.9*par_est$estimate_resc[par_est$description=="maximal efficacy"]/100)[1],dose_ci_from_eff(0.5)[1]),
                                        CI_upper=c(dose_ci_from_eff(0.9*par_est$estimate_resc[par_est$description=="maximal efficacy"]/100)[2],dose_ci_from_eff(0.5)[1])))

# CIs for EC50 and EC90 via bootstrap but with fixed maximal efficacy (uncertainty in EC50 and in the slope parameter only):
n_par <- 3
n_par2 <- 2
covar2 <- solve(dose_response_fit$hessian)[nrow(dose_response_fit$hessian)-c(n_par2:1),ncol(dose_response_fit$hessian)-c(n_par2:1)] # covariance matrix for efficacy parameters
par_bootstr2 <- rmvnorm(n=N,mean=dose_response_par_fit[length(dose_response_par_fit)-c(n_par2:1)],sigma=covar2)

# compute efficacy CIs for a given dose (with fixed maximal efficacy):
eff_ci_from_dose2 <- function(d){quantile(sapply(c(1:nrow(par_bootstr2)),function(x){eff(d,c(par_bootstr2[x,],tail(dose_response_par_fit,1)))}),
                                          probs = c((1-conf.level)/2,(1+conf.level)/2))}
# compute dose CIs for a given efficacy
dose_ci_from_eff2 <- function(ef){quantile(sapply(c(1:nrow(par_bootstr2)),function(x){dose_from_eff(ef,c(par_bootstr2[x,],tail(dose_response_par_fit,1)))}),
                                           probs = c((1-conf.level)/2,(1+conf.level)/2),na.rm = TRUE)}

table_eff_fixed_m <- data.frame(description=table_eff$description,estimate=table_eff$estimate,
                                CI_lower=c(dose_ci_from_eff2(0.5*par_est$estimate_resc[par_est$description=="maximal efficacy"]/100)[1],table_eff$CI_lower[c(2,3)],
                                           dose_ci_from_eff2(0.9*par_est$estimate_resc[par_est$description=="maximal efficacy"]/100)[1],dose_ci_from_eff2(0.5)[1]),
                                CI_upper=c(dose_ci_from_eff2(0.5*par_est$estimate_resc[par_est$description=="maximal efficacy"]/100)[2],table_eff$CI_upper[c(2,3)],
                                           dose_ci_from_eff2(0.9*par_est$estimate_resc[par_est$description=="maximal efficacy"]/100)[2],dose_ci_from_eff2(0.5)[2]))

# save all estimated parameters:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
save(par_est,covar,table_eff,table_eff_fixed_m, file =glue::glue("output/Dose-resp parameters (symp-hosp)_{Sys.Date()}_{cur_time}.RData"))
write_xlsx(table_eff_fixed_m,glue::glue("output/Dose-resp parameter table (symp-hosp)_{Sys.Date()}_{cur_time}.xlsx")) # Table S11


# visualization of the data and the model fit -----------------------------

# model fit:
dose <- 10^seq(-2,4,length.out = 1e3)
model_fit <- sapply(dose,function(x){100*eff(x,dose_response_par_fit[n_studies+(1:3)])})
model_fit_ci <- sapply(dose,function(x){100*eff_ci_from_dose(x)})
model_data <- data.frame(dose,model_fit,t(model_fit_ci))
names(model_data) <- c("dose","model_fit","lower","upper")

cur_time <- format(Sys.time(),"h%H_m%M_s%S")
save(model_data, file =glue::glue("output/Dose-resp curve (symp-hosp)_{Sys.Date()}_{cur_time}.RData"))

# visualize the data:
dose_resp_tmp <- dose_resp_data
dose_resp_tmp$outcome[dose_resp_tmp$outcome<0] <- 0
dose_resp_tmp$lower[dose_resp_tmp$lower<0] <- 0

dose_resp_aes <- data.frame(Trial=sort(unique(dose_resp_data$Trial)))
dose_resp_tmp <- dose_resp_data
dose_resp_tmp$`Patient risk` <- factor(dose_resp_data$`Patient risk`,levels=c("low","mixed","high"))
my_shapes <- c(0,1,2,4,5,6,7,13,15,16,17,18)
dose_resp_aes <- cbind(dose_resp_aes,shape_ind = rep(my_shapes,ceiling(length(dose_resp_aes$Trial)/length(my_shapes)))[1:length(dose_resp_aes$Trial)])
dose_resp_aes <- dplyr::left_join(dose_resp_aes,dose_resp_tmp[,names(dose_resp_tmp)%in%c("Trial","Patient risk")],by="Trial")

my_colours <- c("yellow2","yellow3","goldenrod1","darkorange","red","red3","orangered4","lightpink","hotpink","orchid3",
                "dodgerblue","turquoise1","palegreen2","aquamarine3","mediumpurple","springgreen4","darkviolet","purple4","black")

Fig2 <- ggplot(data=dose_resp_tmp,aes(x=dose_conv,y=outcome,shape=Trial,colour=Trial,size=`Patient risk`)) +
  # model fit:
  geom_ribbon(data=model_data,inherit.aes=FALSE,aes(x=dose,ymin=lower, ymax=upper),fill="black", alpha = 0.1)+
  geom_line(data=model_data,inherit.aes=FALSE,aes(x=dose,y=model_fit),color="black") +
  # trial data:
  geom_point(size=c(1,2,3)[as.numeric(dose_resp_tmp$`Patient risk`)]) + 
  scale_shape_manual(values=dose_resp_aes$shape_ind) +
  geom_segment(aes(x=dose_conv, y=lower,xend=dose_conv, yend=upper)) +
  scale_color_manual(values=c(my_colours)) + # colour by trial
  scale_size_manual(values=c(0.25,0.6,1.3), name="Patient risk:") + 
  # axis, theme, etc.:
  scale_x_log10(breaks=c(1e-2,1e-1,1e0,1e1,1e2,1e3),labels=c("0.01","0.1","1","10","100","1000")) +
  scale_y_continuous(minor_breaks = seq(0, 100, 20),breaks=seq(0, 100, 20),labels=seq(0, 100, 20),expand = c(1e-2,1e-2)) +
  coord_cartesian(ylim=c(0,100),xlim=c(0.09,1500)) +
  labs(x="Dose [fold convalescent]", y="Efficacy [%]", title = "Therapeutic: symptomatic infection to hospitalisation") +
  theme_bw() +
  theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2))

# Fig2

# save the figure:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
pdf(glue::glue("output/Fig2_{Sys.Date()}_{cur_time}.pdf"),height=4,width=8)
print(Fig2)
dev.off()


# leave-one-out analysis --------------------------------------------------

# leave out one study at a time:
studies_names <- sort(c("Chew","Dougan","Eom","Gottlieb","Gupta","Korley","Libster","Millat-Martinez","Montgomery","Sullivan","Weinreich")) # names of studies to be left out
studies_index <- rep(NA,nrow(dose_resp_data))
for(i in 1:length(studies_names)){
  studies_index[which(grepl(studies_names[i],dose_resp_data$Trial))] <- i
}
par_est_leave <- data.frame(Trial_left_out=studies_names,half_max_dose=rep(NA,length(studies_names)),slope=rep(NA,length(studies_names)),
                            max_eff=rep(NA,length(studies_names)))

# efficacy function by dose (d):
# logistic function with 3 parameters: midpoint (=half-maximal efficacy), growth, and maximal efficacy
eff <- function(d,par){par[3]/(1+exp(-exp(par[2])*(log10(d)-log10(exp(par[1])))))}

# fit data with one study left out:
for(i in 1:length(studies_names)){
  data_tmp <- dose_resp_data[-which(studies_index==i),]
  n_studies_tmp <- nrow(data_tmp)
  nllh <- function(par){-sum(log(dbinom(data_tmp$eCont,data_tmp$nCont,par[1:n_studies_tmp])))-
      sum(log(dbinom(data_tmp$eTreat,data_tmp$nTreat,par[1:n_studies_tmp]*(1-eff(data_tmp$dose_conv,par[n_studies_tmp+(1:3)])))))}
  par_init <- c(data_tmp$eCont/data_tmp$nCont,0,log(1),0.9) 
  fit_tmp <- nlm(nllh,par_init,hessian = TRUE)
  par_tmp <- fit_tmp$estimate[n_studies_tmp+c(1:3)]
  par_est_leave[i,2] <- exp(par_tmp[1]) # half-maximal dose
  par_est_leave[i,3] <- exp(par_tmp[2]) # slope
  par_est_leave[i,4] <- 100*par_tmp[3] # maximal efficacy
}

# Visualize results:
par_all_tmp <- rbind(par_est_leave,data.frame(Trial_left_out="None",half_max_dose=table_eff[1,2],slope=table_eff[2,2],max_eff=table_eff[3,2]))
par_all_tmp <- cbind(par_all_tmp,x=c(rep("Leave-one-out",nrow(par_all_tmp)-1),"All data"))
my_shapes_tmp <- c(0,1,2,4,5,6,8,13,15,16,17,18)

p1 <- ggplot(par_all_tmp,aes(x=x,y=half_max_dose,shape=Trial_left_out)) +
  geom_point(position = position_dodge(width = 1.1)) + 
  scale_shape_manual(values=my_shapes_tmp[1:nrow(par_all_tmp)], name="Trial left out:") +
  geom_segment(inherit.aes=FALSE,aes(x="All data",y=table_eff[1,3],xend="All data",yend=table_eff[1,4])) +
  scale_y_log10(breaks=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6),labels=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6), limits=c(0.1,0.65)) +
  labs(x="", y="Dose [fold conv.]",title="EC50 dose") +
  theme_bw() +
  theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2),
        panel.grid.minor = element_blank())

p2 <- ggplot(par_all_tmp,aes(x=x,y=max_eff,shape=Trial_left_out)) +
  geom_point(position = position_dodge(width = 1.1)) + 
  scale_shape_manual(values=my_shapes_tmp[1:nrow(par_all_tmp)], name="Trial left out:") +
  geom_segment(inherit.aes=FALSE,aes(x="All data",y=table_eff$CI_lower[3],xend="All data",yend=table_eff$CI_upper[3])) + 
  ylim(59,76) + 
  labs(x="", y="Efficacy [%]",title = "Maximal efficacy") +
  theme_bw() +
  theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2),
        panel.grid.minor = element_blank())

p3 <- ggplot(par_all_tmp,aes(x=x,y=slope,shape=Trial_left_out)) +
  geom_point(position = position_dodge(width = 1.1)) + 
  scale_shape_manual(values=my_shapes_tmp[1:nrow(par_all_tmp)], name="Trial left out:") +
  geom_segment(inherit.aes=FALSE,aes(x="All data",y=table_eff$CI_lower[2],xend="All data",yend=table_eff$CI_upper[2])) + 
  scale_y_log10(breaks=c(1,2,4,6,8,10,14,18), limits=c(1,10)) +
  labs(x="", y="Estimated slope",title="Slope parameter") +
  theme_bw() +
  theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2),
        panel.grid.minor = element_blank())

FigS12 <- ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")
# FigS12

# save the figure:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
pdf(glue::glue("output/FigS12_{Sys.Date()}_{cur_time}.pdf"),height=3,width=7)
print(FigS12)
dev.off()


# Efficacy of passive Abs in a non-naive population -----------------------

# efficacy function with variable parameters:
eff <- function(d,par){par[3]/(1+exp(-exp(par[2])*(log10(d)-log10(exp(par[1])))))}
# efficacy function with best fit parameters:
eff_par <- function(d){par_est$estimate[par_est$description=="maximal efficacy"]/(1+exp(-exp(par_est$estimate[par_est$description=="slope"])*(log10(d)-log10(exp(par_est$estimate[par_est$description=="half-maximal dose"])))))}
# function of the initial Ab dose (d):
max_eff_non_naive <- function(d){(par_est$estimate[par_est$description=="maximal efficacy"]-eff_par(d))/(1-eff_par(d))}
# function of the initial Ab dose (d) with an x-fold drop in the neutralisation titre against a variant:
# max_eff_non_naive_var <- function(d,x){(par_est$estimate[par_est$description=="maximal efficacy"]-eff_par(d/x))/(1-eff_par(d/x))}
# function to compute dose for given efficacy (same parameters as eff function, with the same scaling):
dose_from_eff <- function(eff,par){exp(par[1])*10^(-log(par[3]/eff-1)/exp(par[2]))}
# dose by maximal efficacy in non-naive patients:
dose_from_max_eff <- function(max_eff,par){dose_from_eff((par[3]-max_eff)/(1-max_eff),par)}


# Visualize:
dose <- 10^seq(-2,4,length.out = 1e3)
model_fit <- sapply(dose,function(x){100*max_eff_non_naive(x)})

# make 95% CIs using parametric bootstrapping:
N <- 1e5 # number of bootstraps
n_par <- nrow(par_est)
par_bootstr <- rmvnorm(n=N,mean=par_est$estimate,sigma=covar)
# compute efficacy CIs for a given dose:
max_eff_non_naive_ci <- function(d){quantile(sapply(c(1:nrow(par_bootstr)),function(x){(par_bootstr[x,3]-eff(d,par_bootstr[x,]))/(1-eff(d,par_bootstr[x,]))}), probs = c((1-conf.level)/2,(1+conf.level)/2))}

model_fit_ci <- sapply(dose,function(x){100*max_eff_non_naive_ci(x)})

# combine all data and save it:
model_data <- data.frame(dose,model_fit,t(model_fit_ci))#,model_fit_var)
names(model_data) <- c("dose","model_fit","model_fit_lower","model_fit_upper")#),"model_fit_var")
# cur_time <- format(Sys.time(),"h%H_m%M_s%S")
# save(model_data, file =glue::glue("output/Max efficacy non-naive population_{Sys.Date()}_{cur_time}.RData"))

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
model_data_below <- model_data[model_data$model_fit_lower<=hline,]
model_data_below$model_fit_upper[model_data_below$model_fit_upper>=hline] <- hline
model_data_above <- model_data[model_data$model_fit_upper>=hline,]
model_data_above$model_fit_lower[model_data_above$model_fit_lower<=hline] <- hline

Fig4B <- ggplot() +
  # maximal efficacy below the threshold:
  geom_ribbon(data=model_data_below,aes(x=dose,ymin=model_fit_lower, ymax=model_fit_upper),fill="black", alpha = 0.1)+ # TODO
  geom_line(data=model_data[model_data$model_fit<=hline,],aes(x=dose,y=model_fit),color="black", size=1) +
  # maximal efficacy above the threshold:
  geom_ribbon(data=model_data_above,aes(x=dose,ymin=model_fit_lower, ymax=model_fit_upper),fill="dodgerblue", alpha = 0.1)+ # TODO
  geom_line(data=model_data[model_data$model_fit>=hline,],aes(x=dose,y=model_fit),color="dodgerblue", size=1) +
  # add horizonal and vertical lines with annotations:
  geom_segment(aes(x=1e-3,xend=1e4,y=hline,yend=hline),color="red",linetype="dashed") +
  geom_segment(aes(x=dose_from_max_eff(hline/100,par_est$estimate),xend=dose_from_max_eff(hline/100,par_est$estimate),
                   y=hline,yend=0),color="red",linetype="dashed") +
  annotate(geom="text",x=0.1,y=hline,label=paste(hline,"% Efficacy \n",sep=""),color="red") + 
  annotate(geom="text",x=dose_from_max_eff(hline/100,par_est$estimate),y=hline/2,color="red",angle=90,
           label=paste("Neut. for ",hline,"% \n efficacy: ",round(dose_from_max_eff(hline/100,par_est$estimate),2),sep="")) + 
  # add points for different variants:
  geom_point(data=add_neuts_data, aes(x=neuts_conv,y=efficacy,color=variant,shape=variant),size=4) +
  scale_color_manual(values=my_colours) +
  scale_shape_manual(values=my_shapes) +
  # other annotations:
  geom_segment(aes(x=add_neuts_data$neuts_conv[add_neuts_data$variant=="wild type"],xend=add_neuts_data$neuts_conv[add_neuts_data$variant=="wild type"],
                   y=0,yend=100), color="black", linetype="dashed") +
  annotate(geom="text",x=50,y=15,color="black",label="< 30% clinical efficacy of passive \n Ab treatment expected") + 
  annotate(geom="text",x=0.1,y=85,color="black",label="Potential benefit of \n passive Ab treatment \n (if an efficacious Ab \n is available") + 
  # axis, theme, etc.:
  scale_x_log10(breaks=c(1e-2,1e-1,1e0,1e1,1e2,1e3),labels=c("0.01","0.1","1","10","100","1000")) +
  scale_y_continuous(minor_breaks = seq(0, 100, 20),breaks=seq(0, 100, 20),labels=seq(0, 100, 20),expand = c(1e-2,1e-2)) +
  coord_cartesian(ylim=c(0,101),xlim=c(1e-2,1e4),expand=c(0)) +
  labs(x="Endogenous neutralisation level [fold convalescent]", y="Efficacy [%]", 
       title = "Maximal efficacy of passive Abs in a non-naive population") +
  theme_bw() +
  theme(legend.key.size = unit(0.6, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2))
# Fig4B

# save the figure:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
pdf(glue::glue("output/Fig4B_{Sys.Date()}_{cur_time}.pdf"),height=4,width=6)
print(Fig4B)
dev.off()


# cleanup -----------------------------------------------------------------

rm(add_neuts_data,covar,covar2,data_tmp,dose_resp_aes,dose_resp_data,dose_resp_tmp,
   dose_response_fit,dose_response_par_fit_ci,Fig2,Fig4B,FigS12,fit_tmp,model_data,
   model_data_above,model_data_below,model_fit_ci,neut.data.2vax.hist,p1,p2,p3,
   par_all_tmp,par_bootstr,par_bootstr2,par_est,par_est_leave,table_eff,
   table_eff_fixed_m,add_var,cur_time,dose,dose_response_par_fit,hline,i,
   initial_stage,model_fit,my_colours,my_shapes,my_shapes_tmp,N,n_par,n_par2,
   n_studies,n_studies_tmp,neuts,neuts_conv,outcome_stage,outcome_stage_col,
   par_init,par_tmp,studies_index,studies_names,var_names,dose_ci_from_eff,
   dose_ci_from_eff2,dose_from_eff,dose_from_max_eff,eff,eff_ci_from_dose,
   eff_ci_from_dose2,eff_par,max_eff_non_naive,max_eff_non_naive_ci,nllh)
