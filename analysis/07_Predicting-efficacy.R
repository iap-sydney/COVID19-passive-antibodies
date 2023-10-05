# -------------------------------------------------------------------------
#' 
#' Predicting the efficacy of mAbs against different SARS-CoV-2 variants
#' 
# -------------------------------------------------------------------------


TreatmentRegimen<-studies_overview[studies_overview$`Treatment type`=="mAb",c("Treatment","dose_mg")]
TreatmentRegimen<-TreatmentRegimen[!is.na(TreatmentRegimen$dose_mg),]
TreatmentRegimen<-unique(TreatmentRegimen)


TableOfEff=Means_table[,c(1,2)]

TableOfEff<-TableOfEff[TableOfEff$Antibodies %in% TreatmentRegimen$Treatment,]

AbList<-unique(TreatmentRegimen$Treatment)
for (i in 1:length(AbList)){
  
  tempnumberdoses<-sum(TreatmentRegimen$Treatment== AbList[i])
  tempdoselist<-TreatmentRegimen$dose_mg[TreatmentRegimen$Treatment== AbList[i]]
  
  
  preTableOfEff_dose<-TableOfEff[TableOfEff$Antibodies== AbList[i],]
  
  for (j in 1:tempnumberdoses){
    tempTableOfEff_dose<-preTableOfEff_dose
    tempTableOfEff_dose$dose_mg<-tempdoselist[j]
    if (i==1 & j==1){
      
      
      TableOfEff_dose=tempTableOfEff_dose
    } else {
      TableOfEff_dose<-rbind(TableOfEff_dose,tempTableOfEff_dose)  
    }
  }
}

# load estimated parameters for the dose-response curve:
load("raw-data/Dose-resp parameters (symp-hosp).RData")

# load meta-analysis fit data:
load("raw-data/metafit_20230220_h22_m34_s06.RData")

eff <- function(d,par){par[3]/(1+exp(-exp(par[2])*(log10(d)-log10(exp(par[1])))))}



###Bootstrap
N=10000
ModelParam<-rmvnorm(N,mean=par_est$estimate,sigma=covar)
IC50_log_rand<- (-rmvnorm(N,mean=fit3$beta,sigma=fit3$varFix))
IC50_rand<-exp(IC50_log_rand)
ConvSamples=1/IC50_rand[,which(Means_table$Antibodies=="convalescent")]


TableOfEff_dose$IC50<-NA
TableOfEff_dose$dose_conv<-NA
TableOfEff_dose$efficacy<-NA
TableOfEff_dose$efficacy_LB<-NA
TableOfEff_dose$efficacy_UB<-NA
for (i in 1:nrow(TableOfEff_dose)) {
  
  best_est_IC50<-Means_table$IC50[which(tolower(Means_table$Group)==tolower(paste(TableOfEff_dose$Antibodies[i],TableOfEff_dose$VariantCategory[i])))]
  
  
  
  if (best_est_IC50<10000) {
    best_est_IC50_conv=best_est_IC50*conv
    best_est_dose_conv <- (TableOfEff_dose$dose_mg[i]/plasma_L)/(best_est_IC50_conv*1e-3)
    tempIC50_conv=IC50_rand[,which(tolower(Means_table$Group)==tolower(paste(TableOfEff_dose$Antibodies[i],TableOfEff_dose$VariantCategory[i])))]*ConvSamples
    temp_dose_conv <- (TableOfEff_dose$dose_mg[i]/plasma_L)/(tempIC50_conv*1e-3)
    temp_eff<-NULL
    for (j in 1:N){
      temp_eff[j]<-eff(temp_dose_conv[j],ModelParam[j,])
    }
    best_est_eff<-eff(best_est_dose_conv,par_est$estimate)
    TableOfEff_dose$efficacy_LB[i]<-100*quantile(temp_eff,0.025)
    TableOfEff_dose$efficacy_UB[i]<-100*quantile(temp_eff,0.975)
  } else {
    best_est_IC50_conv=10000*conv
    best_est_dose_conv <- (TableOfEff_dose$dose_mg[i]/plasma_L)/(best_est_IC50_conv*1e-3)
    tempIC50_conv=10000*ConvSamples
    temp_dose_conv <- (TableOfEff_dose$dose_mg[i]/plasma_L)/(tempIC50_conv*1e-3)
    temp_eff<-NULL
    for (j in 1:N){
      temp_eff[j]<-eff(temp_dose_conv[j],ModelParam[j,])
    }
    best_est_eff<-eff(best_est_dose_conv,par_est$estimate)
    TableOfEff_dose$efficacy_LB[i]<-100*0
    TableOfEff_dose$efficacy_UB[i]<-100*quantile(temp_eff,0.95)
    
  }
  
  
  
  TableOfEff_dose$IC50[i]<-best_est_IC50
  TableOfEff_dose$censIC50[i]<-(best_est_IC50>=10000)
  TableOfEff_dose$dose_conv[i]<-best_est_dose_conv
  TableOfEff_dose$efficacy[i]<-100*best_est_eff
  
}

TableOfEff_dose$censIC50<-(TableOfEff_dose$IC50>=10000)

PredictedEfficacyFig<-ggplot(TableOfEff_dose[grepl("Omicron",TableOfEff_dose$VariantCategory),],aes(x=Antibodies,y=efficacy,color=log2(dose_mg),group=dose_mg,alpha=(!censIC50)))+
  geom_point(position=position_dodge(width=0.5),shape=21,fill=NA)+
  geom_errorbar(aes(ymin=efficacy_LB,ymax=efficacy_UB),width=0.3,position=position_dodge(width=0.5))+
  facet_wrap(~VariantCategory,ncol=3)+
  scale_color_continuous(breaks=c(log2(500),log2(1000),log2(2000),log2(4000),log2(8000)),labels=c(500,1000,2000,4000,8000))+
  scale_y_continuous(limits=c(0,100),breaks=c(0,20,40,60,80,100))+
  scale_alpha_manual(values=c(0.25,1))+
  geom_hline(yintercept=30,linetype=3)+
  theme_classic()+
  labs(y="Efficacy (%)",
       alpha="In vitro\nneutralisation\ndetected",
       color="Dose (mg)")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill=NA,color=NA),
        plot.margin = margin("l"=1.5,"unit"="cm"))

# PredictedEfficacyFig

# save Fig. 3:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
ggsave(PredictedEfficacyFig,file=glue::glue("output/Fig3_PredictedEfficacy__{Sys.Date()}_{cur_time}.pdf"),width=7,height=9)


# cleanup -----------------------------------------------------------------

rm(covar,fit3,IC50_log_rand,IC50_rand,ModelParam,par_est,PredictedEfficacyFig,
   preTableOfEff_dose,table_eff,table_eff_fixed_m,TableOfEff,TableOfEff_dose,
   tempTableOfEff_dose,TreatmentRegimen,AbList,best_est_dose_conv,best_est_eff,
   best_est_IC50,best_est_IC50_conv,ConvSamples,cur_time,i,j,N,temp_dose_conv,
   temp_eff,tempdoselist,tempIC50_conv,tempnumberdoses,eff)
