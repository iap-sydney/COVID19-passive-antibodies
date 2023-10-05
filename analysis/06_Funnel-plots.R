# -------------------------------------------------------------------------
#' 
#' Funnel plots to visually assess evidence for publication bias
#' 
# -------------------------------------------------------------------------


# select rows such that no studies & data are included multiple times (by removing subgroup analysis rows)
# for O'Brien (PCR-) study: use data grouped by time to infection (1 week vs 2-4 weeks)
studies_tmp <- rbind(studies_overview[which(is.na(studies_overview$Subgroup) & studies_overview$Trial!="O'Brien (PCR-)"),],
                     studies_overview[which(studies_overview$Trial=="O'Brien (PCR-)" & grepl("Time",studies_overview$Subgroup)),])

# make new data for pooled studies (pool same initial and outcome stage for mAb treatment and for CP/hIVIG treatment)
stages <- factor(c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death"),
                 levels=c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death"))
pooled_data <- data.frame(initial_stage = stages,"pre-exposure"=NA,"peri-(post-)exposure"=NA,"symptomatic"=NA,"hospitalisation"=NA,
                          "ventilation"=NA,"death"=NA)
pooled_data <- melt(pooled_data,id="initial_stage")
names(pooled_data) <- c("initial_stage","outcome_stage","efficacy")
pooled_data$outcome_stage <- stages[as.numeric(pooled_data$outcome_stage)]
pooled_data <- cbind(pooled_data,"number"=NA,"CI lower"=NA,"CI upper"=NA)
pooled_data <- rbind(pooled_data,pooled_data)
pooled_data <- cbind(pooled_data,"treatment"=c(rep("mAb",nrow(pooled_data)/2),rep("plas",nrow(pooled_data)/2)))
pooled_data <- pooled_data[-which(pooled_data$initial_stage==pooled_data$outcome_stage),]


outcome_stage_col <- sapply(c(1:length(stages)),function(x){ifelse(any(names(studies_tmp)==stages[x]),which(names(studies_tmp)==stages[x]),NA)})
for(i in 1:nrow(pooled_data)){
  data_tmp <- studies_tmp[studies_tmp$`Treatment stage`==pooled_data$initial_stage[i],]
  if(pooled_data$treatment[i]=="mAb"){
    data_tmp <- data_tmp[data_tmp$`Treatment type`=="mAb",]
  }else if(pooled_data$treatment[i]=="plas"){
    data_tmp <- data_tmp[data_tmp$`Treatment type`%in%c("CP","hIVIG"),]
  }
  if(!is.na(outcome_stage_col[as.numeric(pooled_data$outcome_stage[i])]) & nrow(data_tmp)>0){
    data_tmp <- data_tmp[,c(1,outcome_stage_col[as.numeric(pooled_data$outcome_stage[i])]+(3:6))]
    # remove NA rows and proceed only if there are at least two studies
    data_tmp <- data_tmp[!is.na(data_tmp[,2]),]
    if(sum(is.na(data_tmp[,2]))<nrow(data_tmp)){
      if(nrow(data_tmp)>1){
        res<-rma(measure="RR", ai=data_tmp[,2], bi=data_tmp[,3]-data_tmp[,2], ci=data_tmp[,4], di=data_tmp[,5]-data_tmp[,4])
        cur_time <- format(Sys.time(),"h%H_m%M_s%S")
        
        # save funnel plots (Fig. S8):
        pdf(file=glue::glue("output/Funnel_{pooled_data$initial_stage[i]}_{pooled_data$outcome_stage[i]}_{pooled_data$treatment[i]}_{Sys.Date()}_{cur_time}.pdf"),width=4.5,height=4)
        funnel(res, yaxis="seinv", main=paste0(ifelse(pooled_data$treatment[i]=="plas","CP/hIVIG",pooled_data$treatment[i]),"\n",pooled_data$initial_stage[i]," to\nprevent ", pooled_data$outcome_stage[i]))
        dev.off()

      }
      
    }
  }
}


# cleanup -----------------------------------------------------------------

rm(data_tmp,pooled_data,res,studies_tmp,cur_time,i,outcome_stage_col,stages)
