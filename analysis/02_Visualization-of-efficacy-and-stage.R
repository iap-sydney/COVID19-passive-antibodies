# -------------------------------------------------------------------------
#' 
#' Visualisation of treatment efficacy and stage (Figure 1)
#' 
# -------------------------------------------------------------------------


# Visualization of CP & hIG treatment (Fig. 1A) ---------------------------

# select CP & hIG studies (note that hIG studies are referred to hIVIG in the data)
plas_studies <- studies_overview[studies_overview$`Treatment type`%in%c("CP","hIVIG"),]

# select only relevant rows (no subgroups, except seronegative subgroups):
plas_studies <- rbind(plas_studies[which(is.na(plas_studies$Subgroup)),],plas_studies[which(plas_studies$Subgroup=="Serology (negative)"),])
# rename studies to include subgroup:
plas_studies$Trial[plas_studies$Trial=="RECOVERY (CP)"] <- "RECOVERY"
plas_studies$Trial[which(plas_studies$Subgroup=="Serology (negative)")] <- paste(plas_studies$Trial[which(plas_studies$Subgroup=="Serology (negative)")]," (sero-negative)",sep="")


# Make data frame with stage transitions:
plas_stages <- plas_studies[,names(plas_studies)%in%c("Trial","Treatment stage","symptomatic","hospitalisation","ventilation","death")]
plas_stages <- melt(plas_stages,id=c("Trial","Treatment stage"))
names(plas_stages) <- c("Trial","Treatment stage","Outcome stage","efficacy")
plas_stages <- plas_stages[!is.na(plas_stages$efficacy),]
CI_lower <- sapply(c(1:nrow(plas_stages)),function(x){plas_studies[plas_studies$Trial==plas_stages$Trial[x],which(names(plas_studies)==plas_stages$`Outcome stage`[x])+1]})
CI_upper <- sapply(c(1:nrow(plas_stages)),function(x){plas_studies[plas_studies$Trial==plas_stages$Trial[x],which(names(plas_studies)==plas_stages$`Outcome stage`[x])+2]})
plas_stages <- cbind(plas_stages,"CI lower"=CI_lower,"CI upper"=CI_upper)


### Visualizations:

plas_stages_tmp <- plas_stages
# treatment and outcome stage number (1=pre-exposure, 2=peri-(post-)exposure, 3=symptomatic, 4=hospitalisation, 5=ventilation, 6=death):
treatment_index <- as.numeric(factor(plas_stages_tmp$`Treatment stage`,levels = c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death")))
outcome_index <- as.numeric(factor(plas_stages_tmp$`Outcome stage`,levels = c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death")))
plas_stages_tmp <- cbind(plas_stages_tmp,treatment_index,outcome_index)

# "dodge" at treatment stage:
plas_stages_tmp <- plas_stages_tmp[order(plas_stages_tmp$efficacy,decreasing = TRUE),] # re-order by efficacy
for(i in 1:length(unique(plas_stages_tmp$treatment_index))){
  tmp <- which(plas_stages_tmp$treatment_index==unique(treatment_index)[i])
  d <- 0.027 # dodge by fixed distance d
  plas_stages_tmp$treatment_index[tmp] <- plas_stages_tmp$treatment_index[tmp[1]]+d*(0:(length(tmp)-1))-(length(tmp)-1)*d/2 
}

# include a dummy variable for linetype:
tmp <- rep("solid",length(plas_stages_tmp$Trial))
tmp[plas_stages_tmp$Trial=="RECOVERY (sero-negative)" | plas_stages_tmp$Trial=="REMAP-CAP (sero-negative)"] <- "dashed"
plas_stages_tmp <- cbind(plas_stages_tmp,line_type=as.factor(tmp))

# re-scale negative values (divide by 5):
plas_stages_tmp$efficacy[plas_stages_tmp$efficacy<0] <- plas_stages_tmp$efficacy[plas_stages_tmp$efficacy<0]/5
plas_stages_tmp$`CI lower`[plas_stages_tmp$`CI lower`<0] <- plas_stages_tmp$`CI lower`[plas_stages_tmp$`CI lower`<0]/5
plas_stages_tmp$`CI upper`[plas_stages_tmp$`CI upper`<0] <- plas_stages_tmp$`CI upper`[plas_stages_tmp$`CI upper`<0]/5

plas_stages_tmp$Trial <- as.factor(plas_stages_tmp$Trial)

# colour, shape, and linetype by Trial:
tmp_trial <- sort(unique(plas_stages_tmp$Trial))
plas_aes <- data.frame(Trial=tmp_trial)
plas_aes <- dplyr::left_join(plas_aes,plas_stages_tmp[,which(names(plas_stages_tmp)%in%c("Trial","line_type"))],by="Trial")
plas_aes <- unique(plas_aes)
plas_aes <- cbind(plas_aes,colour_ind = c(1:31,31,32,32:40))
my_shapes <- c(0,1,2,4,5,6,7,13,15,16,17,18)
plas_aes <- cbind(plas_aes,shape_ind = rep(my_shapes,ceiling(length(plas_aes$Trial)/length(my_shapes)))[1:length(plas_aes$Trial)])

my_colours <- c("#000000",scales::seq_gradient_pal("#000000","#E69F00","Lab")(0.5),"#E69F00",
                scales::seq_gradient_pal("#E69F00","#56B4E9","Lab")(0.5),"#56B4E9",
                scales::seq_gradient_pal("#56B4E9","#009E73","Lab")(0.5),"#009E73",
                scales::seq_gradient_pal("#009E73","#F0E442","Lab")(0.5),"#F0E442",
                scales::seq_gradient_pal("#F0E442","#0072B2","Lab")(0.2),"#0072B2",
                scales::seq_gradient_pal("#0072B2","#D55E00","Lab")(0.5),"#D55E00",
                scales::seq_gradient_pal("#D55E00","#CC79A7","Lab")(0.5),"#CC79A7","red")

if(max(plas_aes$colour_ind)>length(my_colours)){
  my_colours <- rep(my_colours,ceiling(max(plas_aes$colour_ind)/length(my_colours)))
}

Fig1A <- ggplot(data=plas_stages_tmp,aes(x=treatment_index,y=efficacy,color=Trial,shape=Trial,linetype=Trial,alpha=factor(`CI lower`>0))) +
  geom_rect(fill = "#999999",xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = 0, alpha = 0.006, color =NA) + # grey shaded area for negative efficacy
  geom_hline(yintercept=0, color="#999999") + # bolder line at 0
  # geom_hline(yintercept=-24, color="#999999") + # bolder line at "< -100"
  geom_segment(aes(x = treatment_index, y = `CI lower`, xend = treatment_index, yend = `CI upper`), size=0.2) + # vertical line for CIs
  geom_segment(aes(x = treatment_index, y = efficacy, xend = outcome_index, yend = efficacy), size=0.7) + # horizontal line between stages
  scale_linetype_manual(values=as.character(plas_aes$line_type)) +
  geom_point(size=1.5) + # initial stage points
  geom_point(aes(x=outcome_index,y=efficacy),size=0.9) + # outcome stage points
  scale_shape_manual(values=plas_aes$shape_ind) +
  scale_color_manual(values=my_colours[plas_aes$colour_ind]) +
  scale_alpha_manual(values=c(0.17,1),name=NULL,labels=c("not significant","significant")) +
  labs(x="Stage", y="Efficacy [%]") +
  scale_y_continuous(minor_breaks = c(seq(-24,-4,4),seq(0, 100, 20)), breaks=c(seq(-24,-4,4),seq(0, 100, 20)),
                     labels=c("below -100","","-80","","-40","",seq(0, 100, 20)),expand=c(0,1.5)) +
  scale_x_discrete(breaks=c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death"),
                   limits = c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death"),
                   labels=c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death")) +
  coord_cartesian(ylim=c(-20,100)) +
  guides(size = "none", color = guide_legend(ncol=1,title.position = "top"), shape = guide_legend(ncol=1),
         linetype  = guide_legend(ncol=1,title.position = "top")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2))

# Fig1A

cur_time <- format(Sys.time(),"h%H_m%M_s%S")
pdf(glue::glue("output/Fig1A_{Sys.Date()}_{cur_time}.pdf"),height=4,width=8)
print(Fig1A)
dev.off()

# cleanup
rm(plas_aes,plas_stages,plas_stages_tmp,plas_studies,CI_upper,CI_lower,d,i,my_colours,
   my_shapes,outcome_index,tmp,tmp_trial,treatment_index,cur_time)


# Visualization of Ab treatment (Fig. 1B) ---------------------------------

# select Ab (antibody) studies:
Ab_studies <- studies_overview[studies_overview$`Treatment type`=="mAb",]

# select only relevant rows (no subgroups, except seronegative subgroups and O'Brien by time):
Ab_studies <- rbind(Ab_studies[which(is.na(Ab_studies$Subgroup)),],
                    Ab_studies[which(Ab_studies$Subgroup=="Serology (negative)"),],
                    Ab_studies[which(Ab_studies$Trial=="O'Brien (PCR-)" & grepl("Time",Ab_studies$Subgroup)),])
# rename studies to include mAb and/or subgroup:
Ab_studies$Trial[which(Ab_studies$Subgroup=="Serology (negative)")] <- paste(Ab_studies$Trial[which(Ab_studies$Subgroup=="Serology (negative)")]," (sero-negative)",sep="")
Ab_studies$Trial[which(Ab_studies$Trial=="O'Brien (PCR-)" & grepl("Time",Ab_studies$Subgroup))] <- paste("O'Brien (PCR-)",gsub("Time of symptomatic infection "," ",Ab_studies$Subgroup[which(Ab_studies$Trial=="O'Brien (PCR-)" & grepl("Time",Ab_studies$Subgroup))]),sep="")
Ab_studies$Trial[which(is.na(Ab_studies$Subgroup))] <- paste(Ab_studies$Trial[which(is.na(Ab_studies$Subgroup))]," (",Ab_studies$Treatment[which(is.na(Ab_studies$Subgroup))],")",sep="")
Ab_studies$Trial[Ab_studies$Trial=="ACTIV-3 (evu) (sero-negative)"] <- "ACTIV-3 (cilgavimab + tixagevimab, sero-negative)"
char_remove <- c(" \\(mAb\\) "," \\(amu\\+rom\\) "," \\(bam\\) "," \\(bam \\+ ete\\) "," \\(bam\\+ete, 2.1g\\) "," \\(evu\\) "," \\(sot\\) ")
for(i in 1:length(char_remove)){
  Ab_studies$Trial[which(grepl(char_remove[i],Ab_studies$Trial))] <- gsub(char_remove[i]," ",Ab_studies$Trial[which(grepl(char_remove[i],Ab_studies$Trial))])
}
Ab_studies$Trial[which(grepl("\\) \\(",Ab_studies$Trial))] <- gsub("\\) \\(",", ",Ab_studies$Trial[which(grepl("\\) \\(",Ab_studies$Trial))])


# Make data frame with stage transitions:
Ab_stages <- Ab_studies[,names(Ab_studies)%in%c("Trial","Treatment stage","symptomatic","hospitalisation","ventilation","death")]
Ab_stages <- melt(Ab_stages,id=c("Trial","Treatment stage"))
names(Ab_stages) <- c("Trial","Treatment stage","Outcome stage","efficacy")
Ab_stages <- Ab_stages[!is.na(Ab_stages$efficacy),]
CI_lower <- sapply(c(1:nrow(Ab_stages)),function(x){Ab_studies[Ab_studies$Trial==Ab_stages$Trial[x],which(names(Ab_studies)==Ab_stages$`Outcome stage`[x])+1]})
CI_upper <- sapply(c(1:nrow(Ab_stages)),function(x){Ab_studies[Ab_studies$Trial==Ab_stages$Trial[x],which(names(Ab_studies)==Ab_stages$`Outcome stage`[x])+2]})
Ab_stages <- cbind(Ab_stages,"CI lower"=CI_lower,"CI upper"=CI_upper)


### Visualizations:

Ab_stages_tmp <- Ab_stages
# treatment and outcome stage number (1=pre-exposure, 2=peri-(post-)exposure, 3=symptomatic, 4=hospitalisation, 5=ventilation, 6=death):
treatment_index <- as.numeric(factor(Ab_stages_tmp$`Treatment stage`,levels = c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death")))
outcome_index <- as.numeric(factor(Ab_stages_tmp$`Outcome stage`,levels = c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death")))
Ab_stages_tmp <- cbind(Ab_stages_tmp,treatment_index,outcome_index)

# "dodge" at treatment stage:
Ab_stages_tmp <- Ab_stages_tmp[order(Ab_stages_tmp$efficacy,decreasing = TRUE),] # re-order by efficacy
for(i in 1:length(unique(Ab_stages_tmp$treatment_index))){
  tmp <- which(Ab_stages_tmp$treatment_index==unique(treatment_index)[i])
  d <- 0.05 # dodge by fixed distance d
  Ab_stages_tmp$treatment_index[tmp] <- Ab_stages_tmp$treatment_index[tmp[1]]+d*(0:(length(tmp)-1))-(length(tmp)-1)*d/2 
}

# include a dummy variable for linetype:
tmp <- rep("solid",length(Ab_stages_tmp$Trial))
tmp[grepl("sero-negative",Ab_stages_tmp$Trial) | grepl("week",Ab_stages_tmp$Trial)] <- "dashed"
Ab_stages_tmp <- cbind(Ab_stages_tmp,line_type=as.factor(tmp))

# re-scale negative values (divide by 5):
Ab_stages_tmp$efficacy[Ab_stages_tmp$efficacy<0] <- Ab_stages_tmp$efficacy[Ab_stages_tmp$efficacy<0]/5
Ab_stages_tmp$`CI lower`[Ab_stages_tmp$`CI lower`<0] <- Ab_stages_tmp$`CI lower`[Ab_stages_tmp$`CI lower`<0]/5
Ab_stages_tmp$`CI upper`[Ab_stages_tmp$`CI upper`<0] <- Ab_stages_tmp$`CI upper`[Ab_stages_tmp$`CI upper`<0]/5

Ab_stages_tmp$Trial <- as.factor(Ab_stages_tmp$Trial)

# colour, shape, and linetype by Trial:
tmp_trial <- sort(unique(Ab_stages_tmp$Trial))
Ab_aes <- data.frame(Trial=tmp_trial)
Ab_aes <- dplyr::left_join(Ab_aes,Ab_stages_tmp[,which(names(Ab_stages_tmp)%in%c("Trial","line_type"))],by="Trial")
Ab_aes <- unique(Ab_aes)
Ab_aes <- cbind(Ab_aes,colour_ind = c(1:3,3:15,15,rep(16,4),rep(17,2),rep(18,2),rep(19,2),rep(20,2)))
my_shapes <- c(0,1,2,4,5,6,13,15,16,17,18)
Ab_aes <- cbind(Ab_aes,shape_ind = rep(my_shapes,ceiling(length(Ab_aes$Trial)/length(my_shapes)))[1:length(Ab_aes$Trial)])

my_colours <- c(scales::seq_gradient_pal("#000000","#E69F00","Lab")(0.5),"#E69F00",
                scales::seq_gradient_pal("#E69F00","#56B4E9","Lab")(0.5),"#0072B2","#56B4E9","aquamarine3","forestgreen",
                scales::seq_gradient_pal("#009E73","#F0E442","Lab")(0.5),"#F0E442",
                scales::seq_gradient_pal("#0072B2","#D55E00","Lab")(0.5),"#D55E00",
                scales::seq_gradient_pal("#D55E00","#CC79A7","Lab")(0.5),"#CC79A7","red")

if(max(Ab_aes$colour_ind)>length(my_colours)){
  my_colours <- rep(my_colours,ceiling(max(Ab_aes$colour_ind)/length(my_colours)))
}

Fig1B <- ggplot(data=Ab_stages_tmp,aes(x=treatment_index,y=efficacy,color=Trial,shape=Trial,linetype=Trial,alpha=factor(`CI lower`>0))) +
  geom_rect(fill = "#999999",xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = 0, alpha = 0.006, color =NA) + # grey shaded area for negative efficacy
  geom_hline(yintercept=0, color="#999999") + # bolder line at 0
  geom_segment(aes(x = treatment_index, y = `CI lower`, xend = treatment_index, yend = `CI upper`), size=0.2) + # vertical line for CIs
  geom_segment(aes(x = treatment_index, y = efficacy, xend = outcome_index, yend = efficacy), size=0.7) + # horizontal line between stages
  scale_linetype_manual(values=as.character(Ab_aes$line_type)) +
  geom_point(size=1.5) + # initial stage points
  geom_point(aes(x=outcome_index,y=efficacy),size=0.9) + # outcome stage points
  scale_shape_manual(values=Ab_aes$shape_ind) +
  scale_color_manual(values=my_colours[Ab_aes$colour_ind]) +
  scale_alpha_manual(values=c(0.17,1),name=NULL,labels=c("not significant","significant")) +
  labs(x="Stage", y="Efficacy [%]") +
  scale_y_continuous(minor_breaks = c(seq(-24,-4,4),seq(0, 100, 20)), breaks=c(seq(-24,-4,4),seq(0, 100, 20)),
                     labels=c("below -100","","-80","","-40","",seq(0, 100, 20)),expand=c(0,1.5)) +
  scale_x_discrete(breaks=c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death"),
                   limits = c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death"),
                   labels=c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death")) +
  coord_cartesian(ylim=c(-20,100)) +
  guides(size = "none", color = guide_legend(ncol=1,title.position = "top"), shape = guide_legend(ncol=1),
         linetype  = guide_legend(ncol=1,title.position = "top")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2))

# Fig1B

cur_time <- format(Sys.time(),"h%H_m%M_s%S")
pdf(glue::glue("output/Fig1B_{Sys.Date()}_{cur_time}.pdf"),height=4,width=8)
print(Fig1B)
dev.off()

# cleanup:
rm(Ab_aes,Ab_stages,Ab_stages_tmp,Ab_studies,char_remove,CI_lower,CI_upper,d,i,my_colours,
   my_shapes,outcome_index,tmp,tmp_trial,treatment_index,cur_time)


# Visualization of pooled data (Fig. 1C) ----------------------------------

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
pooled_data <- cbind(pooled_data,"CI lower"=NA,"CI upper"=NA)
pooled_data <- rbind(pooled_data,pooled_data)
pooled_data <- cbind(pooled_data,"treatment"=c(rep("mAb",nrow(pooled_data)/2),rep("plas",nrow(pooled_data)/2)))
pooled_data <- pooled_data[-which(pooled_data$initial_stage==pooled_data$outcome_stage),]

# combine studies and calculate the outcome:
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
        data_tmp <- data.frame(Trial=rep(data_tmp$Trial,2),"events"=c(data_tmp[,2],data_tmp[,4]),
                               "n"=c(data_tmp[,3],data_tmp[,5]),"treatment"=c(rep("treatment",nrow(data_tmp)),rep("control",nrow(data_tmp))))
        data_tmp <- dplyr::mutate(data_tmp,no_events=n-events)
        tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment + (1 | Trial),family=binomial("log"), data=data_tmp)
        pooled_data$efficacy[i] <- exp(fixef(tmp.glmer)[2]) # estimated relative risk for pooled data
        tmp.ci <- exp(confint(tmp.glmer)) # confidence interval for relative risk
        pooled_data$`CI lower`[i] <- tmp.ci[3,1] 
        pooled_data$`CI upper`[i] <- tmp.ci[3,2] 
        # drop1(tmp.glmer, test="Chisq") # test for significance
        # convert to efficacy:
        pooled_data$efficacy[i] <- 100*(1-pooled_data$efficacy[i])
        tmp <- pooled_data$`CI lower`[i]
        pooled_data$`CI lower`[i] <- 100*(1-pooled_data$`CI upper`[i])
        pooled_data$`CI upper`[i] <- 100*(1-tmp)
      }else{
        pooled_data$efficacy[i] <- studies_tmp[studies_tmp$Trial==data_tmp$Trial,which(names(studies_tmp)==pooled_data$outcome_stage[i])]
        pooled_data$`CI lower`[i] <- studies_tmp[studies_tmp$Trial==data_tmp$Trial,which(names(studies_tmp)==pooled_data$outcome_stage[i])+1]
        pooled_data$`CI upper`[i] <- studies_tmp[studies_tmp$Trial==data_tmp$Trial,which(names(studies_tmp)==pooled_data$outcome_stage[i])+2]
      }
    }
  }
}

pooled_data <- pooled_data[!is.na(pooled_data$efficacy),]

# save efficacy estimates for for pooled data: (Table S3)
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(pooled_data,glue::glue("output/Efficacy-pooled-data_{Sys.Date()}_{cur_time}.xlsx"))

### Visualization:

pooled_data_tmp <- pooled_data
# include stage index:
pooled_data_tmp <- dplyr::mutate(pooled_data_tmp,initial_stage_index=as.numeric(initial_stage),outcome_stage_index=as.numeric(outcome_stage))
# "dodge" at initial stage:
pooled_data_tmp <- pooled_data_tmp[order(pooled_data_tmp$efficacy,decreasing = TRUE),]
for(i in 1:length(unique(pooled_data_tmp$initial_stage_index))){
  tmp <- which(pooled_data_tmp$initial_stage_index==unique(as.numeric(pooled_data_tmp$initial_stage))[i])
  d <- 0.06 # dodge by fixed distance d
  pooled_data_tmp$initial_stage_index[tmp] <- pooled_data_tmp$initial_stage_index[tmp[1]]+d*(0:(length(tmp)-1))-(length(tmp)-1)*d/2 
}

# re-scale negative values (divide by 5):
pooled_data_tmp$efficacy[pooled_data_tmp$efficacy<0] <- pooled_data_tmp$efficacy[pooled_data_tmp$efficacy<0]/5
pooled_data_tmp$`CI lower`[pooled_data_tmp$`CI lower`<0] <- pooled_data_tmp$`CI lower`[pooled_data_tmp$`CI lower`<0]/5
pooled_data_tmp$`CI upper`[pooled_data_tmp$`CI upper`<0] <- pooled_data_tmp$`CI upper`[pooled_data_tmp$`CI upper`<0]/5


Fig1C <- ggplot(data=pooled_data_tmp,aes(x=initial_stage_index,y=efficacy,color=treatment,shape=treatment,alpha=factor(`CI lower`>0))) +
  geom_rect(fill = "#999999",xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = 0, alpha = 0.006, color =NA) + # grey shaded area for negative efficacy
  geom_hline(yintercept=0, color="#999999") + # bolder line at 0
  # geom_hline(yintercept=-24, color="#999999") + # bolder line at "< -100"
  geom_segment(aes(x = initial_stage_index, y = `CI lower`, xend = initial_stage_index, yend = `CI upper`), size=0.2) + # vertical line for CIs
  geom_segment(aes(x = initial_stage_index, y = efficacy, xend = outcome_stage_index, yend = efficacy), size=0.7) + # horizontal line between stages
  geom_point(size=1.5) + # initial stage points
  geom_point(aes(x=outcome_stage_index,y=efficacy),size=0.9) + # outcome stage points
  scale_shape_manual(values=c(16,15),name="Treatment",labels=c("mAb","CP and hIVIG")) +
  scale_color_manual(values=c("darkorange2","dodgerblue"),name="Treatment",labels=c("mAb","CP and hIVIG")) +
  scale_alpha_manual(values=c(0.17,1),name=NULL,labels=c("not significant","significant")) +
  labs(x="Stage", y="Efficacy [%]") +
  scale_y_continuous(minor_breaks = c(seq(-24,-4,4),seq(0, 100, 20)), breaks=c(seq(-24,-4,4),seq(0, 100, 20)),
                     labels=c("below -100","","-80","","-40","",seq(0, 100, 20)),expand=c(0,1.5)) +
  scale_x_discrete(breaks=c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death"),
                   limits = c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death"),
                   labels=c("pre-exposure","peri-(post-)exposure","symptomatic","hospitalisation","ventilation","death")) + 
  coord_cartesian(ylim=c(-20,100)) + 
  guides(size = "none", color = guide_legend(ncol=1,title.position = "top"), shape = guide_legend(ncol=1), 
         linetype  = guide_legend(ncol=1,title.position = "top")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2))

Fig1C

cur_time <- format(Sys.time(),"h%H_m%M_s%S")
pdf(glue::glue("output/Fig1C_{Sys.Date()}_{cur_time}.pdf"),height=4,width=8)
print(Fig1C)
dev.off()

# cleanup:
rm(data_tmp,pooled_data,pooled_data_tmp,studies_tmp,tmp.ci,tmp.glmer,d,i,outcome_stage_col,stages,
   tmp,cur_time)


# Combine subfigures and save the figure ----------------------------------

# Combine Fig1A, Fig1B, and Fig1C:
Fig1 <- Fig1A / Fig1B / Fig1C
# Fig1

cur_time <- format(Sys.time(),"h%H_m%M_s%S")
pdf(glue::glue("output/Fig1_{Sys.Date()}_{cur_time}.pdf"),height=12,width=8)
print(Fig1)
dev.off()


# cleanup -----------------------------------------------------------------

rm(Fig1,Fig1A,Fig1B,Fig1C,cur_time)
