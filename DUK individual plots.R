
# Use sex, ethnicity, HbA1c, age and BMI from individuals
# Assume no CKD or renal disease as for primary prevention (also assume no AF - <10% with)
# Use median duration and tds_2011 and look at what popular smoking categories are

########################################################################################################################

library(tidyverse)
library(EHRBiomarkr)


# Import rest of cohort to get median characteristics
## Calibrate to remove calibration subset

# 1 Cohort selection and variable setup

## A Cohort selection (see cohort_definition function for details)

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Raw data/")
load("20230308_t2d_1stinstance.Rda")
load("20230308_t2d_all_drug_periods.Rda")

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Scripts/Functions")
source("cohort_definition.R")

cohort <- define_cohort(t2d_1stinstance, t2d_all_drug_periods)

cohort <- cohort %>% filter(studydrug!="GLP1")

table(cohort$studydrug)
# DPP4SU 95362
# SGLT2 50234


## B Make variables for survival analysis of all endpoints (see survival_variables function for details)

source("survival_variables.R")

cohort <- add_surv_vars(cohort, main_only=TRUE)


## C Just keep variables of interest

cohort <- cohort %>%
  
  select(patid, malesex, ethnicity_5cat_decoded, imd2015_10, tds_2011, regstartdate, gp_record_end, death_date, drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, prehba1c, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preegfr, preckdstage, qrisk2_smoking_cat, contains("cens"), qrisk2_lin_predictor, qrisk2_5yr_score, qdiabeteshf_lin_predictor, qdiabeteshf_5yr_score, starts_with("ckdpc"), last_sglt2_stop)

rm(list=setdiff(ls(), "cohort"))


## D Calibrate score

full_cohort <- cohort

source("calibrate_risk_scores.R")

cohort <- calibrate_risk_score(cohort, risk_score="qdiabeteshf", outcome="hf")

table(cohort$studydrug)
# DPP4SU 95362
# SGLT2 50234

########################################################################################################################

# Get median 

table(cohort$qrisk2_smoking_cat) #use non-smoking
summary(cohort$dstartdate_dm_dur_all) #median=6.6 years
summary(cohort$tds_2011) #median=-1.485

individuals <- data.frame(name=c("Anthony", "Beryl", "Claire"), sex=c("male", "female", "female"), age=c(75, 65, 53), ethrisk=c(0, 0, 5), smoking=c(0, 0, 0), duration_cat=c(2, 4, 0), type1=c(0, 0, 0), cvd=c(0, 0, 0), renal=c(0, 0, 0), af=c(0, 0, 0), hba1c=c(62, 75, 78), cholhdl=c(NA, 11, NA), sbp=c(NA, 190, NA), bmi=c(27, 40, 35), tds_2011=c(-1.485, -1.485, -1.485), surv=c(5, 5, 5))

individuals <- individuals %>%
  calculate_qdiabeteshf(sex=sex, age=age, ethrisk=ethrisk, smoking=smoking, duration=duration_cat, type1=type1, cvd=cvd, renal=renal, af=af, hba1c=hba1c, cholhdl=cholhdl, sbp=sbp, bmi=bmi, town=tds_2011, surv=surv)

individuals <- individuals %>%
  mutate(centred_qdiabeteshf_lin_predictor=ifelse(sex=="male", qdiabeteshf_lin_predictor-0.464641698462435, qdiabeteshf_lin_predictor-0.468985206202883)) %>%
  mutate(`non-SGLT2i`=ifelse(sex=="male", (1-(0.965209400040619^exp(centred_qdiabeteshf_lin_predictor)))*100, (1-(0.965202142600253^exp(centred_qdiabeteshf_lin_predictor)))*100))


individuals <- individuals %>%
  mutate(qdhf_survival_cal=(100-`non-SGLT2i`)/100,
         qdhf_survival_cal_sglt2=qdhf_survival_cal^0.63,
         SGLT2i=100-(qdhf_survival_cal_sglt2*100),
         sglt2_benefit=qdhf_survival_cal_sglt2-qdhf_survival_cal)

individuals <- individuals %>%
  select(name, SGLT2i, `non-SGLT2i`) %>%
  pivot_longer(!name, names_to="drug", values_to="var")

individuals



setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/")


#Anthony
x <- individuals %>% filter(name=="Anthony")

pdis <- 
  ggplot(x, aes(x=drug, y=var, group=drug)) + 
  geom_point(shape=18,size=10,aes(colour=drug)) +
  geom_segment(aes(x=drug,xend=drug,y=0,yend=var,colour=drug,size=6)) +
  ylab("5-year risk of heart failure (%)") + xlab("") + 
  scale_y_continuous(limits=c(0,17),breaks=c(seq(0,17,by=4))) +
  coord_cartesian(ylim = c(0,40)) +
  geom_hline(yintercept=0,linetype="dashed") +
  coord_flip() +
  scale_color_manual(values=c("#998ec3","#f1a340"))+
  theme_bw() +
  theme(axis.text=element_text(size=20)) +
  #theme(axis.text.y= element_text(colour=col)) +
  theme(axis.title=element_text(size=20))+
  theme(legend.position = "none") + 
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(2.2), face = "bold")) 

pdf.options(reset = TRUE, onefile = FALSE)
pdf("Anthony.pdf",width=10,height=1.7)
pdis
dev.off() 

#Beryl
x <- individuals %>% filter(name=="Beryl")

pdis <- 
  ggplot(x, aes(x=drug, y=var, group=drug)) + 
  geom_point(shape=18,size=10,aes(colour=drug)) +
  geom_segment(aes(x=drug,xend=drug,y=0,yend=var,colour=drug,size=6)) +
  ylab("5-year risk of heart failure (%)") + xlab("") + 
  scale_y_continuous(limits=c(0,17),breaks=c(seq(0,17,by=4))) +
  coord_cartesian(ylim = c(0,40)) +
  geom_hline(yintercept=0,linetype="dashed") +
  coord_flip() +
  scale_color_manual(values=c("#998ec3","#f1a340"))+
  theme_bw() +
  theme(axis.text=element_text(size=20)) +
  #theme(axis.text.y= element_text(colour=col)) +
  theme(axis.title=element_text(size=20))+
  theme(legend.position = "none") + 
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(2.2), face = "bold")) 

pdf.options(reset = TRUE, onefile = FALSE)
pdf("Beryl.pdf",width=10,height=1.7)
pdis
dev.off()    

#Claire
x <- individuals %>% filter(name=="Claire")

pdis <- 
  ggplot(x, aes(x=drug, y=var, group=drug)) + 
  geom_point(shape=18,size=10,aes(colour=drug)) +
  geom_segment(aes(x=drug,xend=drug,y=0,yend=var,colour=drug,size=6)) +
  ylab("5-year risk of heart failure (%)") + xlab("") + 
  scale_y_continuous(limits=c(0,17),breaks=c(seq(0,17,by=4))) +
  coord_cartesian(ylim = c(0,40)) +
  geom_hline(yintercept=0,linetype="dashed") +
  coord_flip() +
  scale_color_manual(values=c("#998ec3","#f1a340"))+
  theme_bw() +
  theme(axis.text=element_text(size=20)) +
  #theme(axis.text.y= element_text(colour=col)) +
  theme(axis.title=element_text(size=20))+
  theme(legend.position = "none") + 
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(2.2), face = "bold")) 

pdf.options(reset = TRUE, onefile = FALSE)
pdf("Claire.pdf",width=10,height=1.7)
pdis
dev.off()    



individuals <- individuals %>%
  pivot_wider(id_cols=name, names_from=drug, values_from=var) %>%
  mutate(benefit=`non-SGLT2i`-SGLT2i)

individuals




# example at 20%

x <- data.frame(drug=c("non-SGLT2i", "SGLT2i"), var=c(20,18))

pdis <- 
  ggplot(x, aes(x=drug, y=var, group=drug)) + 
  geom_point(shape=18,size=10,aes(colour=drug)) +
  geom_segment(aes(x=drug,xend=drug,y=0,yend=var,colour=drug,size=6)) +
  ylab("Risk of outcome (%)") + xlab("") + 
  scale_y_continuous(limits=c(0,21),breaks=c(seq(0,21,by=4))) +
  coord_cartesian(ylim = c(0,40)) +
  geom_hline(yintercept=0,linetype="dashed") +
  coord_flip() +
  scale_color_manual(values=c("#998ec3","#f1a340"))+
  theme_bw() +
  theme(axis.text=element_text(size=20)) +
  #theme(axis.text.y= element_text(colour=col)) +
  theme(axis.title=element_text(size=20))+
  theme(legend.position = "none") + 
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(2.2), face = "bold")) 

pdf.options(reset = TRUE, onefile = FALSE)
pdf("example.pdf",width=10,height=1.7)
pdis
dev.off() 


x <- data.frame(drug=c("non-SGLT2i", "SGLT2i"), var=c(4,3.6))

pdis <- 
  ggplot(x, aes(x=drug, y=var, group=drug)) + 
  geom_point(shape=18,size=10,aes(colour=drug)) +
  geom_segment(aes(x=drug,xend=drug,y=0,yend=var,colour=drug,size=6)) +
  ylab("Risk of outcome (%)") + xlab("") + 
  scale_y_continuous(limits=c(0,21),breaks=c(seq(0,21,by=4))) +
  coord_cartesian(ylim = c(0,40)) +
  geom_hline(yintercept=0,linetype="dashed") +
  coord_flip() +
  scale_color_manual(values=c("#998ec3","#f1a340"))+
  theme_bw() +
  theme(axis.text=element_text(size=20)) +
  #theme(axis.text.y= element_text(colour=col)) +
  theme(axis.title=element_text(size=20))+
  theme(legend.position = "none") + 
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(2.2), face = "bold")) 

pdf.options(reset = TRUE, onefile = FALSE)
pdf("example2.pdf",width=10,height=1.7)
pdis
dev.off() 





