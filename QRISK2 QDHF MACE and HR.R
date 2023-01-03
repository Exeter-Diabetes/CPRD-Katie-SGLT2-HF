
# Calibrate QRISK2 to MACE outcome in subset of control arm (+ plot graphs + look at discrimination before and after in calibration cohort [by decile, include no. of events])
# Predict absolute benefit of SGLT2 in all based on 0.94 HR
# Plot graphs + look at discrimination of QRISK2 alone in control + SGLT2 arm + with HR in SGLT2 arm (by decile, include no. of events)
# Compare deciles of predicted benefit to actual

############################################################################################

# Setup
library(tidyverse)
library(gtsummary)
library(survival)
library(survminer)
library(broom)
rm(list=ls())


############################################################################################

# 1 Cohort selection and variable setup



## FILTERING

# Inclusion/exclusion criteria:
## a) T2Ds
## b) With HES linkage
## c) 1st instance
## d) Exclude if start drug within 91 days of registration

## e) Aged 18+
## f) SGLT2/DPP4/SU
## g) Initiated between 01/01/2013 and end of data (31/10/2020)
## h) Exclude if first line
## i) Exclude if also on insulin/GLP1/SGLT2 (except SGLT2 arm)/TZD
## j) No CVD (NICE definition: angina, IHD, MI, PAD, revasc, stroke) or HF before index date
## k) Exclude if CKD (stage 3a-5) before index date
## l) Exclude if missing QRISK2 or QDHF
### NB: this includes if any required variables missing (smoking status, baseline HbA1c) or if values out of range (age<25 or >84, cholHDL<1 or >11, HbA1c<40 or >150, SBP<70 or >210, BMI<20 - QDHF has not been calculated in these cases))
### will also exclude anyone without QRISK2 score (missing if missing smoking status or age/cholHDL/SBP/BMI outside of range [weirdly cholHDL range for QRISK2 is 1-12 vs 1-11 for QDHF])


# Start with saved dataset of 1st instance data for T2Ds with HES data, excluding drugs started within 91 days following registration (a-d above)
# 20221212_t2d_1stinstance is identical to 20221110_t2d_1stinstance except corrected dates for ICD10-only conditions including primary HHF in 20221116 version, and removed QDHF score for those aged>84 or with BMI<20, and added 5 year and 10 year QRISK2, added all-cause emergency hospitalisation as a postdrug comorbitidity, and added HF and CV death, and added statins

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2022/1 SGLT2 CVD project/Raw data/")
load("20221212_t2d_1stinstance.Rda")


# Keep those aged >=18 and within study period and second line or later (e-h above)
cohort <- t2d_1stinstance %>%
  filter(dstartdate_age>=18 &
           (drugclass=="SGLT2" | drugclass=="DPP4" | drugclass=="SU") &
           dstartdate>=as.Date("2013-01-01") &
           drugline_all!=1) %>%
  mutate(studydrug=ifelse(drugclass=="SGLT2", "SGLT2", ifelse(drugclass=="GLP1", "GLP1", "DPP4SU")))

cohort %>% group_by(studydrug) %>% distinct(patid) %>% summarise(count=n())
#SGLT2: 99,030
#DPP4SU: 202,272


# Remove if on insulin at start (i above)
cohort <- cohort %>%
  filter(INS==0)

cohort %>% group_by(studydrug) %>% distinct(patid) %>% summarise(count=n())
#SGLT2: 85,103
#DPP4SU: 191,393


# Remove if on GLP1/SGLT2 (except SGL2 arm)/TZD (i above)
cohort <- cohort %>%
  filter(GLP1==0 & TZD==0 & (drugclass=="SGLT2" | SGLT2==0))

cohort %>% group_by(studydrug) %>% distinct(patid) %>% summarise(count=n())
#SGLT2: 78,236
#DPP4SU: 181,634


# Remove if CVD or HF before index date (j above)
cohort <- cohort %>%
  mutate(predrug_cvd=ifelse(predrug_angina==1 | predrug_ihd==1 | predrug_myocardialinfarction==1 | predrug_pad==1 | predrug_revasc==1 | predrug_stroke==1, 1, 0)) %>%
  filter(predrug_cvd==0 & predrug_heartfailure==0)

cohort %>% group_by(studydrug) %>% distinct(patid) %>% summarise(count=n())
#SGLT2: 61,027
#DPP4SU: 126,819


# Remove if CKD before index date (k above)
cohort <- cohort %>%
  filter(is.na(preckdstage) | (preckdstage!="stage_3a" & preckdstage!="stage_3b" & preckdstage!="stage_4" & preckdstage!="stage_5"))

cohort %>% group_by(studydrug) %>% distinct(patid) %>% summarise(count=n())
#SGLT2: 59,318
#DPP4SU: 111,417


# Remove if don't have QDHF variables (also removes those without QRISK2)
cohort <- cohort %>%
  filter(!is.na(qdiabeteshf_5yr_score))

cohort %>% group_by(studydrug) %>% distinct(patid) %>% summarise(count=n())
#SGLT2: 48,304
#DPP4SU: 90,904


# For DPP4SU arm, some people will be in dataset twice with both DPP4 and SU - choose earliest period
## If start both on same day, choose 1 at random
cohort <- cohort %>%
  group_by(patid, studydrug) %>%
  mutate(earliest_start=min(dstartdate, na.rm=TRUE)) %>%
  filter(dstartdate==earliest_start) %>%
  mutate(id=row_number()) %>%
  filter(id==min(id, na.rm=TRUE)) %>%
  ungroup()



## Use all SGLT2, GLP1 + TZD starts to code up later censoring
### Also get latest SGLT2 stop dates before drug start for DPP4SU arm as need this for sensitivity analysis
### 20221212_t2d_all_drug_periods is identical to 20221110_t2d_all_drug_periods

load("20221212_t2d_all_drug_periods.Rda")

later_sglt2 <- cohort %>%
  select(patid, dstartdate) %>%
  inner_join((t2d_all_drug_periods %>%
                filter(drugclass=="SGLT2") %>%
                select(patid, next_sglt2=dstartdate)), by="patid") %>%
  filter(next_sglt2>dstartdate) %>%
  group_by(patid, dstartdate) %>%
  summarise(next_sglt2_start=min(next_sglt2, na.rm=TRUE)) %>%
  ungroup()


later_glp1 <- cohort %>%
  select(patid, dstartdate) %>%
  inner_join((t2d_all_drug_periods %>%
                filter(drugclass=="GLP1") %>%
                select(patid, next_glp1=dstartdate)), by="patid") %>%
  filter(next_glp1>dstartdate) %>%
  group_by(patid, dstartdate) %>%
  summarise(next_glp1_start=min(next_glp1, na.rm=TRUE)) %>%
  ungroup()


later_tzd <- cohort %>%
  select(patid, dstartdate) %>%
  inner_join((t2d_all_drug_periods %>%
                filter(drugclass=="TZD") %>%
                select(patid, next_tzd=dstartdate)), by="patid") %>%
  filter(next_tzd>dstartdate) %>%
  group_by(patid, dstartdate) %>%
  summarise(next_tzd_start=min(next_tzd, na.rm=TRUE)) %>%
  ungroup()


last_sglt2_stop <- cohort %>%
  select(patid, dstartdate) %>%
  inner_join((t2d_all_drug_periods %>%
                filter(drugclass=="SGLT2") %>%
                select(patid, last_sglt2=dstopdate)), by="patid") %>%
  filter(last_sglt2<dstartdate) %>%
  group_by(patid, dstartdate) %>%
  summarise(last_sglt2_stop=min(last_sglt2, na.rm=TRUE)) %>%
  ungroup()

cohort <- cohort %>%
  left_join(later_sglt2, by=c("patid", "dstartdate")) %>%
  left_join(later_glp1, by=c("patid", "dstartdate")) %>%
  left_join(later_tzd, by=c("patid", "dstartdate")) %>%
  left_join(last_sglt2_stop, by=c("patid", "dstartdate"))



# Keep variables of interest
## Also tidy up gender, ncurrtx, drugline and ethnicity variables
## Add death cause variables
cohort <- cohort %>%
  
  mutate(malesex=ifelse(gender==1, 1, 0),
         ncurrtx=DPP4+GLP1+MFN+SU+SGLT2+TZD+INS,          #INS, GLP1 and TZD should be 0 but include anyway; ignore Acarbose and Glinide
         drugline_all=as.factor(ifelse(drugline_all>=5, 5, drugline_all)),
         drugsubstances=ifelse(grepl("&", drugsubstances), NA, drugsubstances),
         ethnicity_5cat_decoded=case_when(ethnicity_5cat==0 ~"White",
                                          ethnicity_5cat==1 ~"South Asian",
                                          ethnicity_5cat==2 ~"Black",
                                          ethnicity_5cat==3 ~"Other",
                                          ethnicity_5cat==4 ~"Mixed"),
         cv_death_date_any_cause=if_else(!is.na(death_date) & !is.na(cv_death_any_cause) & cv_death_any_cause==1, death_date, as.Date(NA)),
         cv_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(cv_death_primary_cause) & cv_death_primary_cause==1, death_date, as.Date(NA)),
         hf_death_date_any_cause=if_else(!is.na(death_date) & !is.na(hf_death_any_cause) & hf_death_any_cause==1, death_date, as.Date(NA)),
         hf_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(hf_death_primary_cause) & hf_death_primary_cause==1, death_date, as.Date(NA))) %>%
  
  select(patid, malesex, ethnicity_5cat_decoded, imd2015_10, regstartdate, gp_record_end, death_date, drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, prehba1c, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preckdstage, qrisk2_smoking_cat, postdrug_first_myocardialinfarction, postdrug_first_primary_incident_mi, postdrug_first_stroke, postdrug_first_primary_incident_stroke, postdrug_first_heartfailure, postdrug_first_primary_hhf, postdrug_first_all_cause_hosp, next_glp1_start, next_sglt2_start, next_tzd_start, last_sglt2_stop, cv_death_date_primary_cause, cv_death_date_any_cause, hf_death_date_primary_cause, hf_death_date_any_cause, qrisk2_lin_predictor, qrisk2_5yr_score, qrisk2_10yr_score, qdiabeteshf_lin_predictor, qdiabeteshf_5yr_score, contains("statins"))

rm(list=setdiff(ls(), "cohort"))


############################################################################################

# 2 Make outcome and survival analysis variables for MACE and HF outcomes
## Broad MACE (hospitalisation/death any cause or in GP records)

cohort <- cohort %>%
  
  mutate(
    
    # broad_mace = GP (MI/stroke) + expanded HES codelist (MI/stroke; any cause) + CV death (any cause)
    postdrug_broad_mace=pmin(postdrug_first_myocardialinfarction, postdrug_first_stroke, cv_death_date_any_cause, na.rm=TRUE),
    
  )


# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Start SGLT2 (not SGLT2 arm), GLP1, or TZD
## End of GP records
## Broad MACE outcome

cohort <- cohort %>%
  mutate(five_years_post_dstart=dstartdate+(365.25*5),
         
         mace_broad_censdate=if_else(studydrug=="SGLT2",
                                     pmin(five_years_post_dstart,
                                          death_date,
                                          next_glp1_start,
                                          next_tzd_start,
                                          gp_record_end,
                                          postdrug_broad_mace, na.rm=TRUE),
                                     if_else(studydrug=="DPP4SU",
                                             pmin(five_years_post_dstart,
                                                  death_date,
                                                  next_sglt2_start,
                                                  next_glp1_start,
                                                  next_tzd_start,
                                                  gp_record_end,
                                                  postdrug_broad_mace, na.rm=TRUE),
                                             as.Date(NA))),
         
         mace_broad_censvar=ifelse(!is.na(postdrug_broad_mace) & mace_broad_censdate==postdrug_broad_mace, 1, 0),
         
         mace_broad_censtime_yrs=as.numeric(difftime(mace_broad_censdate, dstartdate, unit="days"))/365.25)


############################################################################################

# 3 Calibrate QRISK2 to subset of control arm

# Assign random 20% of DPP4SU arm as calibration cohort and remove from main cohort
set.seed(123)

cal_cohort <- cohort %>%
  filter(studydrug=="DPP4SU") %>%
  slice_sample(prop=0.2)
#18,180

noncal_cohort <- cohort %>%
  anti_join(cal_cohort, by=c("patid", "dstartdate", "studydrug"))
table(noncal_cohort$studydrug)
#SGLT2: 48,304
#DPP4SU: 72,724


# Compare deciles of QRISK2 score (predicted) with observed

## Get mean predicted probabilities for each decile
cal_cohort <- cal_cohort %>%
  mutate(qrisk2_decile=ntile(qrisk2_5yr_score, 10))

predicted <- cal_cohort %>% group_by(qrisk2_decile) %>% summarise(mean_pred=mean(qrisk2_5yr_score)/100)


## Find actual observed probabilities by QDHF decile

observed <- survfit(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ qrisk2_decile, data=cal_cohort) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(estimate=1-estimate,
         lower_ci=1-conf.low,
         upper_ci=1-conf.high) %>%
  select(observed=estimate, lower_ci, upper_ci, strata)


## Plot predicted vs observed
obs_v_pred <- cbind(observed, predicted)

ggplot(obs_v_pred, aes(x=qrisk2_decile)) + 
  geom_point(aes(y = mean_pred*100), color = "darkred") + 
  geom_point(aes(y = observed*100), color="steelblue") +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1,width=0.25,size=1,colour="steelblue") + theme_bw() +
  xlab("Tenth of predicted risk") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Initial calibration")

# Looks good for deciles 1-4, underestimates for rest



# Re-estimate baseline hazard for females
## Original: 0.994671821594238
cal_females <- cal_cohort %>%
  filter(malesex==0)

recal_mod <- coxph(Surv(mace_broad_censtime_yrs,mace_broad_censvar)~offset(qrisk2_lin_predictor), data=cal_females)
female_surv <- summary(survfit(recal_mod),time=5)$surv
female_surv
# 0.9536599


# Re-estimate baseline hazard for males
## Original: 0.981611728668213
cal_males <- cal_cohort %>%
  filter(malesex==1)

recal_mod <- coxph(Surv(mace_broad_censtime_yrs,mace_broad_censvar)~offset(qrisk2_lin_predictor), data=cal_males)
male_surv <- summary(survfit(recal_mod),time=5)$surv
male_surv
# 0.9322214



# Recalculate QRISK2 and compare to observed
cal_cohort <- cal_cohort %>%
  group_by(malesex) %>%
  mutate(centred_lin_predictor=qrisk2_lin_predictor-mean(qrisk2_lin_predictor)) %>%
  ungroup() %>%
  mutate(qrisk2_5yr_score_cal=ifelse(malesex==1, (1-(male_surv^exp(centred_lin_predictor)))*100, (1-(female_surv^exp(centred_lin_predictor)))*100))

cal_cohort <- cal_cohort %>%
  mutate(qrisk2_decile=ntile(qrisk2_5yr_score_cal, 10))

predicted <- cal_cohort %>% group_by(qrisk2_decile) %>% summarise(mean_pred=mean(qrisk2_5yr_score_cal)/100)


## Find actual observed probabilities by QDHF decile

observed <- survfit(Surv(mace_broad_censtime_yrs,mace_broad_censvar) ~ qrisk2_decile, data=cal_cohort) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(estimate=1-estimate,
         lower_ci=1-conf.low,
         upper_ci=1-conf.high) %>%
  select(observed=estimate, lower_ci, upper_ci, strata)


## Plot predicted vs observed
obs_v_pred <- cbind(observed, predicted)

ggplot(obs_v_pred, aes(x=qrisk2_decile)) + 
  geom_point(aes(y = mean_pred*100), color = "darkred") + 
  geom_point(aes(y = observed*100), color="steelblue") +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1,width=0.25,size=1,colour="steelblue") + theme_bw() +
  xlab("Tenth of predicted risk") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Initial calibration")


# Recalculate QRISK2 score for rest of cohort
noncal_cohort <- noncal_cohort %>%
  group_by(malesex) %>%
  mutate(centred_lin_predictor=qrisk2_lin_predictor-mean(qrisk2_lin_predictor)) %>%
  ungroup() %>%
  mutate(qrisk2_5yr_score_cal=ifelse(malesex==1, (1-(male_surv^exp(centred_lin_predictor)))*100, (1-(female_surv^exp(centred_lin_predictor)))*100))


# Calculate the C-statistic

cal_cohort <- cal_cohort %>%
  mutate(qrisk2_5yr_survival=(100-qrisk2_5yr_score)/100,
         qrisk2_5yr_survival_cal=(100-qrisk2_5yr_score_cal)/100)

noncal_cohort <- noncal_cohort %>%
  mutate(qrisk2_5yr_survival=(100-qrisk2_5yr_score)/100,
         qrisk2_5yr_survival_cal=(100-qrisk2_5yr_score_cal)/100)


## Calibration cohort
surv_mod <- coxph(Surv(mace_broad_censtime_yrs,mace_broad_censvar)~qrisk2_5yr_survival,method="breslow",data=cal_cohort)
round(summary(surv_mod)$concordance[1],2)
round(summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2]),2)
round(summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2]),2)
# 0.69 (0.67-0.71)

surv_mod <- coxph(Surv(mace_broad_censtime_yrs,mace_broad_censvar)~qrisk2_5yr_survival_cal,method="breslow",data=cal_cohort)
round(summary(surv_mod)$concordance[1],2)
round(summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2]),2)
round(summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2]),2)
# 0.68 (0.66-0.7)


## Non-calibration cohort
surv_mod <- coxph(Surv(mace_broad_censtime_yrs,mace_broad_censvar)~qrisk2_5yr_survival,method="breslow",data=noncal_cohort)
round(summary(surv_mod)$concordance[1],2)
round(summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2]),2)
round(summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2]),2)
# 0.67 (0.67-0.68)

surv_mod <- coxph(Surv(mace_broad_censtime_yrs,mace_broad_censvar)~qrisk2_5yr_survival_cal,method="breslow",data=noncal_cohort)
round(summary(surv_mod)$concordance[1],2)
round(summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2]),2)
round(summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2]),2)
# 0.67 (0.67-0.68)



############################################################################################

# 4 Look at survival benefit from SGLT2s when add hazard ratio from trials meta-analysis
# https://jamanetwork.com/journals/jamacardiology/fullarticle/2771459


# Define survival, then add HR from trials
noncal_cohort <- noncal_cohort %>%
  mutate(qrisk2_5yr_survival_sglt2=qrisk2_5yr_survival_cal^0.94) 


# Estimate survival benefit
noncal_cohort <- noncal_cohort %>% mutate(sglt2_benefit=100*(qrisk2_5yr_survival_sglt2-qrisk2_5yr_survival_cal))
summary(noncal_cohort$sglt2_benefit)
# mean = 0.4% (absolute reduction in HR risk)


# Find actual predicted survival based on study drug
noncal_cohort <- noncal_cohort %>%
  mutate(actual_model_survival=ifelse(studydrug=="DPP4SU", qrisk2_5yr_survival_cal, qrisk2_5yr_survival_sglt2))
summary((1-noncal_cohort$actual_model_survival)*100)
# Mean=7.1%

## And find observed at 5 years for whole cohort
obs_surv <- summary(survfit(Surv(mace_broad_censtime_yrs,mace_broad_censvar)~1, data=noncal_cohort),time=5)$surv
1-obs_surv
# 6.8%



# Look at calibration for SGLT2 arm with HR from trials
sglt2_cohort <- noncal_cohort %>%
  filter(studydrug=="SGLT2")

surv_mod <- coxph(Surv(mace_broad_censtime_yrs,mace_broad_censvar)~qrisk2_5yr_survival_sglt2,method="breslow",data=sglt2_cohort)
round(summary(surv_mod)$concordance[1],2)
round(summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2]),2)
round(summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2]),2)
# 0.67 (0.66-0.69)







#Paper like plot
library(rms)
attach(sglt2_cohort)
plot(0,0)
recal <- data.frame(groupkm(qrisk2_5yr_survival_sglt2, S = Surv(mace_broad_censtime_yrs,mace_broad_censvar), g=10, u=5, pl=T, add=T, lty=0, cex.subtitle=FALSE))
recal$pred.decile=c(10,9,8,7,6,5,4,3,2,1)
head(recal)

recal.obs <- data.frame(est=(1-recal$KM)*100,
                        conf.low=(1-recal$KM-(1.96*recal$std.err))*100,
                        conf.high=(1-recal$KM+(1.96*recal$std.err))*100,
                        pred.decile=recal$pred.decile,
                        type="Observed")

recal.pred <- data.frame(est=(1-recal$x)*100,
                         conf.low=NA,conf.high=NA,pred.decile=recal$pred.decile,
                         type="Predicted")

recal <- rbind(recal.obs,recal.pred)

ggplot(data=recal, aes(x=pred.decile,y=est,group=type,colour=type)) + 
  geom_point(position=position_dodge(width=-0.3)) + 
  geom_errorbar(aes(ymax=conf.high,ymin=conf.low),position=position_dodge(width=-0.3))+ 
  theme_bw() +
  scale_color_manual(values=c("steelblue","darkred")) +
  xlab("Quintiles of predicted risk") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("SGLT2 group model calibration")

detach("package:rms", unload=TRUE)




Look at hazard ratio

table(cohort$studydrug)
#SGLT2: 48,304
#DPP4SU: 90,904



# MACE

cohort %>% group_by(studydrug) %>% summarise(time=median(mace_broad_censtime_yrs))

cohort %>% group_by(studydrug) %>% summarise(events=sum(mace_broad_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ studydrug, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_10yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

############################################################################################

# 4 Look at MACE by arm vs QRISK2

## Use ~deciles of QRISK2, but convert to actual values
quantile(cohort$qrisk2_5yr_score, probs = seq(.1, .9, by = .1))
# 3, 5, 6.5, 8, 9.5, 11.5, 13.5, 16, 20

cohort <- cohort %>%
  mutate(qrisk2_cat=case_when(
    qrisk2_5yr_score<=3 ~ 1,
    qrisk2_5yr_score>3 & qrisk2_5yr_score<=5 ~ 2,
    qrisk2_5yr_score>5 & qrisk2_5yr_score<=6.5 ~ 3,
    qrisk2_5yr_score>6.5 & qrisk2_5yr_score<=8 ~ 4,
    qrisk2_5yr_score>8 & qrisk2_5yr_score<=9.5 ~ 5,
    qrisk2_5yr_score>9.5 & qrisk2_5yr_score<=11.5 ~ 6,
    qrisk2_5yr_score>11.5 & qrisk2_5yr_score<=13.5 ~ 7,
    qrisk2_5yr_score>13.5 & qrisk2_5yr_score<=16 ~ 8,
    qrisk2_5yr_score>16 & qrisk2_5yr_score<=20 ~ 9,
    qrisk2_5yr_score>20 ~ 10
  ))

table(cohort$qrisk2_cat, useNA="always")



# Overall (both arms together)

## Get mean predicted probabilities for each QRISK2 category
predicted_overall <- cohort %>%
  group_by(qrisk2_cat) %>%
  summarise(mean_pred_overall=mean(qrisk2_5yr_score)/100)

## Find actual observed probabilities by QRISK2 category
observed_overall <- survfit(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ qrisk2_cat, data=cohort) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_overall=1-estimate,
         lower_ci_overall=1-conf.low,
         upper_ci_overall=1-conf.high) %>%
  select(observed_overall, lower_ci_overall, upper_ci_overall, strata)

## Plot
obs_v_pred <- cbind(predicted_overall, observed_overall)


ggplot(obs_v_pred, aes(x=qrisk2_cat)) + 
  geom_point(aes(y = mean_pred_overall*100), color = "darkred") +
  geom_point(aes(y = observed_overall*100), color="steelblue") +
  geom_errorbar(aes(ymax=upper_ci_overall*100,ymin=lower_ci_overall*100),alpha=1,width=0.25,size=1,colour="steelblue") +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Overall")
dev.off()



# By arm

## Get mean predicted probabilities for each QRISK2 category
predicted <- cohort %>%
  group_by(qrisk2_cat, studydrug) %>%
  summarise(mean_pred=mean(qrisk2_5yr_score)/100)

## Find actual observed probabilities by QRISK2 category
observed_dpp4su <- survfit(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_dpp4su=1-estimate,
         lower_ci_dpp4su=1-conf.low,
         upper_ci_dpp4su=1-conf.high) %>%
  select(observed_dpp4su, lower_ci_dpp4su, upper_ci_dpp4su, strata)


observed_sglt2 <- survfit(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_sglt2=1-estimate,
         lower_ci_sglt2=1-conf.low,
         upper_ci_sglt2=1-conf.high) %>%
  select(observed_sglt2, lower_ci_sglt2, upper_ci_sglt2, strata)


obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug=="DPP4SU")), observed_dpp4su),
  cbind((predicted %>% filter(studydrug=="SGLT2")), observed_sglt2)
) %>%
  mutate(observed=coalesce(observed_dpp4su, observed_sglt2),
         lower_ci=coalesce(lower_ci_dpp4su, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_dpp4su, upper_ci_sglt2))



dodge <- position_dodge(width=0.3)  
ggplot(obs_v_pred, aes(x=qrisk2_cat, group=studydrug, color=studydrug)) + 
  geom_point(aes(y = mean_pred*100), position=dodge) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1,width=0.25,size=1, position=dodge) +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,25,by=5)), limits=c(0,26)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("By drug")
dev.off()



############################################################################################

# 5 Look at individual MACE components

cohort <- cohort %>%
  mutate(which_mace=ifelse(!is.na(postdrug_broad_mace) & !is.na(postdrug_first_myocardialinfarction) & !is.na(postdrug_first_stroke) & !is.na(cv_death_date_primary_cause) & postdrug_broad_mace==postdrug_first_myocardialinfarction & postdrug_broad_mace==postdrug_first_stroke & postdrug_broad_mace==cv_death_date_primary_cause, "all three",
                           ifelse(!is.na(postdrug_broad_mace) & !is.na(postdrug_first_myocardialinfarction) & !is.na(postdrug_first_stroke) & postdrug_broad_mace==postdrug_first_myocardialinfarction & postdrug_broad_mace==postdrug_first_stroke, "mi & stroke",
                                  ifelse(!is.na(postdrug_broad_mace) & !is.na(postdrug_first_myocardialinfarction) & !is.na(cv_death_date_primary_cause) & postdrug_broad_mace==postdrug_first_myocardialinfarction & postdrug_broad_mace==cv_death_date_primary_cause, "mi & death",
                                         ifelse(!is.na(postdrug_broad_mace) & !is.na(postdrug_first_stroke) & !is.na(cv_death_date_primary_cause) & postdrug_broad_mace==postdrug_first_stroke & postdrug_broad_mace==cv_death_date_primary_cause, "stroke & death",
                                                ifelse(!is.na(postdrug_broad_mace) & !is.na(postdrug_first_myocardialinfarction) & postdrug_broad_mace==postdrug_first_myocardialinfarction, "mi",
                                                       ifelse(!is.na(postdrug_broad_mace) & !is.na(postdrug_first_stroke) & postdrug_broad_mace==postdrug_first_stroke, "stroke",
                                                              ifelse(!is.na(postdrug_broad_mace) & !is.na(cv_death_date_primary_cause) & postdrug_broad_mace==cv_death_date_primary_cause, "cv death", NA))))))))



prop.table(table(cohort$which_mace))
# 13% CV death, 44% MI, 42% stroke

censvaryes <- cohort %>%
  filter(mace_broad_censvar==1)

prop.table(table(censvaryes$which_mace))
# 13% CV death, 43% MI, 43% stroke


# Make new survival variables; censor at postdrug_broad_mace
cohort <- cohort %>%
  mutate(mi_censvar=ifelse(!is.na(postdrug_first_myocardialinfarction) & mace_broad_censdate==postdrug_first_myocardialinfarction, 1, 0),
         stroke_censvar=ifelse(!is.na(postdrug_first_stroke) & mace_broad_censdate==postdrug_first_stroke, 1, 0),
         cv_death_censvar=ifelse(!is.na(cv_death_date_any_cause) & mace_broad_censdate==cv_death_date_any_cause, 1, 0))



# By arm

## Find actual observed probabilities by QRISK2 category

### MI

mi_dpp4su <- survfit(Surv(mace_broad_censtime_yrs, mi_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_dpp4su=1-estimate,
         lower_ci_dpp4su=1-conf.low,
         upper_ci_dpp4su=1-conf.high) %>%
  select(observed_dpp4su, lower_ci_dpp4su, upper_ci_dpp4su, strata)

mi_sglt2 <- survfit(Surv(mace_broad_censtime_yrs, mi_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_sglt2=1-estimate,
         lower_ci_sglt2=1-conf.low,
         upper_ci_sglt2=1-conf.high) %>%
  select(observed_sglt2, lower_ci_sglt2, upper_ci_sglt2, strata)


obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug=="DPP4SU")), mi_dpp4su),
  cbind((predicted %>% filter(studydrug=="SGLT2")), mi_sglt2)
) %>%
  mutate(observed=coalesce(observed_dpp4su, observed_sglt2),
         lower_ci=coalesce(lower_ci_dpp4su, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_dpp4su, upper_ci_sglt2))



dodge <- position_dodge(width=0.3)  
ggplot(obs_v_pred, aes(x=qrisk2_cat, group=studydrug, color=studydrug)) + 
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1,width=0.25,size=1, position=dodge) +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,10,by=2)), limits=c(0,10))+
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("MI only")
dev.off()



### stroke

stroke_dpp4su <- survfit(Surv(mace_broad_censtime_yrs, stroke_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_dpp4su=1-estimate,
         lower_ci_dpp4su=1-conf.low,
         upper_ci_dpp4su=1-conf.high) %>%
  select(observed_dpp4su, lower_ci_dpp4su, upper_ci_dpp4su, strata)

stroke_sglt2 <- survfit(Surv(mace_broad_censtime_yrs, stroke_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_sglt2=1-estimate,
         lower_ci_sglt2=1-conf.low,
         upper_ci_sglt2=1-conf.high) %>%
  select(observed_sglt2, lower_ci_sglt2, upper_ci_sglt2, strata)


obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug=="DPP4SU")), stroke_dpp4su),
  cbind((predicted %>% filter(studydrug=="SGLT2")), stroke_sglt2)
) %>%
  mutate(observed=coalesce(observed_dpp4su, observed_sglt2),
         lower_ci=coalesce(lower_ci_dpp4su, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_dpp4su, upper_ci_sglt2))



dodge <- position_dodge(width=0.3)  
ggplot(obs_v_pred, aes(x=qrisk2_cat, group=studydrug, color=studydrug)) + 
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1,width=0.25,size=1, position=dodge) +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,10,by=2)), limits=c(0,10))+
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Stroke only")
dev.off()


### CV death

cvdeath_dpp4su <- survfit(Surv(mace_broad_censtime_yrs, cv_death_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_dpp4su=1-estimate,
         lower_ci_dpp4su=1-conf.low,
         upper_ci_dpp4su=1-conf.high) %>%
  select(observed_dpp4su, lower_ci_dpp4su, upper_ci_dpp4su, strata)

cvdeath_sglt2 <- survfit(Surv(mace_broad_censtime_yrs, cv_death_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_sglt2=1-estimate,
         lower_ci_sglt2=1-conf.low,
         upper_ci_sglt2=1-conf.high) %>%
  select(observed_sglt2, lower_ci_sglt2, upper_ci_sglt2, strata)


obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug=="DPP4SU")), cvdeath_dpp4su),
  cbind((predicted %>% filter(studydrug=="SGLT2")), cvdeath_sglt2)
) %>%
  mutate(observed=coalesce(observed_dpp4su, observed_sglt2),
         lower_ci=coalesce(lower_ci_dpp4su, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_dpp4su, upper_ci_sglt2))



dodge <- position_dodge(width=0.3)  
ggplot(obs_v_pred, aes(x=qrisk2_cat, group=studydrug, color=studydrug)) + 
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1,width=0.25,size=1, position=dodge) +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,10,by=2)), limits=c(0,10))+
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("CV death only")
dev.off()



############################################################################################

# 6 Calibration

## DPP4SU
dpp4su <- cohort %>% filter(studydrug=="DPP4SU")

# Re-estimate baseline hazard for females
## Original: 0.994671821594238
dpp4su_females <- dpp4su %>%
  filter(malesex==0)

recal_mod <- coxph(Surv(mace_broad_censtime_yrs,mace_broad_censvar)~offset(qrisk2_lin_predictor), data=dpp4su_females)
female_surv <- summary(survfit(recal_mod),time=5)$surv
female_surv
# 0.9537195


# Re-estimate baseline hazard for males
## Original: 0.989570081233978
dpp4su_males <- dpp4su %>%
  filter(malesex==1)

recal_mod <- coxph(Surv(mace_broad_censtime_yrs,mace_broad_censvar)~offset(qrisk2_lin_predictor), data=dpp4su_males)
male_surv <- summary(survfit(recal_mod),time=5)$surv
male_surv
# 0.9361302


# Recalculate QRISK2 and compare to observed
dpp4su <- dpp4su %>%
  group_by(malesex) %>%
  mutate(centred_lin_predictor=qrisk2_lin_predictor-mean(qrisk2_lin_predictor)) %>%
  ungroup() %>%
  mutate(qrisk2_5yr_score_cal=ifelse(malesex==1, (1-(male_surv^exp(centred_lin_predictor)))*100, (1-(female_surv^exp(centred_lin_predictor)))*100))


predicted_dpp4su_new <- dpp4su %>% group_by(qrisk2_cat) %>% summarise(mean_pred_cal=mean(qrisk2_5yr_score_cal)/100)

## Plot predicted vs observed
obs_v_pred <- cbind((predicted %>% filter(studydrug=="DPP4SU") %>% left_join(predicted_dpp4su_new, by="qrisk2_cat")),
                    observed_dpp4su)

ggplot(obs_v_pred, aes(x=qrisk2_cat)) + 
  geom_point(aes(y = mean_pred*100), color = "darkred") +
  geom_point(aes(y = observed_dpp4su*100), color="steelblue") +
  geom_errorbar(aes(ymax=upper_ci_dpp4su*100,ymin=lower_ci_dpp4su*100),alpha=1,width=0.25,size=1,colour="steelblue") +
  theme_bw() +
  geom_point(aes(y = mean_pred_cal*100), color = "black", shape=4, size=4, stroke=2) +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("DPP4SU")
dev.off()




## SGLT2
sglt2 <- cohort %>% filter(studydrug=="SGLT2")

# Re-estimate baseline hazard for females
## Original: 0.994671821594238
sglt2_females <- sglt2 %>%
  filter(malesex==0)

recal_mod <- coxph(Surv(mace_broad_censtime_yrs,mace_broad_censvar)~offset(qrisk2_lin_predictor), data=sglt2_females)
female_surv <- summary(survfit(recal_mod),time=5)$surv
female_surv
# 0.9623085


# Re-estimate baseline hazard for males
## Original: 0.989570081233978
sglt2_males <- sglt2 %>%
  filter(malesex==1)

recal_mod <- coxph(Surv(mace_broad_censtime_yrs,mace_broad_censvar)~offset(qrisk2_lin_predictor), data=sglt2_males)
male_surv <- summary(survfit(recal_mod),time=5)$surv
male_surv
# 0.9442052


# Recalculate QRISK2 and compare to observed
sglt2 <- sglt2 %>%
  group_by(malesex) %>%
  mutate(centred_lin_predictor=qrisk2_lin_predictor-mean(qrisk2_lin_predictor)) %>%
  ungroup() %>%
  mutate(qrisk2_5yr_score_cal=ifelse(malesex==1, (1-(male_surv^exp(centred_lin_predictor)))*100, (1-(female_surv^exp(centred_lin_predictor)))*100))


predicted_sglt2_new <- sglt2 %>% group_by(qrisk2_cat) %>% summarise(mean_pred_cal=mean(qrisk2_5yr_score_cal)/100)

## Plot predicted vs observed
obs_v_pred <- cbind((predicted %>% filter(studydrug=="SGLT2") %>% left_join(predicted_sglt2_new, by="qrisk2_cat")),
                    observed_sglt2)

ggplot(obs_v_pred, aes(x=qrisk2_cat)) + 
  geom_point(aes(y = mean_pred*100), color = "darkred") +
  geom_point(aes(y = observed_sglt2*100), color="steelblue") +
  geom_errorbar(aes(ymax=upper_ci_sglt2*100,ymin=lower_ci_sglt2*100),alpha=1,width=0.25,size=1,colour="steelblue") +
  theme_bw() +
  geom_point(aes(y = mean_pred_cal*100), color = "black", shape=4, size=4, stroke=2) +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("SGLT2")
dev.off()

