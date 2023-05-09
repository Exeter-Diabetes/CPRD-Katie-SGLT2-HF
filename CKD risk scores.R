
# Exploring 2 x CKDPC risk scores

## 1 Look at counts with missing scores - not in cohort definition at the moment

## 2 (Uncalibrated) plots

## 3 Calibrated plots

## 4 Quartile plots


############################################################################################

# Setup
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(patchwork)
library(rms)
library(cowplot)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

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
  
  select(patid, malesex, ethnicity_5cat_decoded, imd2015_10, regstartdate, gp_record_end, death_date, drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preacr, preacr_from_separate, preweight, prehba1c, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preegfr, preckdstage, smoking_cat, qrisk2_smoking_cat, contains("cens"), qrisk2_lin_predictor, qrisk2_5yr_score, qdiabeteshf_lin_predictor, qdiabeteshf_5yr_score, starts_with("ckdpc"), last_sglt2_stop, contains("statins"), prepcr)

rm(list=setdiff(ls(), "cohort"))


############################################################################################

# 1 Look at missing counts

missing <- cohort %>%
  mutate(uacr=ifelse(!is.na(preacr), preacr, ifelse(!is.na(preacr_from_separate), preacr_from_separate, NA)),
         uacr_val=ifelse(is.na(uacr), "missing",
                         ifelse(uacr==0, "0",
                                ifelse(uacr>0 & uacr<0.6, "below range", "in range or above"))),
         egfr_under_60=ifelse(!is.na(preegfr) & preegfr<60, 1, 0),
         missing_egfr=ifelse(is.na(preegfr), 1, 0),
         hba1c_val=ifelse(is.na(prehba1c), "missing",
                          ifelse(prehba1c<42, "under",
                                 ifelse(prehba1c>97, "over", "in range"))),
         no_smoking_cat=ifelse(is.na(smoking_cat), 1, 0),
         bmi_val=ifelse(is.na(prebmi), "missing", 
                        ifelse(prebmi<20, "under",
                               ifelse(prebmi>40, "over", "in range"))))

prop.table(table(missing$uacr_val))

missing <- missing %>%
  mutate(uacr_val=ifelse(uacr_val=="missing" & !is.na(prepcr), "pcr_instead", uacr_val))


prop.table(table(missing$egfr_under_60))

prop.table(table(missing$missing_egfr))

prop.table(table(missing$hba1c_val))

prop.table(table(missing$no_smoking_cat))

prop.table(table(missing$bmi_val))

missing %>% filter(bmi_val=="over") %>% summarise(median=median(prebmi))
missing %>% filter(hba1c_val=="over") %>% summarise(median=median(prehba1c))


(missing %>% filter(!is.na(ckdpc_egfr60_risk_confirmed_score)) %>% count()) / (cohort %>% count())

(missing %>% filter(!is.na(ckdpc_40egfr_risk_score) & !is.na(preegfr)) %>% count()) / (cohort %>% count())
  

ckd60_cohort <- cohort %>%
  filter(!is.na(ckdpc_egfr60_risk_confirmed_score)) %>%
  select(-c(contains("40"), contains("risk_total")))

ckd40_cohort <- cohort %>%
  filter(!is.na(ckdpc_40egfr_risk_score)) %>%
  select(-c(contains("60"), contains("ckd_345")))


############################################################################################

# 2 How well do uncalibrated risk scores predict incidence?


# CKD <60

ckd60_cohort$risk_decile <- ntile(ckd60_cohort$ckdpc_egfr60_risk_confirmed_score, 10)

## Get mean predicted probabilities by studydrug
predicted <- ckd60_cohort %>%
  group_by(risk_decile, studydrug) %>%
  summarise(mean_ckd60_pred=mean(ckdpc_egfr60_risk_confirmed_score)/100)

## Find actual observed probabilities by risk score category and studydrug
observed_dpp4su <- survfit(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ risk_decile, data=ckd60_cohort[ckd60_cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_dpp4su=1-estimate,
         lower_ci_dpp4su=1-conf.low,
         upper_ci_dpp4su=1-conf.high) %>%
  select(observed_dpp4su, lower_ci_dpp4su, upper_ci_dpp4su, strata)

dpp4su_events <- ckd60_cohort %>%
  filter(studydrug=="DPP4SU" & ckd_345_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4SU=n())

observed_sglt2 <- survfit(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ risk_decile, data=ckd60_cohort[ckd60_cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_sglt2=1-estimate,
         lower_ci_sglt2=1-conf.low,
         upper_ci_sglt2=1-conf.high) %>%
  select(observed_sglt2, lower_ci_sglt2, upper_ci_sglt2, strata)

sglt2_events <- ckd60_cohort %>%
  filter(studydrug=="SGLT2" & ckd_345_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLT2=n())

obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug=="DPP4SU")), observed_dpp4su),
  cbind((predicted %>% filter(studydrug=="SGLT2")), observed_sglt2)
) %>%
  mutate(observed=coalesce(observed_dpp4su, observed_sglt2),
         lower_ci=coalesce(lower_ci_dpp4su, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_dpp4su, upper_ci_sglt2))

events_table <- data.frame(t(dpp4su_events %>%
  inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd60_pred=NA, risk_decile=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=risk_decile, group=studydrug, color=studydrug)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_ckd60_pred*100), position=dodge, shape=17, size=2) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,60,by=10)), limits=c(0,61)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Uncalibrated risk score vs CKD incidence (5 year)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))

## C-stat
ckd60_cohort <- ckd60_cohort %>%
  mutate(ckdpc_egfr60_risk_confirmed_survival=(100-ckdpc_egfr60_risk_confirmed_score)/100)

surv_mod <- coxph(Surv(ckd_345_censtime_yrs, ckd_345_censvar)~ckdpc_egfr60_risk_confirmed_survival, data=ckd60_cohort, method="breslow")
summary(surv_mod)$concordance[1]
summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2])
summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2])
# 0.9035459 (0.897121-0.9099707)



# 40% decline in eGFR / CKD stage 5

ckd40_cohort$risk_decile <- ntile(ckd40_cohort$ckdpc_40egfr_risk_score, 10)

## Get mean predicted probabilities by studydrug
predicted <- ckd40_cohort %>%
  group_by(risk_decile, studydrug) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_risk_score)/100)

## Find actual observed probabilities by risk score category and studydrug
observed_dpp4su <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, data=ckd40_cohort[ckd40_cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_dpp4su=1-estimate,
         lower_ci_dpp4su=1-conf.low,
         upper_ci_dpp4su=1-conf.high) %>%
  select(observed_dpp4su, lower_ci_dpp4su, upper_ci_dpp4su, strata)

dpp4su_events <- ckd40_cohort %>%
  filter(studydrug=="DPP4SU" & ckd_egfr40_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4SU=n())

observed_sglt2 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, data=ckd40_cohort[ckd40_cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_sglt2=1-estimate,
         lower_ci_sglt2=1-conf.low,
         upper_ci_sglt2=1-conf.high) %>%
  select(observed_sglt2, lower_ci_sglt2, upper_ci_sglt2, strata)

sglt2_events <- ckd40_cohort %>%
  filter(studydrug=="SGLT2" & ckd_egfr40_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLT2=n())

obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug=="DPP4SU")), observed_dpp4su),
  cbind((predicted %>% filter(studydrug=="SGLT2")), observed_sglt2)
) %>%
  mutate(observed=coalesce(observed_dpp4su, observed_sglt2),
         lower_ci=coalesce(lower_ci_dpp4su, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_dpp4su, upper_ci_sglt2))

events_table <- data.frame(t(dpp4su_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd40_pred=NA, risk_decile=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=risk_decile, group=studydrug, color=studydrug)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_ckd40_pred*100), position=dodge, shape=17, size=2) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,7,by=1)), limits=c(0,7.5)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Uncalibrated risk score vs incidence (3 year)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))

## C-stat
ckd40_cohort <- ckd40_cohort %>%
  mutate(ckdpc_40egfr_risk_survival=(100-ckdpc_40egfr_risk_score)/100)

surv_mod <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckdpc_40egfr_risk_survival, data=ckd40_cohort, method="breslow")
summary(surv_mod)$concordance[1]
summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2])
summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2])
# 0.7091707 (0.6762553-0.7420862)


############################################################################################

# 3 Does calibrating help?

# CKD < 60

## Assign random 20% of DPP4SU arm as calibration cohort and remove from main cohort
set.seed(123)

cal_ckd60_cohort <- ckd60_cohort %>%
  filter(studydrug=="DPP4SU") %>%
  slice_sample(prop=0.2)
#8,486

noncal_ckd60_cohort <- ckd60_cohort %>%
  anti_join(cal_ckd60_cohort, by=c("patid", "dstartdate", "studydrug"))
table(noncal_ckd60_cohort$studydrug)
#SGLT2: 22,651
#DPP4SU: 33948

## Males and females together
## Original survival constant: 0.9212477
## BUT this is not baseline hazard due to different equation: BH=exp(-5^0.9212477)
exp(-5^0.9212477)
#0.01221876
recal_mod <- coxph(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ offset(ckdpc_egfr60_risk_confirmed_lin_predictor), data=cal_ckd60_cohort)
new_ckd60_surv <- summary(survfit(recal_mod),time=5)$surv
sprintf("%.15f", new_ckd60_surv)
# 0.955351637598774
## new surv constant=(log(-log(a)))/log(5)
new_ckd60_surv2 <- (log(-log(new_ckd60_surv)))/log(5)
new_ckd60_surv2


## Recalculate scores in rest of cohort
noncal_ckd60_cohort <- noncal_ckd60_cohort %>%
  mutate(centred_ckd60_lin_predictor=ckdpc_egfr60_risk_confirmed_lin_predictor-mean(ckdpc_egfr60_risk_confirmed_lin_predictor)) %>%
  ungroup() %>%
  mutate(ckdpc_egfr60_risk_confirmed_score_cal=(1-(new_ckd60_surv^exp(centred_ckd60_lin_predictor)))*100)


## Plot
noncal_ckd60_cohort$risk_decile <- ntile(noncal_ckd60_cohort$ckdpc_egfr60_risk_confirmed_score_cal, 10)

### Get mean predicted probabilities by studydrug
predicted <- noncal_ckd60_cohort %>%
  group_by(risk_decile, studydrug) %>%
  summarise(mean_ckd60_pred=mean(ckdpc_egfr60_risk_confirmed_score_cal)/100)

### Find actual observed probabilities by risk score category and studydrug
observed_dpp4su <- survfit(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ risk_decile, data=noncal_ckd60_cohort[noncal_ckd60_cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_dpp4su=1-estimate,
         lower_ci_dpp4su=1-conf.low,
         upper_ci_dpp4su=1-conf.high) %>%
  select(observed_dpp4su, lower_ci_dpp4su, upper_ci_dpp4su, strata)

dpp4su_events <- noncal_ckd60_cohort %>%
  filter(studydrug=="DPP4SU" & ckd_345_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4SU=n())

observed_sglt2 <- survfit(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ risk_decile, data=noncal_ckd60_cohort[noncal_ckd60_cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_sglt2=1-estimate,
         lower_ci_sglt2=1-conf.low,
         upper_ci_sglt2=1-conf.high) %>%
  select(observed_sglt2, lower_ci_sglt2, upper_ci_sglt2, strata)

sglt2_events <- noncal_ckd60_cohort %>%
  filter(studydrug=="SGLT2" & ckd_345_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLT2=n())

obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug=="DPP4SU")), observed_dpp4su),
  cbind((predicted %>% filter(studydrug=="SGLT2")), observed_sglt2)
) %>%
  mutate(observed=coalesce(observed_dpp4su, observed_sglt2),
         lower_ci=coalesce(lower_ci_dpp4su, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_dpp4su, upper_ci_sglt2))

events_table <- data.frame(t(dpp4su_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd60_pred=NA, risk_decile=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=risk_decile, group=studydrug, color=studydrug)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_ckd60_pred*100), position=dodge, shape=17, size=2) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,60,by=10)), limits=c(0,61)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Calibrated risk score vs CKD incidence (5 year)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))




# 40% decline in eGFR / CKD 5

## Assign random 20% of DPP4SU arm as calibration cohort and remove from main cohort
set.seed(123)

cal_ckd40_cohort <- ckd40_cohort %>%
  filter(studydrug=="DPP4SU") %>%
  slice_sample(prop=0.2)
#8,522

noncal_ckd40_cohort <- ckd40_cohort %>%
  anti_join(cal_ckd40_cohort, by=c("patid", "dstartdate", "studydrug"))
table(noncal_ckd40_cohort$studydrug)
#SGLT2: 22,777
#DPP4SU: 34,091



# First need the linear predictor without the intercept
lp_noint=lin_pred2 + 4.182049

# Fit logistic model using lp_noint as offset (same as forcing beta for lp_noint to be 1)
new_mod <- glm(DAY30~offset(lp_noint),family="binomial")
new_mod$coef









## Males and females together
## Original: 
recal_mod <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ offset(ckdpc_40egfr_risk_lin_predictor), data=cal_ckd40_cohort)
new_ckd40_surv <- summary(survfit(recal_mod),time=3)$surv
sprintf("%.15f", new_ckd40_surv)
# 0.993905745768886

## Recalculate scores in rest of cohort
noncal_ckd40_cohort <- noncal_ckd40_cohort %>%
  mutate(centred_40ckd_lin_predictor=ckdpc_40egfr_risk_lin_predictor-mean(ckdpc_40egfr_risk_lin_predictor)) %>%
  ungroup() %>%
  mutate(ckdpc_40egfr_risk_score_cal=(1-(new_ckd40_surv^exp(centred_40ckd_lin_predictor)))*100)


## Plot
noncal_ckd40_cohort$risk_decile <- ntile(noncal_ckd40_cohort$ckdpc_40egfr_risk_score_cal, 10)

### Get mean predicted probabilities by studydrug
predicted <- noncal_ckd40_cohort %>%
  group_by(risk_decile, studydrug) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_risk_score_cal)/100)

### Find actual observed probabilities by risk score category and studydrug
observed_dpp4su <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, data=noncal_ckd40_cohort[noncal_ckd40_cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_dpp4su=1-estimate,
         lower_ci_dpp4su=1-conf.low,
         upper_ci_dpp4su=1-conf.high) %>%
  select(observed_dpp4su, lower_ci_dpp4su, upper_ci_dpp4su, strata)

dpp4su_events <- noncal_ckd40_cohort %>%
  filter(studydrug=="DPP4SU" & ckd_345_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4SU=n())

observed_sglt2 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, data=noncal_ckd40_cohort[noncal_ckd40_cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_sglt2=1-estimate,
         lower_ci_sglt2=1-conf.low,
         upper_ci_sglt2=1-conf.high) %>%
  select(observed_sglt2, lower_ci_sglt2, upper_ci_sglt2, strata)

sglt2_events <- noncal_ckd40_cohort %>%
  filter(studydrug=="SGLT2" & ckd_345_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLT2=n())

obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug=="DPP4SU")), observed_dpp4su),
  cbind((predicted %>% filter(studydrug=="SGLT2")), observed_sglt2)
) %>%
  mutate(observed=coalesce(observed_dpp4su, observed_sglt2),
         lower_ci=coalesce(lower_ci_dpp4su, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_dpp4su, upper_ci_sglt2))

events_table <- data.frame(t(dpp4su_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd40_pred=NA, risk_decile=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=risk_decile, group=studydrug, color=studydrug)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_ckd40_pred*100), position=dodge, shape=17, size=2) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,7,by=1)), limits=c(0,7.5)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Calibrated risk score vs CKD incidence (5 year)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))











############################################################################################

# 4 Look at survival benefit from SGLT2s when add hazard ratio from trials meta-analysis
# https://jamanetwork.com/journals/jamacardiology/fullarticle/2771459

# Add HR from trials
noncal_cohort <- noncal_cohort %>%
  mutate(qrisk2_survival_cal=(100-qrisk2_5yr_score_cal)/100,
         qrisk2_survival_cal_sglt2=qrisk2_survival_cal^0.63,
         qrisk2_5yr_score_cal_sglt2=100-(qrisk2_survival_cal_sglt2*100),
         
         qdhf_survival_cal=(100-qdiabeteshf_5yr_score_cal)/100,
         qdhf_survival_cal_sglt2=qdhf_survival_cal^0.63,
         qdiabeteshf_5yr_score_cal_sglt2=100-(qdhf_survival_cal_sglt2*100))


# Distribution of predicted benefits
noncal_cohort <- noncal_cohort %>%
  mutate(qrisk2_sglt2_benefit=qrisk2_survival_cal_sglt2-qrisk2_survival_cal,
         qdhf_sglt2_benefit=qdhf_survival_cal_sglt2-qdhf_survival_cal)

describe(noncal_cohort$qrisk2_sglt2_benefit)
#mean=1.8%

describe(noncal_cohort$qdhf_sglt2_benefit)
#mean=1.7%

test <- noncal_cohort %>% select(qrisk2_sglt2_benefit, qdhf_sglt2_benefit) %>% pivot_longer(cols=c(qrisk2_sglt2_benefit, qdhf_sglt2_benefit))

ggplot(test, aes(x=value*100, color=name, fill=name)) + 
  geom_histogram(binwidth=0.05, alpha=0.5, position="identity") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))



# Compare predicted and actual benefits by decile (predicted = mean of decile, as previous)

## Unadjusted

## QRISK2

noncal_cohort$sglt2_benefit_decile <- ntile(noncal_cohort$qrisk2_sglt2_benefit, 10)
noncal_cohort <- noncal_cohort %>% mutate(sglt2_benefit_decile=as.factor(sglt2_benefit_decile))

pred <- noncal_cohort %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(mean_qrisk2_predicted_benefit=mean(qrisk2_sglt2_benefit, na.rm=T))


noncal_cohort <- noncal_cohort %>% mutate(studydrug=as.vector(studydrug))
ddist <- datadist(noncal_cohort)

model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug*sglt2_benefit_decile, data=noncal_cohort, x=TRUE, y=TRUE, surv=TRUE)

survival_est <- survest(model, newdata=expand.grid(studydrug=c("SGLT2","DPP4SU"), sglt2_benefit_decile=c(1:10)), times=5)


obs <- data.frame(surv=unlist(survival_est$surv), studydrug=rep(c("SGLT2","DPP4SU"),5), sglt2_benefit_decile=rep(1:10, rep_len(2, 10)))

obs <- obs %>%
  pivot_wider(id_cols=sglt2_benefit_decile, values_from=surv, names_from=studydrug) %>%
  mutate(surv_diff=SGLT2-DPP4SU) %>%
  select(sglt2_benefit_decile, surv_diff)


se <- data.frame(se=unlist(survival_est$std.err), studydrug=rep(c("SGLT2","DPP4SU"),5), sglt2_benefit_decile=rep(1:10, rep_len(2, 10)))

se <- se %>%
  pivot_wider(id_cols=sglt2_benefit_decile, values_from=se, names_from=studydrug) %>%
  mutate(se=sqrt((SGLT2^2)+(DPP4SU^2))) %>%
  select(sglt2_benefit_decile, se)

obs <- obs %>%
  inner_join(se, by="sglt2_benefit_decile") %>%
  mutate(lower_ci=surv_diff-(1.96*se),
         upper_ci=surv_diff+(1.96*se))

dpp4su_events <- noncal_cohort %>%
  filter(studydrug=="DPP4SU" & mace_censvar==1) %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(DPP4SU=n())

sglt2_events <- noncal_cohort %>%
  filter(studydrug=="SGLT2" & mace_censvar==1) %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(SGLT2=n())

events_table <- data.frame(t(dpp4su_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="sglt2_benefit_decile")

obs <- obs %>% mutate(sglt2_benefit_decile=as.factor(sglt2_benefit_decile))
obs_v_pred <- pred %>% inner_join(obs, by="sglt2_benefit_decile")


empty_tick <- obs_v_pred %>%
  filter(sglt2_benefit_decile==1) %>%
  mutate(mean_predicted_benefit=NA, surv_diff=NA, lower_ci=NA, upper_ci=NA, sglt2_benefit_decile=as.factor(0))

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_qrisk2_predicted_benefit*100)) + 
  geom_point(aes(y = surv_diff*100)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1) +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  xlab("Predicted SGLT2 benefit (mean by predicted decile)") + ylab("Observed SGLT2 benefit (by predicted decile)") +
  scale_x_continuous(breaks=c(seq(0,5,by=1)), limits=c(0,5)) +
  scale_y_continuous(breaks=c(seq(-2,6,by=1)), limits=c(-2,6)) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5)))

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))

#NNT
obs_v_pred <- obs_v_pred %>%
  mutate(nnt_predicted=1/(mean_qdhf_predicted_benefit*100),
         nnt_observed=1/(surv_diff*100))
obs_v_pred



## QDHF

noncal_cohort$sglt2_benefit_decile <- ntile(noncal_cohort$qdhf_sglt2_benefit, 10)
noncal_cohort <- noncal_cohort %>% mutate(sglt2_benefit_decile=as.factor(sglt2_benefit_decile))

pred <- noncal_cohort %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(mean_qdhf_predicted_benefit=mean(qdhf_sglt2_benefit, na.rm=T))


noncal_cohort <- noncal_cohort %>% mutate(studydrug=as.vector(studydrug))
ddist <- datadist(noncal_cohort)

model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug*sglt2_benefit_decile, data=noncal_cohort, x=TRUE, y=TRUE, surv=TRUE)

survival_est <- survest(model, newdata=expand.grid(studydrug=c("SGLT2","DPP4SU"), sglt2_benefit_decile=c(1:10)), times=5)


obs <- data.frame(surv=unlist(survival_est$surv), studydrug=rep(c("SGLT2","DPP4SU"),5), sglt2_benefit_decile=rep(1:10, rep_len(2, 10)))

obs <- obs %>%
  pivot_wider(id_cols=sglt2_benefit_decile, values_from=surv, names_from=studydrug) %>%
  mutate(surv_diff=SGLT2-DPP4SU) %>%
  select(sglt2_benefit_decile, surv_diff)


se <- data.frame(se=unlist(survival_est$std.err), studydrug=rep(c("SGLT2","DPP4SU"),5), sglt2_benefit_decile=rep(1:10, rep_len(2, 10)))

se <- se %>%
  pivot_wider(id_cols=sglt2_benefit_decile, values_from=se, names_from=studydrug) %>%
  mutate(se=sqrt((SGLT2^2)+(DPP4SU^2))) %>%
  select(sglt2_benefit_decile, se)

obs <- obs %>%
  inner_join(se, by="sglt2_benefit_decile") %>%
  mutate(lower_ci=surv_diff-(1.96*se),
         upper_ci=surv_diff+(1.96*se))

dpp4su_events <- noncal_cohort %>%
  filter(studydrug=="DPP4SU" & mace_censvar==1) %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(DPP4SU=n())

sglt2_events <- noncal_cohort %>%
  filter(studydrug=="SGLT2" & mace_censvar==1) %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(SGLT2=n())

events_table <- data.frame(t(dpp4su_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="sglt2_benefit_decile")

obs <- obs %>% mutate(sglt2_benefit_decile=as.factor(sglt2_benefit_decile))
obs_v_pred <- pred %>% inner_join(obs, by="sglt2_benefit_decile")


empty_tick <- obs_v_pred %>%
  filter(sglt2_benefit_decile==1) %>%
  mutate(mean_predicted_benefit=NA, surv_diff=NA, lower_ci=NA, upper_ci=NA, sglt2_benefit_decile=as.factor(0))

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_qdhf_predicted_benefit*100)) + 
  geom_point(aes(y = surv_diff*100)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1) +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  xlab("Predicted SGLT2 benefit (mean by predicted decile)") + ylab("Observed SGLT2 benefit (by predicted decile)") +
  scale_x_continuous(breaks=c(seq(0,6,by=1)), limits=c(0,6)) +
  scale_y_continuous(breaks=c(seq(-2,8,by=1)), limits=c(-2,8)) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5)))

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))

#NNT
obs_v_pred <- obs_v_pred %>%
  mutate(nnt_predicted=1/(mean_qdhf_predicted_benefit*100),
         nnt_observed=1/(surv_diff*100))
obs_v_pred




## Observed adjusted for age, sex, duration, IMD, drugline and ncurrtx

## QRISK2

memory.limit(size=40000)
gc()

noncal_cohort$sglt2_benefit_decile <- ntile(noncal_cohort$qrisk2_sglt2_benefit, 10)
noncal_cohort <- noncal_cohort %>% mutate(sglt2_benefit_decile=as.factor(sglt2_benefit_decile))

# model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + drugline_all + ncurrtx, data=noncal_cohort, x=TRUE, y=TRUE, surv=TRUE)
# 
# obs_SGLT2 <- noncal_cohort %>%
#   select(patid, sglt2_benefit_decile, dstartdate_age, malesex, dstartdate_dm_dur_all, imd2015_10, drugline_all, ncurrtx) %>%
#   mutate(studydrug="SGLT2",
#          rowno=row_number())
# 
# observed_sglt2 <- survfit(model, newdata=as.data.frame(obs_SGLT2)) %>%
#   tidy() %>%
#   filter(time==5) %>%
#   pivot_longer(cols=-c(time, n.risk, n.event, n.censor), names_to = c(".value", "group"), names_pattern = "(.*)\\.(.*)") %>%
#   select(group, estimate_sglt2=estimate) %>%
#   mutate(group=as.numeric(group)) %>%
#   inner_join(obs_SGLT2, by=c("group"="rowno"))
# 
# save(observed_sglt2, file="../../Raw data/HF_adjust_SGLT2.Rda")
# 
# obs_DPP4SU <- noncal_cohort %>%
#   select(patid, sglt2_benefit_decile, dstartdate_age, malesex, dstartdate_dm_dur_all, imd2015_10, drugline_all, ncurrtx) %>%
#   mutate(studydrug="DPP4SU",
#          rowno=row_number())
# 
# observed_dpp4su <- survfit(model, newdata=as.data.frame(obs_DPP4SU)) %>%
#   tidy() %>%
#   filter(time==5) %>%
#   pivot_longer(cols=-c(time, n.risk, n.event, n.censor), names_to = c(".value", "group"), names_pattern = "(.*)\\.(.*)") %>%
#   select(group, estimate_dpp4su=estimate) %>%
#   mutate(group=as.numeric(group)) %>%
#   inner_join(obs_DPP4SU, by=c("group"="rowno"))
# 
# save(observed_dpp4su, file="../../Raw data/HF_adjust_DPP4SU.Rda")

load("../../Raw data/HF_adjust_SGLT2.Rda")
load("../../Raw data/HF_adjust_DPP4SU.Rda")

observed <- observed_sglt2 %>%
  select(group, estimate_sglt2) %>%
  inner_join(observed_dpp4su, by="group") %>%
  mutate(survdiff=estimate_sglt2-estimate_dpp4su) %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(mean_benefit=mean(survdiff),
            median_benefit=median(survdiff),
            lq_benefit=quantile(survdiff, prob=c(.25)),
            uq_benefit=quantile(survdiff, prob=c(.75)))

### Predicted same as previous
pred <- noncal_cohort %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(median_qrisk2_predicted_benefit=median(qrisk2_sglt2_benefit, na.rm=T))


obs_v_pred <- pred %>% inner_join(observed, by="sglt2_benefit_decile")


empty_tick <- obs_v_pred %>%
  filter(sglt2_benefit_decile==1) %>%
  mutate(mean_benefit=NA, median_benefit=NA, lq_benefit=NA, uq_benefit=NA, mean_qrisk2_predicted_benefit=NA, sglt2_benefit_decile=as.factor(0))

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=median_qrisk2_predicted_benefit*100)) + 
  geom_point(aes(y = median_benefit*100)) +
  geom_errorbar(aes(ymax=uq_benefit*100,ymin=lq_benefit*100),alpha=1) +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  xlab("Predicted SGLT2 benefit (mean by predicted decile)") + ylab("Observed SGLT2 benefit (by predicted decile)") +
  scale_x_continuous(breaks=c(seq(0,6,by=1)), limits=c(0,6)) +
  scale_y_continuous(breaks=c(seq(-2,6,by=1)), limits=c(-2,8)) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5)))

p1





noncal_cohort$sglt2_benefit_decile <- ntile(noncal_cohort$qrisk2_sglt2_benefit, 10)
noncal_cohort <- noncal_cohort %>% mutate(sglt2_benefit_decile=as.factor(sglt2_benefit_decile))

pred <- noncal_cohort %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(mean_qrisk2_predicted_benefit=mean(qrisk2_sglt2_benefit, na.rm=T))


noncal_cohort <- noncal_cohort %>% mutate(studydrug=as.vector(studydrug))
ddist <- datadist(noncal_cohort)

model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + drugline_all + ncurrtx + qrisk2_5yr_score, data=noncal_cohort, x=TRUE, y=TRUE, surv=TRUE)

obs_SGLT2 <- noncal_cohort %>%
    select(patid, sglt2_benefit_decile, dstartdate_age, malesex, dstartdate_dm_dur_all, imd2015_10, drugline_all, ncurrtx, qrisk2_5yr_score) %>%
    mutate(studydrug="SGLT2",
           rowno=row_number())


memory.limit(size=40000)
gc()
survival_est2 <- survest(model, newdata=as.data.frame(obs_SGLT2), times=5)







