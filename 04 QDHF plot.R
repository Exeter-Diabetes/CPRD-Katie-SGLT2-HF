
# Plot QDHF vs HF incidence for control group, and KM curves for control vs SGLT2


############################################################################################

# Setup
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(gridExtra)
library(grid)
library(PSweight)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# 1 Cohort selection and variable setup

setwd("/slade/CPRD_data/Katie SGLT2/Processed data/")
load("treatment_outcome_cohort_jun24.rda")
#169,041

table(cohort$studydrug)
# DPP4SU 111673   
# SGLT2 57368  


## B Make variables for survival analysis of all endpoints (see survival_variables function for details)

setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions/")
source("survival_variables.R")

cohort <- add_surv_vars(cohort, main_only=TRUE)
         

## C Just keep variables of interest

cohort <- cohort %>%
  
  select(patid, malesex, ethnicity_qrisk2_decoded, imd2015_10, regstartdate, gp_record_end, death_date, drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx_cat, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, prehba1c2yrs, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preegfr, preckdstage, contains("cens"), qdiabeteshf_5yr_score, qdiabeteshf_lin_predictor, initiation_year, hosp_admission_prev_year_count, hypertension, qrisk2_smoking_cat, drugorder, qrisk2_5yr_score, qrisk2_10yr_score, qrisk2_lin_predictor, predrug_af, tds_2011, uacr)

rm(list=setdiff(ls(), "cohort"))


############################################################################################

# 2 Compare QDHF to actual HF incidence (by QDHF decile) in DPP4SU arm

dpp4su <- cohort %>%
  filter(studydrug=="DPP4SU") %>%
  mutate(qdhf_decile=ntile(qdiabeteshf_5yr_score, 10))

predicted <- dpp4su %>%
  group_by(qdhf_decile) %>%
  summarise(median_qdhf_pred=median(qdiabeteshf_5yr_score)/100,
            mean_qdhf_pred=mean(qdiabeteshf_5yr_score)/100)

observed <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ qdhf_decile, data=dpp4su) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed=1-estimate,
         lower_ci=1-conf.low,
         upper_ci=1-conf.high)



events_table <- t(dpp4su %>%
  filter(hf_censvar==1) %>%
  group_by(qdhf_decile) %>%
  summarise(events=n()))


obs_v_pred <- cbind(predicted, observed) %>%
  select(qdhf_decile, median_qdhf_pred, observed, lower_ci, upper_ci)

dodge <- position_dodge(width=0.3)

tiff("../../Plots/QDHF_median.tiff", width=9, height=7, units = "in", res=600) 

ggplot(data=obs_v_pred, aes(x=qdhf_decile)) +
  geom_point(aes(y = observed*100, color="observed", shape="observed"), size=3, position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100, ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = median_qdhf_pred*100, color="median_qdhf_pred", shape="median_qdhf_pred"), position=dodge, stroke=1.5, size=4) +
  theme_bw() +
  xlab("Predicted 5-year absolute HF risk decile") + ylab("5-year HF risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,14,by=2)), limits=c(0,14)) +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"),
        axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold"),
        plot.margin = margin(10,10,10,10),
        axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        legend.position=c(0.4, 0.9),
        legend.text=element_text(size=16)) +
  scale_color_manual(values = c(observed = "black", median_qdhf_pred = "black"), labels = c(observed = "Observed (with 95% CI)", median_qdhf_pred = "Median predicted risk per decile"), name="") +
  scale_shape_manual(values = c(observed = 16, median_qdhf_pred = 4), labels = c(observed = "Observed (with 95% CI)", median_qdhf_pred = "Median predicted risk per decile"), name="") #+
#ggtitle("5-year QDiabetes-Heart Failure vs heart failure incidence")

dev.off()  



## C-score

dpp4su <- dpp4su %>% mutate(qdhf_survival=(100-qdiabeteshf_5yr_score)/100)

surv_mod <- coxph(Surv(hf_censtime_yrs, hf_censvar)~qdhf_survival, data=dpp4su, method="breslow")
cstat <- round(summary(surv_mod)$concordance[1], 3)
cstat_lower <- round(summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2]), 3)
cstat_upper <- round(summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2]), 3)
#0.72 (0.70, 0.73)


## Alternative method

library(SurvMetrics)

Cindex(Surv(dpp4su$hf_censtime_yrs, dpp4su$hf_censvar), dpp4su$qdhf_survival)
#just ran for ages...


## Brier score

library(riskRegression)

raw_mod <- coxph(Surv(hf_censtime_yrs, hf_censvar) ~ qdhf_survival, 
                 data=dpp4su, x=T)


## At time 5 years
score_raw <- 
  Score(object = list(raw_mod),
        formula = Surv(hf_censtime_yrs, hf_censvar) ~ 1, # null model
        data = dpp4su,
       metrics="brier",
       times=5) 

score_raw$Brier$score$Brier[2]  #[1] is NULL model
#0.036
score_raw$Brier$score$lower[2]
#0.034
score_raw$Brier$score$upper[2]
#0.038



## Integrated Brier score for all time points

score_raw <- 
  Score(object = list(raw_mod),
        formula = Surv(hf_censtime_yrs, hf_censvar) ~ 1, # null model
        data = dpp4su,
        summary="ibs") 

score_raw$Brier$score$Brier[2]
#0.017
score_raw$Brier$score$lower[2]
#0.016
score_raw$Brier$score$upper[2]
#0.018
