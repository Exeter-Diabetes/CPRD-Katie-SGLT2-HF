
# HRs for SGLT2 vs DPP4SU for:
## 3-point MACE (ischaemic stroke, MI, CV death)
## Expanded MACE (ischaemic stroke, MI, CV death, revasc, HES unstable angina)
## HF
## all-cause hospitalisation
## all-cause mortality

# Unadjusted and adjusted for age + sex + duration + IMD + QRisk(5yr) + drugline + ncurrtx (adding ethnicity makes ~no difference)

# Not including GLP1 currently
 
############################################################################################

# Setup
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(flextable)
library(PSweight)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# 1 Cohort selection and variable setup

## A Cohort selection (see cohort_definition function for details)

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Raw data/")
load("20231108_t2d_1stinstance.Rda")
load("20231108_t2d_all_drug_periods.Rda")

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Scripts/Functions")
source("cohort_definition.R")

cohort <- define_cohort(t2d_1stinstance, t2d_all_drug_periods)

cohort <- cohort %>% filter(drugclass!="GLP1") %>% mutate(studydrug=relevel(factor(ifelse(drugclass=="SGLT2", "SGLT2", "DPP4SU")), ref="DPP4SU"))
#cohort <- cohort %>% filter(drugclass!="GLP1" & drugclass!="DPP4") %>% mutate(studydrug=relevel(factor(drugclass), ref="SU"))
#cohort <- cohort %>% filter(drugclass!="GLP1" & drugclass!="SGLT2") %>% mutate(studydrug=relevel(factor(drugclass), ref="SU"))

table(cohort$studydrug)
# DPP4SU 98397
# SGLT2 45781

# SU 36378
# SGLT2 45781

# SU 36378
# DPP4 58964 


## B Make variables for survival analysis of all endpoints (see survival_variables function for details)

source("survival_variables.R")

cohort <- add_surv_vars(cohort, main_only=TRUE)
         

## C Just keep variables of interest

cohort <- cohort %>%
  
  select(patid, malesex, ethnicity_5cat_decoded, imd2015_10, regstartdate, gp_record_end, death_date, drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, prehba1c, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preegfr, preckdstage, contains("cens"), starts_with("ckdpc"), qrisk2_5yr_score, initiation_year, hosp_admission_prev_year, predrug_hypertension, qrisk2_smoking_cat, drugorder)

rm(list=setdiff(ls(), "cohort"))


############################################################################################

# 2 Look at hazard ratios

test <- cohort %>%
  select(patid, ethnicity_5cat_decoded, imd2015_10, qrisk2_5yr_score, drugline_all, ncurrtx, initiation_year, prebmi, prehba1c, presbp, qrisk2_smoking_cat, predrug_hypertension, hosp_admission_prev_year) %>%
  filter(is.na(ethnicity_5cat_decoded) | is.na(imd2015_10) | is.na(qrisk2_5yr_score) | is.na(drugline_all) | is.na(ncurrtx) | is.na(initiation_year) | is.na(prebmi) | is.na(prehba1c) | is.na(presbp) | is.na(qrisk2_smoking_cat) | is.na(predrug_hypertension) | is.na(hosp_admission_prev_year))

cohort <- data.frame(cohort) %>% filter(!is.na(presbp))

ps.formula <- formula(paste("studydrug ~ dstartdate_age + ethnicity_5cat_decoded + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx + initiation_year + prebmi + prehba1c + presbp + qrisk2_smoking_cat + predrug_hypertension + hosp_admission_prev_year + drugorder"))

overlap <- SumStat(ps.formula=ps.formula, data = cohort, weight = c("IPW", "overlap"))

cohort2$overlap_weights <- overlap$ps.weights$overlap
cohort2$iptw_weights <- overlap$ps.weights$IPW

f <- as.formula(paste("Surv(hf_censtime_yrs, hf_censvar) ~  studydrug"))

unadjusted <- coxph(f, cohort2, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  #filter(term=="studydrugDPP4") %>%
  mutate(unadjusted_HR=paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>%
  select(unadjusted_HR)










main_outcomes <- c("mace", "hf", "hosp", "death")


all_hrs <- data.frame()

for (i in main_outcomes) {
  
  censvar_var=paste0(i, "_censvar")
  censtime_var=paste0(i, "_censtime_yrs")
  
  count <- cohort %>%
    group_by(studydrug) %>%
    summarise(count=n()) %>%
    pivot_wider(names_from=studydrug,
                names_glue="{studydrug}_count",
                values_from=count)
  
  followup <- cohort %>%
    group_by(studydrug) %>%
    summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
    pivot_wider(names_from=studydrug,
                names_glue="{studydrug}_followup",
                values_from=time)
  
  events <- cohort %>%
    group_by(studydrug) %>%
    summarise(event_count=sum(!!sym(censvar_var)),
              drug_count=n()) %>%
    mutate(events_perc=round(event_count*100/drug_count, 1),
           events=paste0(event_count, " (", events_perc, "%)")) %>%
    select(studydrug, events) %>%
    pivot_wider(names_from=studydrug,
                names_glue="{studydrug}_events",
                values_from=events)
  
  
  f <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug"))
  
  unadjusted <- coxph(f, cohort) %>%
    tidy(conf.int=TRUE, exponentiate=TRUE) %>%
    filter(term=="studydrugSGLT2") %>%
    #filter(term=="studydrugDPP4") %>%
    mutate(unadjusted_HR=paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>%
    select(unadjusted_HR)
  
  f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug + dstartdate_age + ethnicity_5cat_decoded + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx + initiation_year + prebmi + prehba1c + presbp + qrisk2_smoking_cat + predrug_hypertension + hosp_admission_prev_year + drugorder"))
  
  adjusted <- coxph(f_adjusted, cohort) %>%
    tidy(conf.int=TRUE, exponentiate=TRUE) %>%
    filter(term=="studydrugSGLT2") %>%
    #filter(term=="studydrugDPP4") %>%
    mutate(adjusted_HR=paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>%
    select(adjusted_HR)
  
  outcome_hr <- cbind(outcome=i, count, followup, events, unadjusted, adjusted)
  
  all_hrs <- rbind(all_hrs, outcome_hr)
  
}



flextable(all_hrs)


