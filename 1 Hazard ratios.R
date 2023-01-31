
# HRs for SGLT2 vs DPP4SU for:
## 3-point MACE (ischaemic stroke, MI, CV death)
## Expanded MACE (ischaemic stroke, MI, CV death, revasc, HES unstable angina)
## HF
## CKD: onset of CKD stage 3a-5
## CKD: fall in eGFR of <=40% from baseline or onset of CKD stage 5
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

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# 1 Cohort selection and variable setup

## A Cohort selection (see cohort_definition function for details)

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Raw data/")
load("20230116_t2d_1stinstance.Rda")
load("20230116_t2d_all_drug_periods.Rda")

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Scripts/Functions")
source("cohort_definition.R")

cohort <- define_cohort(t2d_1stinstance, t2d_all_drug_periods)

cohort <- cohort %>% filter(studydrug!="GLP1")

table(cohort$studydrug)
# DPP4SU 90848
# SGLT2 48280


## B Make variables for survival analysis of all endpoints (see survival_variables function for details)

source("survival_variables.R")

cohort <- cohort %>%
  mutate(postdrug_first_ischaemic_stroke=postdrug_first_stroke,
         postdrug_first_unstable_angina=postdrug_first_angina)

cohort <- add_surv_vars(cohort, main_only=TRUE)
         

## C Just keep variables of interest

cohort <- cohort %>%
  
  select(patid, malesex, ethnicity_5cat_decoded, imd2015_10, regstartdate, gp_record_end, death_date, drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, prehba1c, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preckdstage, contains("cens"), qrisk2_5yr_score, last_sglt2_stop)

rm(list=setdiff(ls(), "cohort"))


############################################################################################

# 2 Look at hazard ratios

main_outcomes <- c("mace", "expanded_mace", "hf", "ckd_345", "ckd_egfr40", "hosp", "death")


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
    mutate(unadjusted_HR=paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>%
    select(unadjusted_HR)
  
  f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx"))
  
  adjusted <- coxph(f_adjusted, cohort) %>%
    tidy(conf.int=TRUE, exponentiate=TRUE) %>%
    filter(term=="studydrugSGLT2") %>%
    mutate(adjusted_HR=paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>%
    select(adjusted_HR)
  
  outcome_hr <- cbind(outcome=i, count, followup, events, unadjusted, adjusted)
  
  all_hrs <- rbind(all_hrs, outcome_hr)
  
}

flextable(all_hrs)
