
# HRs for SGLT2 vs DPP4SU for MACE, HF, CKD, outcomes, hospitalisation, mortality

# Unadjusted and adjusted for age + sex + duration + IMD + QRisk(5yr) + drugline + ncurrtx (adding ethnicity makes ~no difference)


# 5x sensitivity analysis for MACE and HF, only final 4 for expanded MACE, 2 x CKD outcomes and hospitalisation, and only middle 3 for death
## Using narrower definition for MACE/HF as per papers
## Using as-treated (per protocol - use drug stop date + 6 months) rather than ITT
## Only in 1 arm or the other
## Requiring 1 year of previous data and excluding from DPP4SU arm if SGLT2 use in this year
## Death as competing risk (slow to run)

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

cohort <- add_surv_vars(cohort)


## C Just keep variables of interest

cohort <- cohort %>%
  
  select(patid, malesex, ethnicity_5cat_decoded, imd2015_10, regstartdate, gp_record_end, death_date, drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, prehba1c, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preckdstage, contains("cens"), qrisk2_5yr_score, last_sglt2_stop)

rm(list=setdiff(ls(), "cohort"))


############################################################################################

# 2 Define cohorts for sensitivity analyses

## Only in 1 arm or the other
## Requiring 1 year of previous data and excluding from DPP4SU arm if SGLT2 use in this year


## Only in 1 arm or the other

unique_patid_cohort <- cohort %>%
  group_by(patid) %>%
  mutate(earliest_drug_start=min(dstartdate, na.rm=TRUE)) %>%
  ungroup() %>%
  filter(dstartdate==earliest_drug_start)
#118,585

unique_patid_cohort %>% distinct(patid) %>% count()
#118,585

table(unique_patid_cohort$studydrug)
#SGLT2: 29,671
#DPP4SU: 88,914



## Need 1 year of previous data, and exclude from DPP4SU arm if took SGLT2 in this year

one_year_reg_cohort <- cohort %>%
  mutate(reg_before_drug_start=as.numeric(difftime(dstartdate, regstartdate, units="days")),
         time_since_sglt2_stop=as.numeric(difftime(dstartdate, last_sglt2_stop, units="days"))) %>%
  filter(reg_before_drug_start>=365 & (is.na(last_sglt2_stop) | last_sglt2_stop>=365))

table(one_year_reg_cohort$studydrug)
#SGLT2: 47,127
#DPP4SU: 87,209


############################################################################################

# 3 Look at hazard ratios for all sensitivity analyses except death as competing risk (see next section)

outcomes <- c("mace", "narrow_mace", "mace_pp", "expanded_mace", "expanded_mace_pp", "hf", "narrow_hf", "hf_pp", "ckd_345", "ckd_345_pp", "ckd_egfr40", "ckd_egfr40_pp", "hosp", "hosp_pp", "death", "death_pp")

cohorts <- c("cohort", "unique_patid_cohort", "one_year_reg_cohort")

all_hrs <- data.frame()

for (i in outcomes) {
  
  censvar_var=paste0(i, "_censvar")
  censtime_var=paste0(i, "_censtime_yrs")
  
  for (j in cohorts) {
    
    count <- get(j) %>%
      group_by(studydrug) %>%
      summarise(count=n()) %>%
      pivot_wider(names_from=studydrug,
                  names_glue="{studydrug}_count",
                  values_from=count)
    
    followup <- get(j) %>%
      group_by(studydrug) %>%
      summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
      pivot_wider(names_from=studydrug,
                  names_glue="{studydrug}_followup",
                  values_from=time)

    events <- get(j) %>%
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

    unadjusted <- coxph(f, get(j)) %>%
      tidy(conf.int=TRUE, exponentiate=TRUE) %>%
      filter(term=="studydrugSGLT2") %>%
      mutate(unadjusted_HR=paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>%
      select(unadjusted_HR)

    f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx"))

    adjusted <- coxph(f_adjusted, get(j)) %>%
      tidy(conf.int=TRUE, exponentiate=TRUE) %>%
      filter(term=="studydrugSGLT2") %>%
      mutate(adjusted_HR=paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>%
      select(adjusted_HR)

    outcome_hr <- cbind(outcome=i, cohort=j, count, followup, events, unadjusted, adjusted)

    all_hrs <- rbind(all_hrs, outcome_hr)
    
  }
  
}

all_hrs <- all_hrs %>%
  filter(!((grepl("narrow", outcome) | grepl("pp", outcome)) & cohort!="cohort"))

flextable(all_hrs)


############################################################################################

# 3 Death as competing risk for MACE, HF and hospitalisation

library(tidycmprsk)
library(rms)
library(ggsurvfit)
library(gtsummary)


# A) MACE

# Add death as a competing risk
## New variable = 0 if censored for other reasons than death, 1 if outcome, 2 if death
cohort <- cohort %>% 
  mutate(mace_censvar.cr=as.factor(ifelse(mace_censvar==0 & !is.na(death_date) & death_date==mace_censdate, 2, mace_censvar)))

describe(cohort$mace_censvar.cr)
## 3.3% getting outcome, 1.4% dying of other causes
table(cohort$studydrug, cohort$mace_censvar.cr)
prop.table(table(cohort$studydrug, cohort$mace_censvar.cr), margin=1)


# # Plot cumulative incidence adjusted for competing risk of death
# cuminc(Surv(mace_censtime_yrs, mace_censvar.cr) ~ studydrug, data = cohort) %>% 
#   ggcuminc() + 
#   labs(
#     x = "Years"
#   ) + 
#   add_confidence_interval() +
#   add_risktable()
# 
# 
# # Just look at over 70s and compare to survival without competing risk
# agecut <- 70
# ggsurvplot(survfit(Surv(mace_censtime_yrs, mace_censvar) ~ studydrug, data = cohort, subset=dstartdate_age>agecut),fun = function(x) {100 - x*100})
# 
# cohort.agecut <- cohort %>% 
#   filter(dstartdate_age>agecut)
# 
# cuminc(Surv(mace_censtime_yrs, mace_censvar.cr) ~ studydrug, data = cohort.agecut) %>% 
#   ggcuminc() + 
#   labs(
#     x = "Years"
#   ) + 
#   add_confidence_interval() +
#   add_risktable()


# 'Hazard ratios'for whole cohort from competing risk regression

# Unadjusted
crr(Surv(mace_censtime_yrs, mace_censvar.cr) ~ studydrug, data = cohort) %>% 
  tbl_regression(exp = TRUE)
##  SGLT2	0.81	0.76, 0.87	<0.001

# Adjusted
crr(Surv(mace_censtime_yrs, mace_censvar.cr) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tbl_regression(exp = TRUE)
##SGLT2	0.86	0.80, 0.93	<0.001




# B) Expanded MACE

# Add death as a competing risk
## New variable = 0 if censored for other reasons than death, 1 if outcome, 2 if death
cohort <- cohort %>% 
  mutate(expanded_mace_censvar.cr=as.factor(ifelse(expanded_mace_censvar==0 & !is.na(death_date) & death_date==expanded_mace_censdate, 2, expanded_mace_censvar)))

describe(cohort$expanded_mace_censvar.cr)
## 4.3% getting outcome, 1.4% death
table(cohort$studydrug, cohort$expanded_mace_censvar.cr)
prop.table(table(cohort$studydrug, cohort$expanded_mace_censvar.cr), margin=1)


# 'Hazard ratios'for whole cohort from competing risk regression

# Unadjusted
crr(Surv(expanded_mace_censtime_yrs, expanded_mace_censvar.cr) ~ studydrug, data = cohort) %>% 
  tbl_regression(exp = TRUE)
##SGLT2	0.83	0.78, 0.88	<0.001

# Adjusted
crr(Surv(expanded_mace_censtime_yrs, expanded_mace_censvar.cr) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tbl_regression(exp = TRUE)
##SGLT2	0.87	0.81, 0.93	<0.001



# C) Heart failure

# Add death as a competing risk
## New variable = 0 if censored for other reasons than death, 1 if outcome, 2 if death
cohort <- cohort %>% 
  mutate(hf_censvar.cr=as.factor(ifelse(hf_censvar==0 & !is.na(death_date) & death_date==hf_censdate, 2, hf_censvar)))

describe(cohort$hf_censvar.cr)
## 1.8% getting outcome, 2.0% death
table(cohort$studydrug, cohort$hf_censvar.cr)
prop.table(table(cohort$studydrug, cohort$hf_censvar.cr), margin=1)


# 'Hazard ratios'for whole cohort from competing risk regression

# Unadjusted
crr(Surv(hf_censtime_yrs, hf_censvar.cr) ~ studydrug, data = cohort) %>% 
  tbl_regression(exp = TRUE)
##SGLT2	0.71	0.65, 0.78	<0.001

# Adjusted
crr(Surv(hf_censtime_yrs, hf_censvar.cr) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tbl_regression(exp = TRUE)
##SGLT2	0.77	0.69, 0.86	<0.001



# D) CKD stage 3a-5

# Add death as a competing risk
## New variable = 0 if censored for other reasons than death, 1 if outcome, 2 if death
cohort <- cohort %>% 
  mutate(ckd_345_censvar.cr=as.factor(ifelse(ckd_345_censvar==0 & !is.na(death_date) & death_date==ckd_345_censdate, 2, ckd_345_censvar)))

describe(cohort$ckd_345_censvar.cr)
## 2.8% getting outcome, 2.3% death
table(cohort$studydrug, cohort$ckd_345_censvar.cr)
prop.table(table(cohort$studydrug, cohort$ckd_345_censvar.cr), margin=1)


# 'Hazard ratios'for whole cohort from competing risk regression

# Unadjusted
crr(Surv(ckd_345_censtime_yrs, ckd_345_censvar.cr) ~ studydrug, data = cohort) %>% 
  tbl_regression(exp = TRUE)
##SGLT2	0.42	0.39, 0.46	<0.001

# Adjusted
crr(Surv(ckd_345_censtime_yrs, ckd_345_censvar.cr) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tbl_regression(exp = TRUE)
##SGLT2	0.46	0.42, 0.51	<0.001



# C) CKD eGFR decline <=40%/CKD stage 5

# Add death as a competing risk
## New variable = 0 if censored for other reasons than death, 1 if outcome, 2 if death
cohort <- cohort %>% 
  mutate(ckd_egfr40_censvar.cr=as.factor(ifelse(ckd_egfr40_censvar==0 & !is.na(death_date) & death_date==ckd_egfr40_censdate, 2, ckd_egfr40_censvar)))

describe(cohort$ckd_egfr40_censvar.cr)
## 0.6% getting outcome, 2.3% death
table(cohort$studydrug, cohort$ckd_egfr40_censvar.cr)
prop.table(table(cohort$studydrug, cohort$ckd_egfr40_censvar.cr), margin=1)


# 'Hazard ratios'for whole cohort from competing risk regression

# Unadjusted
crr(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar.cr) ~ studydrug, data = cohort) %>% 
  tbl_regression(exp = TRUE)
##SGLT2	0.53	0.45, 0.64	<0.001

# Adjusted
crr(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar.cr) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tbl_regression(exp = TRUE)
##SGLT2	0.51	0.42, 0.64	<0.001



# F) Hospitalisation

# Add death as a competing risk
## New variable = 0 if censored for other reasons than death, 1 if outcome, 2 if death
cohort <- cohort %>% 
  mutate(hosp_censvar.cr=as.factor(ifelse(hosp_censvar==0 & !is.na(death_date) & death_date==hosp_censdate, 2, hosp_censvar)))

describe(cohort$hosp_censvar.cr)
## 18.7% getting outcome, 0.5% death
table(cohort$studydrug, cohort$hosp_censvar.cr)
prop.table(table(cohort$studydrug, cohort$hosp_censvar.cr), margin=1)


# 'Hazard ratios'for whole cohort from competing risk regression

# Unadjusted
crr(Surv(hosp_censtime_yrs, hosp_censvar.cr) ~ studydrug, data = cohort) %>% 
  tbl_regression(exp = TRUE)
## SGLT2	0.86	0.84, 0.88	<0.001

# Adjusted
crr(Surv(hosp_censtime_yrs, hosp_censvar.cr) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tbl_regression(exp = TRUE)
##SGLT2	0.87	0.84, 0.90	<0.001
