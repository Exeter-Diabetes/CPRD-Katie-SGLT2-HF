
# HRs for SGLT2 vs DPP4SU for MACE, HF, all-cause hospitalisation and all-cause mortality

# Unadjusted and adjusted for age + sex + duration + IMD + QRisk(5yr) + drugline + ncurrtx (adding ethnicity makes ~no difference)

# 5x sensitivity analysis for MACE and HF, only final 4 for hospitalisation and middle 3 for death
## Using narrower definition for MACE/HF as per papers
## Using as-treated (per protocol - use drug stop date + 6 months) rather than ITT
## Only in 1 arm or the other
## Requiring 1 year of previous data and excluding from DPP4SU arm if SGLT2 use in this year
## Death as competing risk (slow to run)

# Uses Rda dataset from Slade + cohort_definition_function as starting point

# Not including GLP1 for now


############################################################################################

# Setup
library(tidyverse)
library(gtsummary)
library(survival)
library(survminer)
library(broom)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# 1 Cohort selection and variable setup

## See cohort_definition_function scripts for details

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Raw data/")
load("20221212_t2d_1stinstance.Rda")
load("20221212_t2d_all_drug_periods.Rda")

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Scripts/")
source("cohort_definition_function.R")

cohort <- define_cohort(t2d_1stinstance, t2d_all_drug_periods)

write.table(count(cohort, studydrug), quote=F, col.names=F, row.names=F, sep=" ")
# DPP4SU 90904
# GLP1 12607
# SGLT2 48304

cohort <- cohort %>%
  filter(studydrug!="GLP1") %>%
  select(patid, malesex, ethnicity_5cat_decoded, imd2015_10, regstartdate, gp_record_end, death_date, drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, prehba1c, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preckdstage, qrisk2_smoking_cat, postdrug_first_myocardialinfarction, postdrug_first_primary_incident_mi, postdrug_first_stroke, postdrug_first_primary_incident_stroke, postdrug_first_heartfailure, postdrug_first_primary_hhf, postdrug_first_all_cause_hosp, next_glp1_start, next_sglt2_start, next_tzd_start, last_sglt2_stop, cv_death_date_primary_cause, cv_death_date_any_cause, hf_death_date_primary_cause, hf_death_date_any_cause, qrisk2_lin_predictor, qrisk2_5yr_score, qrisk2_10yr_score, qdiabeteshf_lin_predictor, qdiabeteshf_5yr_score)

rm(list=setdiff(ls(), "cohort"))


############################################################################################

# 2 Make outcome and survival analysis variables for main analysis
## A) Broad MACE (hospitalisation/death any cause or in GP records) intention to treat (ITT)
## B) Specific MACE (incident MI and stroke codes and CV death in HES or ONS death primary cause position, no GP records) ITT
## C) Broad MACE per-protocol (PP)
## D) Broad heart failure (hospitalisation/death any cause or in GP records) ITT
## E) Specific heart failure (hospitalisation/death primary cause; no GP records) ITT
## F) Broad heart failure PP
## G) All-cause hospitalisation ITT
## H) All-cause hospitalisation PP
## I) All-cause mortality ITT
## J) All-cause mortality PP

## For per-protocol: use drug stop date + 6 months


cohort <- cohort %>%
  
  mutate(
    
    # broad_mace = GP (MI/stroke) + expanded HES codelist (MI/stroke; any cause) + CV death (any cause)
    postdrug_broad_mace=pmin(postdrug_first_myocardialinfarction, postdrug_first_stroke, cv_death_date_any_cause, na.rm=TRUE),
    
    # specific_mace = condensed HES codelist (MI/stroke; primary cause) + CV death (primary cause)
    postdrug_specific_mace=pmin(postdrug_first_primary_incident_mi, postdrug_first_primary_incident_stroke, cv_death_date_primary_cause, na.rm=TRUE),
    
    # broad_hf = GP + HES (any cause) + death (any cause)
    postdrug_broad_hf=pmin(postdrug_first_heartfailure, hf_death_date_any_cause, na.rm=TRUE),
    
    # specific_hf = HES (primary cause) + death (primary cause)
    postdrug_specific_hf=pmin(postdrug_first_primary_hhf, hf_death_date_primary_cause, na.rm=TRUE)
    
  )


# A)  Broad MACE ITT

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



# B)  Specific MACE ITT

# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Start SGLT2 (not SGLT2 arm), GLP1, or TZD
## End of GP records
## Specific MACE outcome

cohort <- cohort %>%
  mutate(mace_specific_censdate=if_else(studydrug=="SGLT2",
                                     pmin(five_years_post_dstart,
                                          death_date,
                                          next_glp1_start,
                                          next_tzd_start,
                                          gp_record_end,
                                          postdrug_specific_mace, na.rm=TRUE),
                                     if_else(studydrug=="DPP4SU",
                                             pmin(five_years_post_dstart,
                                                  death_date,
                                                  next_sglt2_start,
                                                  next_glp1_start,
                                                  next_tzd_start,
                                                  gp_record_end,
                                                  postdrug_specific_mace, na.rm=TRUE),
                                             as.Date(NA))),
         
         mace_specific_censvar=ifelse(!is.na(postdrug_specific_mace) & mace_specific_censdate==postdrug_specific_mace, 1, 0),
         
         mace_specific_censtime_yrs=as.numeric(difftime(mace_specific_censdate, dstartdate, unit="days"))/365.25)



# C)  Broad MACE PP

# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Start SGLT2 (not SGLT2 arm), GLP1, or TZD
## End of GP records
## Stop taking study drug (requires gap of at least 6 months) + 6 months
## Broad MACE outcome

cohort <- cohort %>%
  mutate(mace_broad_pp_censdate=if_else(studydrug=="SGLT2",
                                     pmin(five_years_post_dstart,
                                          death_date,
                                          next_glp1_start,
                                          next_tzd_start,
                                          gp_record_end,
                                          postdrug_broad_mace,
                                          (dstopdate+183), na.rm=TRUE),
                                     if_else(studydrug=="DPP4SU",
                                             pmin(five_years_post_dstart,
                                                  death_date,
                                                  next_sglt2_start,
                                                  next_glp1_start,
                                                  next_tzd_start,
                                                  gp_record_end,
                                                  postdrug_broad_mace,
                                                  (dstopdate+183), na.rm=TRUE),
                                             as.Date(NA))),
         
         mace_broad_pp_censvar=ifelse(!is.na(postdrug_broad_mace) & mace_broad_pp_censdate==postdrug_broad_mace, 1, 0),
         
         mace_broad_pp_censtime_yrs=as.numeric(difftime(mace_broad_pp_censdate, dstartdate, unit="days"))/365.25)



# D) Broad heart failure ITT

# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Start SGLT2 (not SGLT2 arm), GLP1, or TZD
## End of GP records
## Broad heart failure outcome

cohort <- cohort %>%
  mutate(hf_broad_censdate=if_else(studydrug=="SGLT2",
                                   pmin(five_years_post_dstart,
                                        death_date,
                                        next_glp1_start,
                                        next_tzd_start,
                                        gp_record_end,
                                        postdrug_broad_hf, na.rm=TRUE),
                                   if_else(studydrug=="DPP4SU",
                                           pmin(five_years_post_dstart,
                                                death_date,
                                                next_sglt2_start,
                                                next_glp1_start,
                                                next_tzd_start,
                                                gp_record_end,
                                                postdrug_broad_hf, na.rm=TRUE),
                                           as.Date(NA))),
         
         hf_broad_censvar=ifelse(!is.na(postdrug_broad_hf) & hf_broad_censdate==postdrug_broad_hf, 1, 0),
         
         hf_broad_censtime_yrs=as.numeric(difftime(hf_broad_censdate, dstartdate, unit="days"))/365.25)



# E) Specific heart failure ITT

# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Start SGLT2 (not SGLT2 arm), GLP1, or TZD
## End of GP records
## Specific heart failure outcome

cohort <- cohort %>%
  mutate(hf_specific_censdate=if_else(studydrug=="SGLT2",
                                   pmin(five_years_post_dstart,
                                        death_date,
                                        next_glp1_start,
                                        next_tzd_start,
                                        gp_record_end,
                                        postdrug_specific_hf, na.rm=TRUE),
                                   if_else(studydrug=="DPP4SU",
                                           pmin(five_years_post_dstart,
                                                death_date,
                                                next_sglt2_start,
                                                next_glp1_start,
                                                next_tzd_start,
                                                gp_record_end,
                                                postdrug_specific_hf, na.rm=TRUE),
                                           as.Date(NA))),
         
         hf_specific_censvar=ifelse(!is.na(postdrug_specific_hf) & hf_specific_censdate==postdrug_specific_hf, 1, 0),
         
         hf_specific_censtime_yrs=as.numeric(difftime(hf_specific_censdate, dstartdate, unit="days"))/365.25)



# F) Broad heart failure PP

# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Start SGLT2 (not SGLT2 arm), GLP1, or TZD
## End of GP records
## Stop taking study drug (requires gap of at least 6 months) + 6 months
## Broad heart failure outcome

cohort <- cohort %>%
  mutate(hf_broad_pp_censdate=if_else(studydrug=="SGLT2",
                                   pmin(five_years_post_dstart,
                                        death_date,
                                        next_glp1_start,
                                        next_tzd_start,
                                        gp_record_end,
                                        postdrug_broad_hf, 
                                        (dstopdate+183), na.rm=TRUE),
                                   if_else(studydrug=="DPP4SU",
                                           pmin(five_years_post_dstart,
                                                death_date,
                                                next_sglt2_start,
                                                next_glp1_start,
                                                next_tzd_start,
                                                gp_record_end,
                                                postdrug_broad_hf,
                                                (dstopdate+183), na.rm=TRUE),
                                           as.Date(NA))),
         
         hf_broad_pp_censvar=ifelse(!is.na(postdrug_broad_hf) & hf_broad_pp_censdate==postdrug_broad_hf, 1, 0),
         
         hf_broad_pp_censtime_yrs=as.numeric(difftime(hf_broad_pp_censdate, dstartdate, unit="days"))/365.25)



# G) All-cause hospitalisation ITT

# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Start SGLT2 (not SGLT2 arm), GLP1, or TZD
## End of GP records
## Any cause hospitalisation

cohort <- cohort %>%
  mutate(hosp_censdate=if_else(studydrug=="SGLT2",
                               pmin(five_years_post_dstart,
                                    death_date,
                                    next_glp1_start,
                                    next_tzd_start,
                                    gp_record_end,
                                    postdrug_first_all_cause_hosp, na.rm=TRUE),
                               if_else(studydrug=="DPP4SU",
                                       pmin(five_years_post_dstart,
                                            death_date,
                                            next_sglt2_start,
                                            next_glp1_start,
                                            next_tzd_start,
                                            gp_record_end,
                                            postdrug_first_all_cause_hosp, na.rm=TRUE),
                                       as.Date(NA))),
         
         hosp_censvar=ifelse(!is.na(postdrug_first_all_cause_hosp) & hosp_censdate==postdrug_first_all_cause_hosp, 1, 0),
         
         hosp_censtime_yrs=as.numeric(difftime(hosp_censdate, dstartdate, unit="days"))/365.25)



# H) All-cause hospitalisation PP

# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Start SGLT2 (not SGLT2 arm), GLP1, or TZD
## End of GP records
## Stop taking study drug (requires gap of at least 6 months) + 6 months
## Any cause hospitalisation

cohort <- cohort %>%
  mutate(hosp_pp_censdate=if_else(studydrug=="SGLT2",
                               pmin(five_years_post_dstart,
                                    death_date,
                                    next_glp1_start,
                                    next_tzd_start,
                                    gp_record_end,
                                    postdrug_first_all_cause_hosp,
                                    (dstopdate+183), na.rm=TRUE),
                               if_else(studydrug=="DPP4SU",
                                       pmin(five_years_post_dstart,
                                            death_date,
                                            next_sglt2_start,
                                            next_glp1_start,
                                            next_tzd_start,
                                            gp_record_end,
                                            postdrug_first_all_cause_hosp,
                                            (dstopdate+183), na.rm=TRUE),
                                       as.Date(NA))),
         
         hosp_pp_censvar=ifelse(!is.na(postdrug_first_all_cause_hosp) & hosp_pp_censdate==postdrug_first_all_cause_hosp, 1, 0),
         
         hosp_pp_censtime_yrs=as.numeric(difftime(hosp_pp_censdate, dstartdate, unit="days"))/365.25)



# I) All-cause mortality ITT

# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Start SGLT2 (not SGLT2 arm), GLP1, or TZD
## End of GP records

cohort <- cohort %>%
  mutate(death_censdate=if_else(studydrug=="SGLT2",
                                pmin(five_years_post_dstart,
                                     death_date,
                                     next_glp1_start,
                                     next_tzd_start,
                                     gp_record_end, na.rm=TRUE),
                                if_else(studydrug=="DPP4SU",
                                        pmin(five_years_post_dstart,
                                             death_date,
                                             next_sglt2_start,
                                             next_glp1_start,
                                             next_tzd_start,
                                             gp_record_end, na.rm=TRUE),
                                        as.Date(NA))),
         
         death_censvar=ifelse(!is.na(death_date) & death_censdate==death_date, 1, 0),
         
         death_censtime_yrs=as.numeric(difftime(death_censdate, dstartdate, unit="days"))/365.25)



# J) All-cause mortality PP

# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Start SGLT2 (not SGLT2 arm), GLP1, or TZD
## Stop taking study drug (requires gap of at least 6 months) + 6 months
## End of GP records

cohort <- cohort %>%
  mutate(death_pp_censdate=if_else(studydrug=="SGLT2",
                                pmin(five_years_post_dstart,
                                     death_date,
                                     next_glp1_start,
                                     next_tzd_start,
                                     gp_record_end,
                                     (dstopdate+183), na.rm=TRUE),
                                if_else(studydrug=="DPP4SU",
                                        pmin(five_years_post_dstart,
                                             death_date,
                                             next_sglt2_start,
                                             next_glp1_start,
                                             next_tzd_start,
                                             gp_record_end,
                                             (dstopdate+183), na.rm=TRUE),
                                        as.Date(NA))),
         
         death_pp_censvar=ifelse(!is.na(death_date) & death_pp_censdate==death_date, 1, 0),
         
         death_pp_censtime_yrs=as.numeric(difftime(death_pp_censdate, dstartdate, unit="days"))/365.25)



############################################################################################

# 3 Define cohorts for sensitivity analyses

## Only in 1 arm or the other
## Requiring 1 year of previous data and excluding from DPP4SU arm if SGLT2 use in this year


## Only in 1 arm or the other

unique_patid_cohort <- cohort %>%
  group_by(patid) %>%
  mutate(earliest_drug_start=min(dstartdate, na.rm=TRUE)) %>%
  ungroup() %>%
  filter(dstartdate==earliest_drug_start)
#118,653

unique_patid_cohort %>% distinct(patid) %>% count()
#118,653

table(unique_patid_cohort$studydrug)
#SGLT2: 29,684
#DPP4SU: 88,969



## Need 1 year of previous data, and exclude from DPP4SU arm if took SGLT2 in this year

one_year_reg_cohort <- cohort %>%
  mutate(reg_before_drug_start=as.numeric(difftime(dstartdate, regstartdate, units="days")),
         time_since_sglt2_stop=as.numeric(difftime(dstartdate, last_sglt2_stop, units="days"))) %>%
  filter(reg_before_drug_start>=365 & (is.na(last_sglt2_stop) | last_sglt2_stop>=365))

table(one_year_reg_cohort$studydrug)
#SGLT2: 47,151
#DPP4SU: 87,261


############################################################################################

# 4 Look at hazard ratios for all sensitivity analyses except death as competing risk (see next section)

# MACE

## 'Specific' outcome definition
cohort %>% group_by(studydrug) %>% summarise(time=median(mace_specific_censtime_yrs))

cohort %>% group_by(studydrug) %>% summarise(events=sum(mace_specific_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(mace_specific_censtime_yrs, mace_specific_censvar) ~ studydrug, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(mace_specific_censtime_yrs, mace_specific_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))


## Per-protocol
cohort %>% group_by(studydrug) %>% summarise(time=median(mace_broad_pp_censtime_yrs))

cohort %>% group_by(studydrug) %>% summarise(events=sum(mace_broad_pp_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(mace_broad_pp_censtime_yrs, mace_broad_pp_censvar) ~ studydrug, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(mace_broad_pp_censtime_yrs, mace_broad_pp_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))


## Patients can only be in 1 arm
unique_patid_cohort %>% group_by(studydrug) %>% summarise(time=median(mace_broad_censtime_yrs))

unique_patid_cohort %>% group_by(studydrug) %>% summarise(events=sum(mace_broad_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ studydrug, data = unique_patid_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = unique_patid_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))


## Require 1 year previous data and exclude from DPP4SU if SGLT2 in year prior
one_year_reg_cohort %>% group_by(studydrug) %>% summarise(time=median(mace_broad_censtime_yrs))

one_year_reg_cohort %>% group_by(studydrug) %>% summarise(events=sum(mace_broad_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ studydrug, data = one_year_reg_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = one_year_reg_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))



# HF

## 'Specific' outcome definition
cohort %>% group_by(studydrug) %>% summarise(time=median(hf_specific_censtime_yrs))

cohort %>% group_by(studydrug) %>% summarise(events=sum(hf_specific_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(hf_specific_censtime_yrs, hf_specific_censvar) ~ studydrug, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(hf_specific_censtime_yrs, hf_specific_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))


## Per-protocol
cohort %>% group_by(studydrug) %>% summarise(time=median(hf_broad_pp_censtime_yrs))

cohort %>% group_by(studydrug) %>% summarise(events=sum(hf_broad_pp_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(hf_broad_pp_censtime_yrs, hf_broad_pp_censvar) ~ studydrug, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(hf_broad_pp_censtime_yrs, hf_broad_pp_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))


## Patients can only be in 1 arm
unique_patid_cohort %>% group_by(studydrug) %>% summarise(time=median(hf_broad_censtime_yrs))

unique_patid_cohort %>% group_by(studydrug) %>% summarise(events=sum(hf_broad_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(hf_broad_censtime_yrs, hf_broad_censvar) ~ studydrug, data = unique_patid_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(hf_broad_censtime_yrs, hf_broad_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = unique_patid_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))


## Require 1 year previous data and exclude from DPP4SU if SGLT2 in year prior
one_year_reg_cohort %>% group_by(studydrug) %>% summarise(time=median(hf_broad_censtime_yrs))

one_year_reg_cohort %>% group_by(studydrug) %>% summarise(events=sum(hf_broad_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(hf_broad_censtime_yrs, hf_broad_censvar) ~ studydrug, data = one_year_reg_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(hf_broad_censtime_yrs, hf_broad_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = one_year_reg_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))



# Hospitalisation

## Per-protocol
cohort %>% group_by(studydrug) %>% summarise(time=median(hosp_pp_censtime_yrs))

cohort %>% group_by(studydrug) %>% summarise(events=sum(hosp_pp_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(hosp_pp_censtime_yrs, hosp_pp_censvar) ~ studydrug, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(hosp_pp_censtime_yrs, hosp_pp_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))


## Patients can only be in 1 arm
unique_patid_cohort %>% group_by(studydrug) %>% summarise(time=median(hosp_censtime_yrs))

unique_patid_cohort %>% group_by(studydrug) %>% summarise(events=sum(hosp_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(hosp_censtime_yrs, hosp_censvar) ~ studydrug, data = unique_patid_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(hosp_censtime_yrs, hosp_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = unique_patid_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))


## Require 1 year previous data and exclude from DPP4SU if SGLT2 in year prior
one_year_reg_cohort %>% group_by(studydrug) %>% summarise(time=median(hosp_censtime_yrs))

one_year_reg_cohort %>% group_by(studydrug) %>% summarise(events=sum(hosp_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(hosp_censtime_yrs, hosp_censvar) ~ studydrug, data = one_year_reg_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(hosp_censtime_yrs, hosp_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = one_year_reg_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))



# Death

## Per-protocol
cohort %>% group_by(studydrug) %>% summarise(time=median(death_pp_censtime_yrs))

cohort %>% group_by(studydrug) %>% summarise(events=sum(death_pp_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(death_pp_censtime_yrs, death_pp_censvar) ~ studydrug, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(death_pp_censtime_yrs, death_pp_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))


## Patients can only be in 1 arm
unique_patid_cohort %>% group_by(studydrug) %>% summarise(time=median(death_censtime_yrs))

unique_patid_cohort %>% group_by(studydrug) %>% summarise(events=sum(death_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(death_censtime_yrs, death_censvar) ~ studydrug, data = unique_patid_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(death_censtime_yrs, death_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = unique_patid_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))


## Require 1 year previous data and exclude from DPP4SU if SGLT2 in year prior
one_year_reg_cohort %>% group_by(studydrug) %>% summarise(time=median(death_censtime_yrs))

one_year_reg_cohort %>% group_by(studydrug) %>% summarise(events=sum(death_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(death_censtime_yrs, death_censvar) ~ studydrug, data = one_year_reg_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(death_censtime_yrs, death_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = one_year_reg_cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))


############################################################################################

# 5 Death as competing risk for MACE, HF and hospitalisation

library(tidycmprsk)
library(rms)
library(ggsurvfit)


# A) MACE

# Add death as a competing risk
## New variable = 0 if censored for other reasons than death, 1 if outcome, 2 if death
cohort <- cohort %>% 
  mutate(mace_broad_censvar.cr=ifelse(mace_broad_censvar==0 & death_date==mace_broad_censdate, 2, mace_broad_censvar))

cohort <- cohort %>%
  mutate(mace_broad_censvar.cr=ifelse(is.na(mace_broad_censvar.cr), 0, mace_broad_censvar.cr))

table(cohort$mace_broad_censvar.cr)
describe(cohort$mace_broad_censvar.cr)
## 3.3% getting outcome, 1.4% dying of other causes


cohort <- cohort %>%
  mutate(mace_broad_censvar.cr=as.factor(mace_broad_censvar.cr))


# Plot cumulative incidence adjusted for competing risk of death
cuminc(Surv(mace_broad_censtime_yrs, mace_broad_censvar.cr) ~ studydrug, data = cohort) %>% 
  ggcuminc() + 
  labs(
    x = "Years"
  ) + 
  add_confidence_interval() +
  add_risktable()


# Just look at over 70s and compare to survival without competing risk
agecut <- 70
ggsurvplot(survfit(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ studydrug, data = cohort, subset=dstartdate_age>agecut),fun = function(x) {100 - x*100})

cohort.agecut <- cohort %>% 
  filter(dstartdate_age>agecut)

cuminc(Surv(mace_broad_censtime_yrs, mace_broad_censvar.cr) ~ studydrug, data = cohort.agecut) %>% 
  ggcuminc() + 
  labs(
    x = "Years"
  ) + 
  add_confidence_interval() +
  add_risktable()



# Adjusted 'hazard ratios'

# competing risk regression
crr(Surv(hf_broad_itt_censtime_yrs, hf_broad_itt_censvar.cr) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tbl_regression(exp = TRUE)

studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx

##SGLT2	0.76	0.68, 0.85	<0.001

# coxph comparison
coxph(Surv(hf_broad_itt_censtime_yrs, hf_broad_itt_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## studydrugSGLT2           0.759    0.680     0.848




# B) Heart failure

# Add death as a competing risk
## New variable = 0 if censored for other reasons than death, 1 if outcome, 2 if death
cohort <- cohort %>% 
  mutate(hf_broad_censvar.cr=ifelse(hf_broad_censvar==0 & death_date==hf_broad_censdate, 2, hf_broad_censvar))

cohort <- cohort %>%
  mutate(hf_broad_censvar.cr=ifelse(is.na(hf_broad_censvar.cr), 0, hf_broad_censvar.cr))

table(cohort$hf_broad_censvar.cr)
describe(cohort$hf_broad_censvar.cr)
## 3.3% getting outcome, 1.4% dying of other causes


cohort <- cohort %>%
  mutate(hf_broad_censvar.cr=as.factor(hf_broad_censvar.cr))


# Plot cumulative incidence adjusted for competing risk of death
cuminc(Surv(hf_broad_censtime_yrs, hf_broad_censvar.cr) ~ studydrug, data = cohort) %>% 
  ggcuminc() + 
  labs(
    x = "Years"
  ) + 
  add_confidence_interval() +
  add_risktable()


# Just look at over 70s and compare to survival without competing risk
agecut <- 70
ggsurvplot(survfit(Surv(hf_broad_censtime_yrs, hf_broad_censvar) ~ studydrug, data = cohort, subset=dstartdate_age>agecut),fun = function(x) {100 - x*100})

cohort.agecut <- cohort %>% 
  filter(dstartdate_age>agecut)

cuminc(Surv(hf_broad_censtime_yrs, hf_broad_censvar.cr) ~ studydrug, data = cohort.agecut) %>% 
  ggcuminc() + 
  labs(
    x = "Years"
  ) + 
  add_confidence_interval() +
  add_risktable()



# Adjusted 'hazard ratios'

# competing risk regression
crr(Surv(hf_broad_itt_censtime_yrs, hf_broad_itt_censvar.cr) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tbl_regression(exp = TRUE)

##SGLT2	0.76	0.68, 0.85	<0.001

# coxph comparison
coxph(Surv(hf_broad_itt_censtime_yrs, hf_broad_itt_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## studydrugSGLT2           0.759    0.680     0.848




# C) Hospitalisation

# Add death as a competing risk
## New variable = 0 if censored for other reasons than death, 1 if outcome, 2 if death
cohort <- cohort %>% 
  mutate(hosp_censvar.cr=ifelse(hosp_censvar==0 & death_date==hosp_censdate, 2, hosp_censvar))

cohort <- cohort %>%
  mutate(hosp_censvar.cr=ifelse(is.na(hosp_censvar.cr), 0, hosp_censvar.cr))

table(cohort$hosp_censvar.cr)
describe(cohort$hosp_censvar.cr)
## 3.3% getting outcome, 1.4% dying of other causes


cohort <- cohort %>%
  mutate(hosp_censvar.cr=as.factor(hosp_censvar.cr))


# Plot cumulative incidence adjusted for competing risk of death
cuminc(Surv(hosp_censtime_yrs, hosp_censvar.cr) ~ studydrug, data = cohort) %>% 
  ggcuminc() + 
  labs(
    x = "Years"
  ) + 
  add_confidence_interval() +
  add_risktable()


# Just look at over 70s and compare to survival without competing risk
agecut <- 70
ggsurvplot(survfit(Surv(hosp_censtime_yrs, hosp_censvar) ~ studydrug, data = cohort, subset=dstartdate_age>agecut),fun = function(x) {100 - x*100})

cohort.agecut <- cohort %>% 
  filter(dstartdate_age>agecut)

cuminc(Surv(hosp_censtime_yrs, hosp_censvar.cr) ~ studydrug, data = cohort.agecut) %>% 
  ggcuminc() + 
  labs(
    x = "Years"
  ) + 
  add_confidence_interval() +
  add_risktable()



# Adjusted 'hazard ratios'

# competing risk regression
crr(Surv(hosp_itt_censtime_yrs, hosp_itt_censvar.cr) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tbl_regression(exp = TRUE)

##SGLT2	0.76	0.68, 0.85	<0.001

# coxph comparison
coxph(Surv(hosp_itt_censtime_yrs, hosp_itt_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## studydrugSGLT2           0.759    0.680     0.848
