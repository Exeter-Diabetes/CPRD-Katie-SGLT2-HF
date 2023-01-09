
# HRs for SGLT2 vs DPP4SU for 3-point MACE+HF, MACE, HF, all-cause hospitalisation and all-cause mortality

# Unadjusted and adjusted for age + sex + duration + IMD + QRisk(5yr) + drugline + ncurrtx (adding ethnicity makes ~no difference)

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
# DPP4SU 90853
# GLP1 12604
# SGLT2 48279

cohort <- cohort %>%
  filter(studydrug!="GLP1") %>%
  select(patid, malesex, ethnicity_5cat_decoded, imd2015_10, regstartdate, gp_record_end, death_date, drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, prehba1c, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preckdstage, qrisk2_smoking_cat, postdrug_first_myocardialinfarction, postdrug_first_primary_incident_mi, postdrug_first_stroke, postdrug_first_primary_incident_stroke, postdrug_first_heartfailure, postdrug_first_primary_hhf, postdrug_first_all_cause_hosp, next_glp1_start, next_sglt2_start, next_tzd_start, last_sglt2_stop, cv_death_date_primary_cause, cv_death_date_any_cause, hf_death_date_primary_cause, hf_death_date_any_cause, qrisk2_lin_predictor, qrisk2_5yr_score, qrisk2_10yr_score, qdiabeteshf_lin_predictor, qdiabeteshf_5yr_score, contains("statins"))

rm(list=setdiff(ls(), "cohort"))


############################################################################################

# 2 Make outcome and survival analysis variables for main analysis

## A) Combined broad MACE (hospitalisation/death primary or secondary cause or in GP records) and broad heart failure (hospitalisation/death primary or secondary cause or in GP records)
## B) Broad MACE (hospitalisation/death any cause or in GP records)
## C) Broad heart failure (hospitalisation/death any cause or in GP records)
## D) All-cause hospitalisation
## E) All-cause mortality


cohort <- cohort %>%
  
  mutate(
    
    # broad_mace = GP (MI/stroke) + expanded HES codelist (MI/stroke; any cause) + CV death (any cause)
    postdrug_broad_mace=pmin(postdrug_first_myocardialinfarction, postdrug_first_stroke, cv_death_date_any_cause, na.rm=TRUE),
    
    # broad_hf = GP + HES (any cause) + death (any cause)
    postdrug_broad_hf=pmin(postdrug_first_heartfailure, hf_death_date_any_cause, na.rm=TRUE)
    
  )


# A)  Combined broad MACE (hospitalisation/death primary or secondary cause or in GP records) and broad heart failure (hospitalisation/death primary or secondary cause or in GP records)

# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Start SGLT2 (not SGLT2 arm), GLP1, or TZD
## End of GP records
## Broad MACE outcome
## Broad HF outcomes

cohort <- cohort %>%
  mutate(five_years_post_dstart=dstartdate+(365.25*5),
         
         combined_mace_hf_broad_censdate=if_else(studydrug=="SGLT2",
                                                 pmin(five_years_post_dstart,
                                                      death_date,
                                                      next_glp1_start,
                                                      next_tzd_start,
                                                      gp_record_end,
                                                      postdrug_broad_mace,
                                                      postdrug_broad_hf, na.rm=TRUE),
                                                 if_else(studydrug=="DPP4SU",
                                                         pmin(five_years_post_dstart,
                                                              death_date,
                                                              next_sglt2_start,
                                                              next_glp1_start,
                                                              next_tzd_start,
                                                              gp_record_end,
                                                              postdrug_broad_mace,
                                                              postdrug_broad_hf, na.rm=TRUE),
                                                         as.Date(NA))),
         
         combined_mace_hf_broad_censvar=ifelse((!is.na(postdrug_broad_mace) & combined_mace_hf_broad_censdate==postdrug_broad_mace) | (!is.na(postdrug_broad_hf) & combined_mace_hf_broad_censdate==postdrug_broad_hf), 1, 0),
         
         combined_mace_hf_broad_censtime_yrs=as.numeric(difftime(combined_mace_hf_broad_censdate, dstartdate, unit="days"))/365.25)



# B)  Broad MACE

# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Start SGLT2 (not SGLT2 arm), GLP1, or TZD
## End of GP records
## Broad MACE outcome

cohort <- cohort %>%
  mutate(mace_broad_censdate=if_else(studydrug=="SGLT2",
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



# C) Broad heart failure

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



# D) All-cause hospitalisation

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



# E) All-cause mortality

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



############################################################################################

# 3 Look at hazard ratios


# A) MACE+HF

cohort %>% group_by(studydrug) %>% summarise(time=median(combined_mace_hf_broad_censtime_yrs))

cohort %>% group_by(studydrug) %>% summarise(events=sum(combined_mace_hf_broad_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(combined_mace_hf_broad_censtime_yrs, combined_mace_hf_broad_censvar) ~ studydrug, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(combined_mace_hf_broad_censtime_yrs, combined_mace_hf_broad_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted with ethnicity
coxph(Surv(combined_mace_hf_broad_censtime_yrs, combined_mace_hf_broad_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx + ethnicity_5cat_decoded, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))



# B) MACE

cohort %>% group_by(studydrug) %>% summarise(time=median(mace_broad_censtime_yrs))

cohort %>% group_by(studydrug) %>% summarise(events=sum(mace_broad_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ studydrug, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))



# C) Heart failure

cohort %>% group_by(studydrug) %>% summarise(time=median(hf_broad_censtime_yrs))

cohort %>% group_by(studydrug) %>% summarise(events=sum(hf_broad_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(hf_broad_censtime_yrs, hf_broad_censvar) ~ studydrug, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(hf_broad_censtime_yrs, hf_broad_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))



# D) All-cause hospitalisation

cohort %>% group_by(studydrug) %>% summarise(time=median(hosp_censtime_yrs))

cohort %>% group_by(studydrug) %>% summarise(events=sum(hosp_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(hosp_censtime_yrs, hosp_censvar) ~ studydrug, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(hosp_censtime_yrs, hosp_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))



# E) All-cause death

cohort %>% group_by(studydrug) %>% summarise(time=median(death_censtime_yrs))

cohort %>% group_by(studydrug) %>% summarise(events=sum(death_censvar), drug_count=n()) %>% mutate(events_perc=round(events*100/drug_count,1))

## Unadjusted
coxph(Surv(death_censtime_yrs, death_censvar) ~ studydrug, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

## Adjusted
coxph(Surv(death_censtime_yrs, death_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))

coxph(Surv(death_censtime_yrs, death_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qrisk2_5yr_score + drugline_all + ncurrtx + ethnicity_5cat_decoded, data = cohort) %>% 
  tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, starts_with("conf"))
