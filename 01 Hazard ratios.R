
# HRs for SGLT2 vs DPP4SU for:
## 3-point MACE (ischaemic stroke, MI, CV death)
## HF
## all-cause hospitalisation
## all-cause mortality

# Unadjusted and adjusted for lots of things

# Not including GLP1 currently
 
############################################################################################

# Setup
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(flextable)
library(PSweight)
library(gtsummary)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# 1 Cohort selection and variable setup

## A Cohort selection (see cohort_definition function for details)

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Raw data/")
load("20231121_t2d_1stinstance.Rda")
load("20231121_t2d_all_drug_periods.Rda")

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Scripts/Functions")
source("cohort_definition.R")

cohort <- define_cohort(t2d_1stinstance, t2d_all_drug_periods)

cohort <- cohort %>% filter(drugclass!="GLP1") %>% mutate(studydrug=factor(ifelse(drugclass=="SGLT2", "SGLT2", "DPP4SU"), levels=c("DPP4SU", "SGLT2")))
#cohort <- cohort %>% filter(drugclass!="GLP1" & drugclass!="DPP4") %>% mutate(studydrug=factor(drugclass, levels=c("SU", "SGLT2")))
#cohort <- cohort %>% filter(drugclass!="GLP1" & drugclass!="SGLT2") %>% mutate(studydrug=factor(drugclass, levels=c("SU", "DPP4")))

table(cohort$studydrug)
# DPP4SU 96,702
# SGLT2 45,139
 
# SU 36935 
# SGLT2 45139

# SU 36935
# DPP4 59767 

## B Make variables for survival analysis of all endpoints (see survival_variables function for details)

source("survival_variables.R")

cohort <- add_surv_vars(cohort, main_only=FALSE)
         

## C Just keep variables of interest

cohort <- cohort %>%
  
  select(patid, studydrug, drugclass, drugsubstances, dstartdate, dstopdate, dstartdate_age, malesex, ethnicity_qrisk2_decoded, imd2015_10, qrisk2_smoking_cat, drugline_all, ncurrtx, initiation_year, drugorder, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_dm_dur_all, prehba1c, prebmi, pretotalcholesterol, prehdl, preegfr, preckdstage, prealt, presbp, contains("last_6_months"), hosp_admission_prev_year, hosp_admission_prev_year_count, predrug_cld, predrug_hypertension, predrug_af, predrug_copd, predrug_haem_cancer, predrug_solid_cancer, predrug_diabeticnephropathy, predrug_neuropathy, predrug_retinopathy, predrug_dementia, predrug_otherneuroconditions, qdiabeteshf_5yr_score, qdiabeteshf_lin_predictor, contains("cens"), regstartdate, gp_record_end, death_date, primary_death_cause)

rm(list=setdiff(ls(), "cohort"))


# Check follow-up periods for potential overlap

test <- cohort %>%
  select(patid, studydrug, drugclass, dstartdate, dstopdate, hf_censdate) %>%
  group_by(patid) %>%
  arrange(patid, dstartdate) %>%
  mutate(overlap=ifelse(!is.na(lag(hf_censdate)) & dstartdate<lag(hf_censdate), 1, 0))

table(test$overlap)
#0


############################################################################################

# 2 Look at hazard ratios

main_outcomes <- c("mace", "hf", "hosp", "death")

variables <- "dstartdate_age + malesex + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + drugline_all + ncurrtx + initiation_year + dstartdate_dm_dur_all + prehba1c + prebmi + pretotalcholesterol + prehdl + preegfr + prealt + presbp + statin_last_6_months + acei_last_6_months + arb_last_6_months + bb_last_6_months + ccb_last_6_months + thiazide_diuretics_last_6_months + loop_diuretics_last_6_months + ksparing_diuretics_last_6_months + hosp_admission_prev_year_count + predrug_cld + predrug_hypertension + predrug_af + predrug_copd + predrug_haem_cancer + predrug_solid_cancer + predrug_diabeticnephropathy + predrug_neuropathy + predrug_retinopathy + predrug_dementia + predrug_otherneuroconditions + qdiabeteshf_5yr_score"

weight_cohort <- data.frame(cohort) %>% filter(!is.na(dstartdate_age) & !is.na(malesex) & !is.na(ethnicity_qrisk2_decoded) & !is.na(imd2015_10) & !is.na(qrisk2_smoking_cat) & !is.na(drugline_all) & !is.na(ncurrtx) & !is.na(initiation_year) & !is.na(dstartdate_dm_dur_all) & !is.na(prehba1c) & !is.na(prebmi) & !is.na(pretotalcholesterol) & !is.na(prehdl) & !is.na(preegfr) & !is.na(prealt) & !is.na(presbp) & !is.na(statin_last_6_months) & !is.na(acei_last_6_months) & !is.na(arb_last_6_months) & !is.na(bb_last_6_months) & !is.na(ccb_last_6_months) & !is.na(thiazide_diuretics_last_6_months) & !is.na(loop_diuretics_last_6_months) & !is.na(ksparing_diuretics_last_6_months) & !is.na(hosp_admission_prev_year_count) & !is.na(predrug_cld) & !is.na(predrug_hypertension) & !is.na(predrug_af) & !is.na(predrug_copd) & !is.na(predrug_haem_cancer) & !is.na(predrug_solid_cancer) & !is.na(predrug_diabeticnephropathy) & !is.na(predrug_neuropathy) & !is.na(predrug_retinopathy) & !is.na(predrug_dementia) & !is.na(predrug_otherneuroconditions) & !is.na(qdiabeteshf_5yr_score))

all_hrs <- data.frame()

for (i in main_outcomes) {
  
  #activedrug="SGLT2"
  activedrug="DPP4"
  
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

  unadjusted <- coxph(as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug")), cohort) %>%
    tidy(conf.int=TRUE, exponentiate=TRUE) %>%
    filter(term==paste0("studydrug", activedrug)) %>%
    mutate(unadjusted_HR=paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>%
    select(unadjusted_HR)
  
  adjusted <- coxph(as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~ studydrug + ", variables)), cohort) %>%
    tidy(conf.int=TRUE, exponentiate=TRUE) %>%
    filter(term==paste0("studydrug", activedrug)) %>%
    mutate(adjusted_HR=paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>%
    select(adjusted_HR)
  
  
  overlap <- SumStat(ps.formula=formula(paste("studydrug ~ ", variables)), data = weight_cohort, weight = c("IPW", "overlap"))
  
  weight_cohort$overlap_weights <- overlap$ps.weights$overlap
  weight_cohort$iptw_weights <- overlap$ps.weights$IPW
  
  
  unadjusted_overlap_weight <- coxph(as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug")), weight_cohort, weights=overlap_weights) %>%
    tidy(conf.int=TRUE, exponentiate=TRUE) %>%
    filter(term==paste0("studydrug", activedrug)) %>%
    mutate(unadjusted_overlap_weight_HR=paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>%
    select(unadjusted_overlap_weight_HR)
  
  unadjusted_iptw_weight <- coxph(as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug")), weight_cohort, weights=iptw_weights) %>%
    tidy(conf.int=TRUE, exponentiate=TRUE) %>%
    filter(term==paste0("studydrug", activedrug)) %>%
    mutate(unadjusted_iptw_weight_HR=paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>%
    select(unadjusted_iptw_weight_HR)
  
  adjusted_overlap_weight <- coxph(as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~ studydrug + ", variables)), weight_cohort, weights=overlap_weights) %>%
    tidy(conf.int=TRUE, exponentiate=TRUE) %>%
    filter(term==paste0("studydrug", activedrug)) %>%
    mutate(adjusted_overlap_weight_HR=paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>%
    select(adjusted_overlap_weight_HR)
  
  adjusted_iptw_weight <- coxph(as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~ studydrug + ", variables)), weight_cohort, weights=iptw_weights) %>%
    tidy(conf.int=TRUE, exponentiate=TRUE) %>%
    filter(term==paste0("studydrug", activedrug)) %>%
    mutate(adjusted_iptw_weight_HR=paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>%
    select(adjusted_iptw_weight_HR)
  
  outcome_hr <- cbind(outcome=i, count, followup, events, unadjusted, adjusted, unadjusted_overlap_weight, unadjusted_iptw_weight, adjusted_overlap_weight, adjusted_iptw_weight)
  
  all_hrs <- rbind(all_hrs, outcome_hr)
  
}



flextable(all_hrs)



### Baseline characteristics

cohort %>%
  select(drugclass, dstartdate_age, malesex, ethnicity_qrisk2_decoded, imd2015_10, qrisk2_smoking_cat, drugline_all, ncurrtx, initiation_year, dstartdate_dm_dur_all, prehba1c, prebmi, pretotalcholesterol, prehdl, preegfr, prealt, presbp, statin_last_6_months, acei_last_6_months, arb_last_6_months, bb_last_6_months, ccb_last_6_months, thiazide_diuretics_last_6_months, loop_diuretics_last_6_months, ksparing_diuretics_last_6_months, hosp_admission_prev_year_count, predrug_cld, predrug_hypertension, predrug_af, predrug_copd, predrug_haem_cancer, predrug_solid_cancer, predrug_diabeticnephropathy, predrug_neuropathy, predrug_retinopathy, predrug_dementia, predrug_otherneuroconditions, qdiabeteshf_5yr_score) %>%
  tbl_summary(by=drugclass)



  