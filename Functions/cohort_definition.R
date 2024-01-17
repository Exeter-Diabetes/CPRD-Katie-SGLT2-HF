
# Inclusion/exclusion criteria:
## a) T2Ds
## b) With HES linkage
## c) 1st instance
## d) Exclude if start drug within 90 days of registration

## e) Aged 18+
## f) SGLT2/DPP4/SU
## g) Initiated between 01/01/2013 and end of data (31/10/2020)
## h) Exclude if first line
## i) Exclude if also on insulin/GLP1/SGLT2 (except SGLT2 arm)/TZD
## j) No CVD (I'm using a broad definition: angina, IHD, MI, PAD, revasc, stroke, TIA [as per NICE but with TIA])
## k) No HF before index date
## l) No CKD (stage 3a-5) before index date (assume coded on all, so if missing assume negative)
## m) Exclude if missing QRISK2 or QDHF
### NB: this includes if any required variables missing (smoking status, baseline HbA1c) or if values out of range (age<25 or >84, cholHDL<1 or >11, HbA1c<40 or >150, SBP<70 or >210, BMI<20 - QDHF has not been calculated in these cases))
### will also exclude anyone without QRISK2 score (missing if missing smoking status or age/cholHDL/SBP/BMI outside of range [weirdly cholHDL range for QRISK2 is 1-12 vs 1-11 for QDHF])
## n) Exclude if no IMD - as adjusting for this/using for PS
  

# Use "t2d_1stinstance" cohort_dataset which already has a)-d) applied
# all_drug_periods_dataset is used to define later start dates of meds for censoring


define_cohort <- function(cohort_dataset, all_drug_periods_dataset) {
  
  
  # Keep those aged >=18 and within study period and second line or later (e-h above)
  cohort <- cohort_dataset %>%
    filter(dstartdate_age>=18 &
             (drugclass=="SGLT2" | drugclass=="DPP4" | drugclass=="SU") &
             dstartdate>=as.Date("2013-01-01") &
             drugline_all!=1)
  
  
  # Remove if on insulin at start (i above)
  cohort <- cohort %>%
    filter(INS==0)
  

  # Remove if on GLP1/SGLT2 (except SGLT2 arm)/TZD (i above)
  cohort <- cohort %>%
    filter(GLP1==0 & TZD==0 & (drugclass=="SGLT2" | SGLT2==0))
  
    
  # Remove if CVD before index date (j above)
  cohort <- cohort %>%
    mutate(predrug_cvd=ifelse(predrug_angina==1 | predrug_ihd==1 | predrug_myocardialinfarction==1 | predrug_pad==1 | predrug_revasc==1 | predrug_stroke==1 | predrug_tia==1, 1, 0)) %>%
    filter(predrug_cvd==0)
  
  
  # Remove if HF before index date (k above)
  cohort <- cohort %>%
    filter(predrug_heartfailure==0)
  

  # Remove if CKD before index date (l above)
  cohort <- cohort %>%
    filter(is.na(preckdstage) | (preckdstage!="stage_3a" & preckdstage!="stage_3b" & preckdstage!="stage_4" & preckdstage!="stage_5"))

  
  # Remove if don't have QDHF variables (m above; also removes those without QRISK2)
  cohort <- cohort %>%
    filter(!is.na(qdiabeteshf_5yr_score))
  
  
  # Remove if don't have IMD (n above)
  cohort <- cohort %>%
    filter(!is.na(imd2015_10))
  

  ## Use all DPP4, GLP1, SGLT2, SU + TZD starts to code up later censoring
  
  later_dpp4 <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="DPP4") %>%
                  select(patid, next_dpp4=dstartdate)), by="patid") %>%
    filter(next_dpp4>=dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_dpp4_start=min(next_dpp4, na.rm=TRUE)) %>%
    ungroup()
  
  later_glp1 <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="GLP1") %>%
                  select(patid, next_glp1=dstartdate)), by="patid") %>%
    filter(next_glp1>=dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_glp1_start=min(next_glp1, na.rm=TRUE)) %>%
    ungroup()
  
  later_sglt2 <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="SGLT2") %>%
                  select(patid, next_sglt2=dstartdate)), by="patid") %>%
    filter(next_sglt2>=dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_sglt2_start=min(next_sglt2, na.rm=TRUE)) %>%
    ungroup()
  
  later_su <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="SU") %>%
                  select(patid, next_su=dstartdate)), by="patid") %>%
    filter(next_su>=dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_su_start=min(next_su, na.rm=TRUE)) %>%
    ungroup()
  

  later_tzd <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="TZD") %>%
                  select(patid, next_tzd=dstartdate)), by="patid") %>%
    filter(next_tzd>=dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_tzd_start=min(next_tzd, na.rm=TRUE)) %>%
    ungroup()

    cohort <- cohort %>%
      left_join(later_dpp4, by=c("patid", "dstartdate")) %>%
      left_join(later_glp1, by=c("patid", "dstartdate")) %>%
      left_join(later_sglt2, by=c("patid", "dstartdate")) %>%
      left_join(later_su, by=c("patid", "dstartdate")) %>%
      left_join(later_tzd, by=c("patid", "dstartdate"))
  
  
  # Tidy up variables needed for adjustment/weighting
  ## Also code up death cause variables
  
  cohort <- cohort %>%
    
    mutate(malesex=ifelse(gender==1, 1, 0),
           ethnicity_qrisk2_decoded=case_when(is.na(ethnicity_qrisk2) ~"missing",
                                              ethnicity_qrisk2==1 ~"White",
                                              ethnicity_qrisk2==2 ~"Indian",
                                              ethnicity_qrisk2==3 ~"Pakistani",
                                              ethnicity_qrisk2==4 ~"Bangladeshi",
                                              ethnicity_qrisk2==5 ~"Other Asian",
                                              ethnicity_qrisk2==6 ~"Black Caribbean",
                                              ethnicity_qrisk2==7 ~"Black African",
                                              ethnicity_qrisk2==8 ~"Chinese",
                                              ethnicity_qrisk2==9 ~"Other"),
           imd2015_10=as.factor(imd2015_10),
           qrisk2_smoking_cat=as.factor(qrisk2_smoking_cat),
           
           drugline_all=as.factor(ifelse(drugline_all>=5, 5, drugline_all)),
           ncurrtx=as.factor(DPP4+GLP1+MFN+SU+SGLT2+TZD+INS),          #INS, GLP1 and TZD should be 0 but include anyway; ignore Acarbose and Glinide
           initiation_year=as.factor(format(dstartdate,"%Y")),
           
           statin_last_6_months=ifelse(!is.na(predrug_latest_statins) & as.numeric(difftime(dstartdate, predrug_latest_statins, unit="days"))<=183, 1, 0),
           acei_last_6_months=ifelse(!is.na(predrug_latest_ace_inhibitors) & as.numeric(difftime(dstartdate, predrug_latest_ace_inhibitors, unit="days"))<=183, 1, 0),
           arb_last_6_months=ifelse(!is.na(predrug_latest_arb) & as.numeric(difftime(dstartdate, predrug_latest_arb, unit="days"))<=183, 1, 0),
           bb_last_6_months=ifelse(!is.na(predrug_latest_beta_blockers) & as.numeric(difftime(dstartdate, predrug_latest_beta_blockers, unit="days"))<=183, 1, 0),
           ccb_last_6_months=ifelse(!is.na(predrug_latest_calcium_channel_blockers) & as.numeric(difftime(dstartdate, predrug_latest_calcium_channel_blockers, unit="days"))<=183, 1, 0),
           thiazide_diuretics_last_6_months=ifelse(!is.na(predrug_latest_thiazide_diuretics) & as.numeric(difftime(dstartdate, predrug_latest_thiazide_diuretics, unit="days"))<=183, 1, 0),
           loop_diuretics_last_6_months=ifelse(!is.na(predrug_latest_loop_diuretics) & as.numeric(difftime(dstartdate, predrug_latest_loop_diuretics, unit="days"))<=183, 1, 0),
           ksparing_diuretics_last_6_months=ifelse(!is.na(predrug_latest_ksparing_diuretics) & as.numeric(difftime(dstartdate, predrug_latest_ksparing_diuretics, unit="days"))<=183, 1, 0),
           
           hosp_admission_prev_year=as.factor(hosp_admission_prev_year),
           hosp_admission_prev_year_count=as.factor(ifelse(hosp_admission_prev_year_count==0, 0,
                                                           ifelse(hosp_admission_prev_year_count<=2, 1, 2))),
           cv_death_date_any_cause=if_else(!is.na(death_date) & !is.na(cv_death_any_cause) & cv_death_any_cause==1, death_date, as.Date(NA)),
           cv_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(cv_death_primary_cause) & cv_death_primary_cause==1, death_date, as.Date(NA)),
           hf_death_date_any_cause=if_else(!is.na(death_date) & !is.na(hf_death_any_cause) & hf_death_any_cause==1, death_date, as.Date(NA)),
           hf_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(hf_death_primary_cause) & hf_death_primary_cause==1, death_date, as.Date(NA)))
  #note re death_date not always being from ONS death
  
  return(cohort)

}