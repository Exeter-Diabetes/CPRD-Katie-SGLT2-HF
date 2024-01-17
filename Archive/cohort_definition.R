
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
             drugline_all!=1) %>%
    mutate(studydrug=ifelse(drugclass=="SGLT2", "SGLT2", ifelse(drugclass=="SU", "SU", "DPP4")))
  
  
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
    mutate(predrug_heartfailure==0)
  

  # Remove if CKD before index date (l above)
  cohort <- cohort %>%
    filter(is.na(preckdstage) | (preckdstage!="stage_3a" & preckdstage!="stage_3b" & preckdstage!="stage_4" & preckdstage!="stage_5"))

  
  # Remove if don't have QDHF variables (m above; also removes those without QRISK2)
  cohort <- cohort %>%
    filter(!is.na(qdiabeteshf_5yr_score))
  
  
  # Remove if don't have IMD (n above)
  cohort <- cohort %>%
    filter(!is.na(imd2015_10))
  

  ## Use all SGLT2, GLP1 + TZD starts to code up later censoring
  ### Also get latest SGLT2 stop dates - may want to do sensitivity analysis where SU arm have to have not taken SGLT2 in the previous year (as per some papers)
  
  later_sglt2 <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="SGLT2") %>%
                  select(patid, next_sglt2=dstartdate)), by="patid") %>%
    filter(next_sglt2>dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_sglt2_start=min(next_sglt2, na.rm=TRUE)) %>%
    ungroup()
  
  
  later_glp1 <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="GLP1") %>%
                  select(patid, next_glp1=dstartdate)), by="patid") %>%
    filter(next_glp1>dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_glp1_start=min(next_glp1, na.rm=TRUE)) %>%
    ungroup()
  
  
  later_tzd <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="TZD") %>%
                  select(patid, next_tzd=dstartdate)), by="patid") %>%
    filter(next_tzd>dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_tzd_start=min(next_tzd, na.rm=TRUE)) %>%
    ungroup()
  
  
  last_sglt2_stop <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
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
  
  
  # Tidy up gender, ncurrtx, drugline and ethnicity variables
  ## Also code up death cause variables
  
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
           hf_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(hf_death_primary_cause) & hf_death_primary_cause==1, death_date, as.Date(NA)))
  #note re death_date not always being from ONS death
  
  return(cohort)

}