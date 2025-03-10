
# Produce survival variables for all endpoints (including for sensitivity analysis)
## All censored at 5 years post drug start / end of GP records / death
## Additional censoring for sensitivity analyses

# Main analysis:
## 'mace': stroke, MI, CV death
## 'expanded_mace': stroke, MI, CV death, revasc, HES unstable angina
## 'hf'
## 'hosp': all-cause hospitalisation
## 'death': all-cause mortality

# Sensitivity analysis:
## 'narrow_mace': hospitalisation for incident MI (subset of HES codes), incident stroke (subset of HES codes, includes ischaemic only), CV death - all as primary cause for hospitalisation/death only
## 'narrow_hf': hospitalisation or death with HF as primary cause
## Also includes code for alternative censoring strategies


add_surv_vars <- function(cohort_dataset, main_only=FALSE) {
  
  # Add survival variables for outcomes for main analysis
  main_outcomes <- c("mace", "expanded_mace", "hf", "hosp", "death")
  
  cohort <- cohort_dataset %>%
    
    mutate(cens_main=pmin(dstartdate+(365.25*5),
                         gp_record_end,
                         death_date,
                         next_tzd_start,
                         next_glp1_start,
                         if_else(drugclass!="SGLT2", next_sglt2_start, as.Date("2050-01-01")),
                         if_else(drugclass!="DPP4", dpp4_start_later, as.Date("2050-01-01")),
                         if_else(drugclass!="SU", su_start_later, as.Date("2050-01-01")),
                         na.rm=TRUE),
           
           mace_outcome=pmin(postdrug_first_myocardialinfarction,
                             postdrug_first_stroke,
                             cv_death_date_any_cause,
                             na.rm=TRUE),
          
           expanded_mace_outcome=pmin(postdrug_first_myocardialinfarction,
                                      postdrug_first_stroke,
                                      cv_death_date_any_cause,
                                      postdrug_first_revasc,
                                      postdrug_first_unstableangina,
                                      na.rm=TRUE),
             
           hf_outcome=pmin(postdrug_first_heartfailure,
                           hf_death_date_any_cause,
                           na.rm=TRUE),
           
           hosp_outcome=postdrug_first_emergency_hosp,
           
           death_outcome=death_date)
  
  
  for (i in main_outcomes) {

    outcome_var=paste0(i, "_outcome")
    censdate_var=paste0(i, "_censdate")
    censvar_var=paste0(i, "_censvar")
    censtime_var=paste0(i, "_censtime_yrs")

    cohort <- cohort %>%
      mutate({{censdate_var}}:=pmin(!!sym(outcome_var), cens_main, na.rm=TRUE),
             {{censvar_var}}:=ifelse(!is.na(!!sym(outcome_var)) & !!sym(censdate_var)==!!sym(outcome_var), 1, 0),
             {{censtime_var}}:=as.numeric(difftime(!!sym(censdate_var), dstartdate, unit="days"))/365.25)
    
    
  }
  
  if (main_only==TRUE) {
    message(paste("survival variables for", paste(main_outcomes, collapse=", "), "added"))
    }
  
  
  # Add survival variables for outcomes for sensitivity analyses
 
  else {
    
    sensitivity_outcomes <- c("narrow_mace", "narrow_hf")
    
    cohort <- cohort %>%
      
      mutate(cens_itt=pmin(dstartdate+(365.25*5),
                           gp_record_end,
                           death_date,
                           na.rm=TRUE),
             
             cens_pp=pmin(dstartdate+(365.25*5),
                          gp_record_end,
                          death_date,
                          dstartdate+timetochange, #dstopdate can be before dstartdate+timetochange as timetochange is time to new drug combo, where dstopdate is last prescription - have used timetochange instead - max 6 months
                          na.rm=TRUE),
             
             narrow_mace_outcome=pmin(postdrug_first_incident_mi,
                                      postdrug_first_incident_stroke,
                                      cv_death_date_primary_cause,
                                      na.rm=TRUE),
             
             narrow_hf_outcome=pmin(postdrug_first_primary_hhf,
                                    hf_death_date_primary_cause,
                                    na.rm=TRUE),
             
             hf_itt_censdate=pmin(hf_outcome, cens_itt, na.rm=TRUE),
             hf_itt_censvar=ifelse(!is.na(hf_outcome) & hf_itt_censdate==hf_outcome, 1, 0),
             hf_itt_censtime_yrs=as.numeric(difftime(hf_itt_censdate, dstartdate, unit="days"))/365.25,
             
             hf_pp_censdate=pmin(hf_outcome, cens_pp, na.rm=TRUE),
             hf_pp_censvar=ifelse(!is.na(hf_outcome) & hf_pp_censdate==hf_outcome, 1, 0),
             hf_pp_censtime_yrs=as.numeric(difftime(hf_pp_censdate, dstartdate, unit="days"))/365.25,
             
             mace_itt_censdate=pmin(mace_outcome, cens_itt, na.rm=TRUE),
             mace_itt_censvar=ifelse(!is.na(mace_outcome) & mace_itt_censdate==mace_outcome, 1, 0),
             mace_itt_censtime_yrs=as.numeric(difftime(mace_itt_censdate, dstartdate, unit="days"))/365.25,
             
             mace_pp_censdate=pmin(mace_outcome, cens_pp, na.rm=TRUE),
             mace_pp_censvar=ifelse(!is.na(mace_outcome) & mace_pp_censdate==mace_outcome, 1, 0),
             mace_pp_censtime_yrs=as.numeric(difftime(mace_pp_censdate, dstartdate, unit="days"))/365.25)
    
    for (i in sensitivity_outcomes) {
      
      outcome_var=paste0(i, "_outcome")
      censdate_var=paste0(i, "_censdate")
      censvar_var=paste0(i, "_censvar")
      censtime_var=paste0(i, "_censtime_yrs")
      
      cohort <- cohort %>%
        mutate({{censdate_var}}:=pmin(!!sym(outcome_var), cens_itt, na.rm=TRUE),
               {{censvar_var}}:=ifelse(!is.na(!!sym(outcome_var)) & !!sym(censdate_var)==!!sym(outcome_var), 1, 0),
               {{censtime_var}}:=as.numeric(difftime(!!sym(censdate_var), dstartdate, unit="days"))/365.25)
      }

    if (main_only==FALSE) {
      message(paste0("survival variables for ", paste0(main_outcomes, collapse=", "), ", ", paste0(unlist(sensitivity_outcomes), collapse=", "), " added"))
      }
    
    }
   
return(cohort) 
  
  }
  