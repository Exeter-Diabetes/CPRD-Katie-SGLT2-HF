
# Takes random 20% of DPP4SU arm and re-estimates baseline hazard to outcome, then recalibrates rest of cohort
## (Except 40% eGFR decline risk score which is logistic regression: need to add in) - also need to add in other CKD risk score


calibrate_risk_score <- function(cohort_dataset, risk_score, outcome) {
  
  timevar=paste0(outcome, "_censtime_yrs")
  statusvar=paste0(outcome, "_censvar")
  
  
  ## Assign random 20% of DPP4SU arm as calibration cohort and remove from main cohort
  set.seed(123)
  
  cal_cohort <- cohort_dataset %>%
    filter(studydrug=="DPP4SU") %>%
    slice_sample(prop=0.2)
  
  cal_cohort_count <- count(cal_cohort)
  
  print(paste0("calibration subset: n=", cal_cohort_count, " will be removed from cohort dataset"))
  
  noncal_cohort <- cohort_dataset %>%
    anti_join(cal_cohort, by=c("patid", "dstartdate", "studydrug"))


  if (risk_score=="qrisk2") {
    
    # C-statistic in whole cohort for uncalibrated score
    cohort_dataset <- cohort_dataset %>%
      mutate(qrisk2_survival=(100-qrisk2_5yr_score)/100)
    
    surv_mod <- coxph(as.formula(paste0("Surv(", timevar, ", ", statusvar, ")~qrisk2_survival")), data=cohort_dataset, method="breslow")
    cstat <- round(summary(surv_mod)$concordance[1], 3)
    cstat_lower <- round(summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2]), 3)
    cstat_upper <- round(summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2]), 3)
    
    print(paste0("C-statistic for uncalibrated score in whole cohort: ", cstat, " (", cstat_lower, ", ", cstat_upper, ")"))
    
    
    # Recalibrate in calibration cohort
    cal_females <- cal_cohort %>% filter(malesex==0)
    cal_males <- cal_cohort %>% filter(malesex==1)

    female_recal_mod <- coxph(as.formula(paste0("Surv(", timevar, ", ", statusvar, ")~offset(qrisk2_lin_predictor)")), data=cal_females)
    female_qrisk2_surv <- summary(survfit(female_recal_mod),time=5)$surv

    male_recal_mod <- coxph(as.formula(paste0("Surv(", timevar, ", ", statusvar, ")~offset(qrisk2_lin_predictor)")), data=cal_males)
    male_qrisk2_surv <- summary(survfit(male_recal_mod),time=5)$surv

    
    # Recalibrate non-calibration cohort
    noncal_cohort <- noncal_cohort %>%
      group_by(malesex) %>%
      mutate(centred_qrisk2_lin_predictor=qrisk2_lin_predictor-mean(qrisk2_lin_predictor)) %>%
      ungroup() %>%
      mutate(qrisk2_5yr_score_cal=ifelse(malesex==1, (1-(male_qrisk2_surv^exp(centred_qrisk2_lin_predictor)))*100, (1-(female_qrisk2_surv^exp(centred_qrisk2_lin_predictor)))*100))

    print(paste0("Mean linear predictor male: ", noncal_cohort %>% filter(malesex==1) %>% summarise(mean_qrisk2_lin_predictor=mean(qrisk2_lin_predictor))))
    print(paste0("Mean linear predictor female: ", noncal_cohort %>% filter(malesex==0) %>% summarise(mean_qrisk2_lin_predictor=mean(qrisk2_lin_predictor))))
    print(paste0("New male surv: ", male_qrisk2_surv))
    print(paste0("New female surv: ", female_qrisk2_surv))
    
    
    # C-statistic of calibrated score in non-calibration cohort
    noncal_cohort <- noncal_cohort %>%
      mutate(qrisk2_survival_cal=(100-qrisk2_5yr_score_cal)/100)

    surv_mod <- coxph(as.formula(paste0("Surv(", timevar, ", ", statusvar, ")~qrisk2_survival_cal")), data=noncal_cohort, method="breslow")
    cstat <- round(summary(surv_mod)$concordance[1], 3)
    cstat_lower <- round(summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2]), 3)
    cstat_upper <- round(summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2]), 3)

    print(paste0("C-statistic for calibrated score in non-calibration cohort: ", cstat, " (", cstat_lower, ", ", cstat_upper, ")"))
    
  }
  
  
  
  if (risk_score=="qdiabeteshf") {
    
    # C-statistic in whole cohort for uncalibrated score
    cohort_dataset <- cohort_dataset %>%
      mutate(qdiabeteshf_survival=(100-qdiabeteshf_5yr_score)/100)
    
    surv_mod <- coxph(as.formula(paste0("Surv(", timevar, ", ", statusvar, ")~qdiabeteshf_survival")), data=cohort_dataset, method="breslow")
    cstat <- round(summary(surv_mod)$concordance[1], 3)
    cstat_lower <- round(summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2]), 3)
    cstat_upper <- round(summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2]), 3)
    
    print(paste0("C-statistic for uncalibrated score in whole cohort: ", cstat, " (", cstat_lower, ", ", cstat_upper, ")"))
    
    
    # Recalibrate in calibration cohort
    cal_females <- cal_cohort %>% filter(malesex==0)
    cal_males <- cal_cohort %>% filter(malesex==1)
    
    female_recal_mod <- coxph(as.formula(paste0("Surv(", timevar, ", ", statusvar, ")~offset(qdiabeteshf_lin_predictor)")), data=cal_females)
    female_qdiabeteshf_surv <- summary(survfit(female_recal_mod),time=5)$surv
    
    male_recal_mod <- coxph(as.formula(paste0("Surv(", timevar, ", ", statusvar, ")~offset(qdiabeteshf_lin_predictor)")), data=cal_males)
    male_qdiabeteshf_surv <- summary(survfit(male_recal_mod),time=5)$surv
    
    
    # Recalibrate non-calibration cohort
    noncal_cohort <- noncal_cohort %>%
      group_by(malesex) %>%
      mutate(centred_qdiabeteshf_lin_predictor=qdiabeteshf_lin_predictor-mean(qdiabeteshf_lin_predictor)) %>%
      ungroup() %>%
      mutate(qdiabeteshf_5yr_score_cal=ifelse(malesex==1, (1-(male_qdiabeteshf_surv^exp(centred_qdiabeteshf_lin_predictor)))*100, (1-(female_qdiabeteshf_surv^exp(centred_qdiabeteshf_lin_predictor)))*100))
    
    print(paste0("Mean linear predictor male: ", noncal_cohort %>% filter(malesex==1) %>% summarise(mean_qdiabeteshf_lin_predictor=mean(qdiabeteshf_lin_predictor))))
    print(paste0("Mean linear predictor female: ", noncal_cohort %>% filter(malesex==0) %>% summarise(mean_qdiabeteshf_lin_predictor=mean(qdiabeteshf_lin_predictor))))
    print(paste0("New male surv: ", male_qdiabeteshf_surv))
    print(paste0("New female surv: ", female_qdiabeteshf_surv))
    
    
    # C-statistic of calibrated score in non-calibration cohort
    noncal_cohort <- noncal_cohort %>%
      mutate(qdiabeteshf_survival_cal=(100-qdiabeteshf_5yr_score_cal)/100)
    
    surv_mod <- coxph(as.formula(paste0("Surv(", timevar, ", ", statusvar, ")~qdiabeteshf_survival_cal")), data=noncal_cohort, method="breslow")
    cstat <- round(summary(surv_mod)$concordance[1], 3)
    cstat_lower <- round(summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2]), 3)
    cstat_upper <- round(summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2]), 3)
    
    print(paste0("C-statistic for calibrated score in non-calibration cohort: ", cstat, " (", cstat_lower, ", ", cstat_upper, ")"))
    
  }
  
  return(noncal_cohort)

}