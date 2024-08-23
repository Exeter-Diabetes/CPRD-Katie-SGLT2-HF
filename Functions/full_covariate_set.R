

return_cov_set <- function(which_cov_set){
  
  if (which_cov_set=="weight") {
  
  cov_set <- "malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score"
  
  }
  
  if (which_cov_set=="adjust") {
  
  cov_set <- "malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score"
  
  }
  
  return(cov_set)

}
