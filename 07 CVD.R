
# At the moment have removed QRISK2 calibration - so weighting and adjustment using uncalibrated score (shouldn't make much difference although does make a little as men and women calibrated separately and differently)

# HRs calculated via different methods

# HR plot by CVD risk (uncalibrated QRISK2)


############################################################################################

# Setup
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(gridExtra)
library(grid)
library(PSweight)
library(rms)
library(ggthemes)
library(forestplot)
library(tidycmprsk)
library(table1)
library(cowplot)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# 1 Cohort selection and variable setup

setwd("/slade/CPRD_data/Katie SGLT2/Processed data")
load("treatment_outcome_cohort_jun24.rda")
#169,041

table(cohort$studydrug)
# DPP4SU 111673   
# SGLT2 57368  


## B Make variables for survival analysis of all endpoints (see survival_variables function for details)

setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions")
source("survival_variables.R")

cohort <- add_surv_vars(cohort, main_only=FALSE)


## Make sure don't accidentally use HF variables
cohort <- cohort %>% select(-c(hf_censtime_yrs, narrow_hf_censtime_yrs, qdiabeteshf_5yr_score))
         

############################################################################################

# 2 Table of weighted characteristics

cat <- function(x, ...) {
  c("", sapply(stats.apply.rounding(stats.default(x)), function(y) with(y, sprintf("%s (%s%%)", prettyNum(FREQ, big.mark=","), PCT))))
}

cont <- function(x) {
  with(stats.apply.rounding(stats.default(x)), c("Median (IQR)"=sprintf("%s (%s-%s)", round_pad(as.numeric(MEDIAN),1), round_pad(as.numeric(Q1),1), round_pad(as.numeric(Q3),1))))
}

missing <- function(x, ...) {
  with(stats.apply.rounding(stats.default(x)), c("Missing"=sprintf("%s", prettyNum(NMISS, big.mark=","))))
}

strat <- function (label, n, ...) {
  sprintf("<span class='stratlabel'>%s<br><span class='stratn'>(N=%s)</span></span>", 
          label, prettyNum(n, big.mark=","))
}

rndr <- function(x, name, ...) {
  y <- render.default(x, name, ...)
  if (is.logical(x)) {
    y[2]
  } else {
    y
  }
}

# Reduce ethnicity to 5 category, smoking to 3 category, and IMD to quintiles for this

cohort <- cohort %>%
  mutate(ethnicity_decoded=factor(case_when(ethnicity_qrisk2_decoded=="missing" ~"Missing",
                                            ethnicity_qrisk2_decoded=="White" ~"White",
                                            ethnicity_qrisk2_decoded=="Indian" | ethnicity_qrisk2_decoded=="Pakistani" |  ethnicity_qrisk2_decoded=="Bangladeshi" | ethnicity_qrisk2_decoded=="Other Asian" ~"Asian",
                                            ethnicity_qrisk2_decoded=="Black Caribbean" | ethnicity_qrisk2_decoded=="Black African" ~"Black",
                                            ethnicity_qrisk2_decoded=="Chinese" | ethnicity_qrisk2_decoded=="Other" ~"Other"), levels=c("White", "Asian", "Black", "Other", "Missing")),
         smoker_decoded=factor(case_when(is.na(qrisk2_smoking_cat) ~as.character(NA),
                                         qrisk2_smoking_cat==0 ~"Non-smoker",
                                         qrisk2_smoking_cat==1 ~"Ex-smoker",
                                         qrisk2_smoking_cat==2 | qrisk2_smoking_cat==3 | qrisk2_smoking_cat==4 ~"Active smoker"), levels=c("Non-smoker", "Active smoker", "Ex-smoker")),
         imd_quintiles=as.factor(case_when(is.na(imd2015_10) ~as.character(NA),
                                           imd2015_10==1 | imd2015_10==2 ~"1 (least deprived)",
                                           imd2015_10==3 | imd2015_10==4 ~"2",
                                           imd2015_10==5 | imd2015_10==6 ~"3",
                                           imd2015_10==7 | imd2015_10==8 ~"4",
                                           imd2015_10==9 | imd2015_10==10 ~"5 (most deprived)")))


# Add labels
label(cohort$malesex)               <- "Sex (% male)"
label(cohort$dstartdate_age)        <- "Age at drug initiation (years)"
label(cohort$dstartdate_dm_dur_all) <- "Diabetes duration (years)"
label(cohort$ethnicity_decoded)     <- "Ethnicity"
label(cohort$imd_quintiles)         <- "Index of Multiple Deprivation quintile"
label(cohort$smoker_decoded)        <- "Smoking status"
label(cohort$hypertension)          <- "Hypertension"
label(cohort$predrug_af)            <- "Atrial fibrillation"
label(cohort$hosp_admission_prev_year_count) <- "Number of hospital admissions in previous year"
label(cohort$prebmi)                <- "BMI (kg/m2)"
label(cohort$prehba1c2yrs)          <- "HbA1c (mmol/mol)"
label(cohort$presbp)                <- "SBP (mmHg)"
label(cohort$precholhdl)            <- "Cholesterol:HDL"
label(cohort$drugline_all)          <- "Drug line"
label(cohort$ncurrtx_cat)           <- "Number of other current non-insulin glucose-lowering medications"
label(cohort$INS)                   <- "Current insulin use"
label(cohort$initiation_year)       <- "Year of drug initiation"
label(cohort$qrisk2_5yr_score) <- "QRISK2 5-year score (%)"


setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions")
source("full_covariate_set.R")

return_cov_set("weight")
## replace qdhf with qrisk2
"malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + INS + initiation_year + qrisk2_5yr_score"


# Table with overlap weighting
# Did initially manage to code up for proportions, but easier to do this way for medians

ps.formula <- formula("studydrug ~ malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + INS + initiation_year + qrisk2_5yr_score")

overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort), weight="overlap")
cohort$overlap_weight <- overlap$ps.weights$overlap
rm(overlap)
gc()

weighted_cohort <- cohort %>%
  select(studydrug, malesex, dstartdate_age, dstartdate_dm_dur_all, ethnicity_decoded, imd_quintiles, smoker_decoded, hypertension, predrug_af, hosp_admission_prev_year_count, prebmi, prehba1c2yrs, presbp, precholhdl, drugline_all, ncurrtx_cat, INS, initiation_year, qrisk2_5yr_score, overlap_weight) %>%
  mutate(count=round(overlap_weight*10000000),0) %>% 
  uncount(count)


strat <- function (label, n, ...) {
  sprintf("<span class='stratlabel'>%s</span>", 
          label, prettyNum(n, big.mark=","))
}

cat <- function(x, ...) {
  c("", sapply(stats.apply.rounding(stats.default(x)), function(y) with(y, sprintf("%s%%", PCT))))
}


# Same variables as above
table1(~ malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_decoded + imd_quintiles + smoker_decoded + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + precholhdl + drugline_all + ncurrtx_cat + INS + initiation_year + qrisk2_5yr_score | studydrug, data=weighted_cohort, overall=F, render=rndr, render.categorical=cat, render.continuous=cont, render.strat=strat)
#need to remove missing counts from Chol:HDL afterwards
#superscript kg/m2



############################################################################################

# 3 Calculate hazard ratios -with all sensitivity analyses

setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions")
source("full_covariate_set.R")
# now using functions to define covariates adjusted and weighted for

return_cov_set("weight")
## replace qdhf with uncalibrated qrisk2
weight <- "malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + INS + initiation_year + qrisk2_5yr_score"

return_cov_set("adjust")
## replace qdhf with uncalibrated qrisk2
adjust <- "malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year + qrisk2_5yr_score"



## Overlap and IPTW weighting

ps.formula <- formula(paste0("studydrug ~ ", weight))

f_adjusted <- as.formula(paste0("Surv(mace_censtime_yrs, mace_censvar) ~ studydrug  +  ", adjust))

f_unadjusted <- as.formula(paste0("Surv(mace_censtime_yrs, mace_censvar) ~ studydrug"))


## Main analysis

overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(cohort), weight = c("IPW", "overlap"))

cohort$overlap_weights <- overlap$ps.weights$overlap
cohort$iptw_weights <- overlap$ps.weights$IPW

counts <- cohort %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(mace_censvar),
            total_time=sum(mace_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

adjoverlap <- coxph(f_adjusted, cohort, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="overlap weight")



## No DPP4i

no_dpp4i <- cohort %>% filter(drugclass!="DPP4")

overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(no_dpp4i), weight = "overlap")

no_dpp4i$overlap_weights <- overlap$ps.weights$overlap

counts_no_dpp4i <- no_dpp4i %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(mace_censvar),
            total_time=sum(mace_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))


adjoverlap_no_dpp4i <- coxph(f_adjusted, no_dpp4i, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="overlap weight no dpp4i")


## No SU

no_su <- cohort %>% filter(drugclass!="SU")

overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(no_su), weight = "overlap")

no_su$overlap_weights <- overlap$ps.weights$overlap

counts_no_su <- no_su %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(mace_censvar),
            total_time=sum(mace_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))


adjoverlap_no_su <- coxph(f_adjusted, no_su, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="overlap weight no su")



## Different weighting/adjustment


unadjoverlap <- coxph(f_unadjusted, cohort, weights=iptw_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="unadjusted overlap weight")

adjusted <- coxph(f_adjusted, cohort) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="no weight")

adjiptw <- coxph(f_adjusted, cohort, weights=iptw_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="iptw weight")

all_counts <- counts %>%
  bind_rows(counts_no_dpp4i) %>%
  bind_rows(counts_no_su) %>%
  bind_rows(counts) %>%
  bind_rows(counts) %>%
  bind_rows(counts)

all_hrs <- adjoverlap %>%
  bind_rows(adjoverlap_no_dpp4i) %>%
  bind_rows(adjoverlap_no_su) %>%
  bind_rows(unadjoverlap) %>%
  bind_rows(adjusted) %>%
  bind_rows(adjiptw)



## Narrow MACE outcome

f_adjusted <- as.formula(paste0("Surv(narrow_mace_censtime_yrs, narrow_mace_censvar) ~  studydrug  +  ", adjust))

counts <- cohort %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(narrow_mace_censvar),
            total_time=sum(narrow_mace_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

adjoverlap <- coxph(f_adjusted, cohort, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="narrow mace")
#Gives log likelihood error but coefficients seem fine - none are infinite

all_counts <- all_counts %>% bind_rows(counts)

all_hrs <- all_hrs %>% bind_rows(adjoverlap)


## Death as competing risk

cohort <- cohort %>% 
  mutate(mace_censvar.cr=as.factor(ifelse(mace_censvar==0 & !is.na(death_date) & death_date==mace_censdate, 2, mace_censvar)))

counts <- cohort %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(mace_censvar.cr==1),
            total_time=sum(mace_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

# death counts
cohort %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(mace_censvar.cr==2),
            total_time=sum(mace_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

# Adjusted with overlap weighting
## Very slow to run

f_adjusted <- as.formula(paste0("Surv(mace_censtime_yrs, mace_censvar.cr) ~  studydrug  +  ", adjust))

adjoverlap <- crr(f_adjusted, data = cohort, weights=overlap) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="death competing")


all_counts <- all_counts %>% bind_rows(counts)

all_hrs <- all_hrs %>% bind_rows(adjoverlap)



## Different followup

cohort %>% group_by(patid) %>% summarise(count_dpp4=sum(drugclass=="DPP4"), count_su=sum(drugclass=="SU")) %>% filter(count_dpp4==1 & count_su==1) %>% ungroup() %>% count()
#15626 patids with both DPP4 and SU

cohort %>% filter(studydrug=="DPP4SU") %>% group_by(patid) %>% filter(n()>1 & max(dstartdate)==min(dstartdate)) %>% ungroup() %>% distinct(patid) %>% count()
#Includes 268 patids which start DPP4 and SU on same day - have checked and follow up time is 0

169041-15626+268 #153683

cohort %>%
  group_by(patid, studydrug) %>%
  filter(dstartdate==min(dstartdate)) %>%
  ungroup() %>%
  count()
#153683 (haven't eliminated if start on same day)


169041-15626 #153415

sensitivity_cohort <- cohort %>%
  group_by(patid, studydrug) %>%
  filter(dstartdate==min(dstartdate) & (n()==1 | dstartdate!=max(dstartdate) | drugclass=="DPP4")) %>%
  ungroup()
#153415

overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(sensitivity_cohort), weight = c("IPW", "overlap"))
sensitivity_cohort$overlap_weights <- overlap$ps.weights$overlap

f_adjusted <- as.formula(paste0("Surv(mace_can_overlap_censtime_yrs, mace_can_overlap_censvar) ~  studydrug  +  ", adjust))

counts <- sensitivity_cohort %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(mace_can_overlap_censvar),
            total_time=sum(mace_can_overlap_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

adjoverlap <- coxph(f_adjusted, sensitivity_cohort, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="no control duplicate IDs")

all_counts <- all_counts %>% bind_rows(counts)

all_hrs <- all_hrs %>% bind_rows(adjoverlap)




## Subset on insulin

on_ins <- cohort %>% filter(INS==1)

return_cov_set("weight")
# remove insulin

ps.formula <- formula(paste("studydrug ~ malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + initiation_year + qrisk2_5yr_score"))

return_cov_set("adjust")
# remove insulin

f_adjusted <- as.formula(Surv(mace_censtime_yrs, mace_censvar) ~  studydrug  +  malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + initiation_year + qrisk2_5yr_score)


overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(on_ins), weight = c("IPW", "overlap"))

on_ins$overlap_weights <- overlap$ps.weights$overlap

counts <- on_ins %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(mace_censvar),
            total_time=sum(mace_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

adjoverlap <- coxph(f_adjusted, on_ins, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="on insulin")
#log likelihoods converge but looks ok


all_counts <- all_counts %>% bind_rows(counts)

all_hrs <- all_hrs %>% bind_rows(adjoverlap)




# All methods

hrs_to_plot <- cbind(all_counts, all_hrs) %>%
  rbind(data.frame(DPP4SU_count=NA, DPP4SU_ir=NA, SGLT2_count=NA, SGLT2_ir=NA, analysis="trial", estimate=0.94, conf.low=0.83, conf.high=1.07))

setwd("/slade/CPRD_data/Katie SGLT2/Processed data")

save(hrs_to_plot, file="hr_data_cvd.Rda")

load("hr_data_cvd.Rda")





# Format for plot

text <- hrs_to_plot %>%
  mutate(hr_text=paste0(round_pad(estimate, 2), " (", round_pad(conf.low, 2), "-", round_pad(conf.high, 2), ")")) %>%
  select(analysis, DPP4SU_count, DPP4SU_ir, SGLT2_count, SGLT2_ir, mean=estimate, lower=conf.low, upper=conf.high, hr_text)


tiff("/slade/CPRD_data/Katie SGLT2/Plots/HR_methods_CVD.tiff", width=18, height=11, units = "in", res=400) 

styles <- fpShapesGp(
  lines = list(
    gpar(col = "black"),
    gpar(col = "black"),
    gpar(col = "black"),
    gpar(col = "black"),
    gpar(col = "black"),
    gpar(col = "black"),
    gpar(col = "black"),
    gpar(col = "black"),
    gpar(col = "black"),
    gpar(col = "black"),
    gpar(col = "black"),
    gpar(col = "red")
  ),
  box = list(
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "red", col="red")
  )
)


text %>% forestplot(labeltext = list(analysis_text=list("Primary analysis: adjusted\nwith overlap weighting", "SU-only comparator group", "DPP4i-only comparator group", "Unadjusted with overlap weighting", "Adjusted with no weighting", "Adjusted with IPTW", expression("More specific MACE definition"^1), expression("Patients on either DPP4i or SU"^2), expression("Death as a competing risk"^3), "Insulin-treated subgroup only", "Trial meta-analysis"), SGLT2_count=text$SGLT2_count, SGLT2_ir=text$SGLT2_ir, DPP4SU_count=text$DPP4SU_count, DPP4SU_ir=text$DPP4SU_ir, hr_text=text$hr_text),
                    xlog = T,
                    ci.vertices = TRUE,
                    xlab = "Hazard ratio ",
                    fn.ci_norm = fpDrawCircleCI,
                    xticks = c(log(0.6), log(0.7), log(0.8), log(0.9), log(1.0), log(1.1), log(1.2)),
                    boxsize = 0.25,
                    align = c("l", "l", "l", "l", "l"),
                    graph.pos=6,
                    lwd.zero = 1.2,
                    shapes_gp = styles,
                    col = fpColors(zero = c("grey50"), text=c('black','black','black','black','black','black','black','black','black','black', 'black', 'red')),
                    lwd.ci=1,
                    ci.vertices.height = 0.2,
                    colgap = unit(6, "mm"),
                    txt_gp = fpTxtGp(xlab=gpar(cex=1.8, fontface="bold"), ticks=gpar(cex=1.7), label=gpar(cex=1.7))) %>%
  fp_add_header(analysis_text = "", SGLT2_count = "SGLT2i\nEvents", SGLT2_ir = "\nIR", DPP4SU_count = "Comparator\nEvents", DPP4SU_ir = "\nIR", hr_text = "Hazard ratio\n(95% CI)")


dev.off()


############################################################################################

# 4 HR by CVD risk

## Keep variables of interest only so don't accidentally use HF ones
cohort <- cohort %>%
  select(patid, studydrug, malesex, dstartdate_age, dstartdate_dm_dur_all, ethnicity_qrisk2_decoded, imd2015_10, qrisk2_smoking_cat, hypertension, predrug_af, hosp_admission_prev_year_count, prebmi, prehba1c2yrs, presbp, drugline_all, ncurrtx_cat, INS, initiation_year, qrisk2_5yr_score, mace_censtime_yrs, mace_censvar)


# Overall HR:

setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions")
source("full_covariate_set.R")

return_cov_set("weight")
"malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score"
# replace qdiabeteshf_5yr_score with qrisk2_5yr_score_cal

ps.formula <- formula("studydrug ~ malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + INS + initiation_year + qrisk2_5yr_score")

overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(cohort), weight = c("IPW", "overlap"))

cohort$overlap_weights <- overlap$ps.weights$overlap

return_cov_set("adjust")
"malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score"

f_adjusted <- as.formula("Surv(mace_censtime_yrs, mace_censvar) ~ studydrug  +  malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year + qrisk2_5yr_score")

coxph(f_adjusted, cohort, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high)
#0.84 (0.77-0.91)




# Adjusted plot
## Adjusting for same as previous except not cal QRISK2 - put in as interaction term with studydrug

cohort <- cohort %>% mutate(studydrug=as.vector(studydrug))

# Remove constant variables otherwise datadist has issues

ddist <- datadist(cohort)
options(datadist='ddist')

setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions")
source("full_covariate_set.R")
return_cov_set("adjust")
#"malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score"
#remove qdiabeteshf_5yr_score and add qrisk2 as interaction term


model <- cph(Surv(mace_censtime_yrs, mace_censvar) ~  studydrug*rcs(qrisk2_5yr_score,3) + malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year, data=cohort, x=T, y=T)

anova(model)
# no interaction - good (p=0.12)

c1 <- quantile(cohort$qrisk2_5yr_score, .01, na.rm=TRUE)
c99 <- quantile(cohort$qrisk2_5yr_score, .99, na.rm=TRUE)


contrast_spline <- contrast(model, list(studydrug = "SGLT2", qrisk2_5yr_score = seq(c1, c99, by=0.05)), list(studydrug = "DPP4SU", qrisk2_5yr_score = seq(c1, c99, by=0.05)))

contrast_spline_df <- as.data.frame(contrast_spline[c('qrisk2_5yr_score','Contrast','Lower','Upper')])



# Add histogram
marginal_distribution <- function(x,var) {
  ggplot(x, aes_string(x = var)) +
    geom_histogram(bins = 64, alpha = 0.4, position = "identity") +
    guides(fill = FALSE) +
    theme(legend.title = element_blank(), panel.background = element_rect( fill = "white",color = "grey50")) +
    scale_x_continuous(breaks = seq(0,30,5)) +
    xlab(expression(paste("5-year absolute cardiovascular disease risk (%)"))) +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(color="grey50"),
          text = element_text(size = 18),
          plot.margin = unit(c(0.3, 0, 0, 0), "cm"))
  
}


spline_plot <- ggplot(data=contrast_spline_df, aes(x=qrisk2_5yr_score, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=qrisk2_5yr_score, y=exp(Contrast)), size=1) +
  xlab(expression(paste("5-year absolute cardiovascular disease risk (%)"))) +
  ylab("CVD HR for SGLT2i vs comparator") +
  ylim(0.5, 2) +
  coord_trans(y = "log10") +
  scale_x_continuous(breaks = seq(0,30,5)) +
  geom_ribbon(data=contrast_spline_df, aes(x=qrisk2_5yr_score, ymin=exp(Lower), ymax=exp(Upper)), alpha=0.5) +
  
  geom_hline(yintercept = 1, linetype = "dashed")  +
  
  geom_hline(yintercept = 0.94, linetype = "dashed", size=0.7, color="red")  +
  geom_hline(yintercept = 0.83, linetype = "dashed", size=0.7, color="red")  +
  geom_hline(yintercept = 1.07, linetype = "dashed", size=0.7, color="red")  +
  
  geom_segment(aes(x = 3, xend = 4.8, y = 1.6, yend = 1.6), linetype = "dashed", linewidth=0.7, color="red") +
  
  theme_bw() +
  theme(text = element_text(size = 18),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0.5, 0, 0, 0), "cm"),
        legend.position = "none")


spline_plot <- spline_plot +
  annotate(geom="text", x=5, y=2, label="Overall HR in study cohort (95% CI): 0.84 (0.77-0.91)", color="black", hjust=0, size = 6) +
  annotate(geom="text", x=5, y=1.8, label="p-value for CVD risk * drug arm interaction on CVD outcome: 0.12", color="black", hjust=0, size = 6) +
  annotate(geom="text", x=5, y=1.6, label="Trial meta-analysis HR (95% CI): 0.94 (0.83-1.07)", color="red", hjust=0, size = 6) 




hist.dta <- cohort %>% filter(qrisk2_5yr_score>=c1 &  qrisk2_5yr_score <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "qrisk2_5yr_score")



setwd("/slade/CPRD_data/Katie SGLT2/Plots/")
tiff("HR_by_CVD_risk.tiff", width=9, height=7, units = "in", res=800) 

plot_grid(spline_plot, x_hist, ncol = 1,align = 'v',
          rel_heights = c(1,0.5), rel_widths = c(1,1))

dev.off()
