
# HRs for SGLT2 vs DPP4SU for HF

# All sensitivity analyses: different adjustment/weight methods, no DPP4i / no SU, narrow HF, different censoring, death as competing risk, insulin, sex

 
############################################################################################

# Setup
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(PSweight)
library(forestplot)
library(rms)
library(tidycmprsk)
library(table1)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# 1 Cohort selection (see script 00)

setwd("/slade/CPRD_data/Katie SGLT2/Processed data/")
load("treatment_outcome_cohort_jun24.rda")

# Add survival variables
setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions/")
source("survival_variables.R")

cohort <- add_surv_vars(cohort)


## Overall median followup
summary(cohort$hf_censtime_yrs)
#2.0 (0.9-3.6)


############################################################################################

# 2 Calculate hazard ratios

setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions/")
source("full_covariate_set.R")
# now using functions to define covariates adjusted and weighted for

## PS scores and adjustment for main analyses

ps.formula <- formula(paste0("studydrug ~ ", return_cov_set("weight")))

f_adjusted <- as.formula(paste0("Surv(hf_censtime_yrs, hf_censvar) ~ studydrug  +  ", return_cov_set("adjust")))

f_unadjusted <- as.formula(paste0("Surv(hf_censtime_yrs, hf_censvar) ~ studydrug"))



## Main analysis

overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(cohort), weight = c("IPW", "overlap"))

cohort$overlap_weights <- overlap$ps.weights$overlap
cohort$iptw_weights <- overlap$ps.weights$IPW

counts <- cohort %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(hf_censvar),
            total_time=sum(hf_censtime_yrs),
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
            event_count=sum(hf_censvar),
            total_time=sum(hf_censtime_yrs),
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
            event_count=sum(hf_censvar),
            total_time=sum(hf_censtime_yrs),
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

unadjoverlap <- coxph(f_unadjusted, cohort, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="unadjusted overlap weight")



adjusted <- coxph(f_adjusted, cohort) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="no weight")



#violated for predrug_af, hosp_admission_prev_year_count, prebmi, presbp, drugline_all and initiation_year - mostly for SGLT2i group

myfit <- survfit(Surv(hf_censtime_yrs,hf_censvar) ~ studydrug + predrug_af,data=cohort)
ggsurvplot(myfit, conf.int=TRUE)

myfit <- survfit(Surv(hf_censtime_yrs,hf_censvar) ~ studydrug + predrug_af, weights=cohort$overlap_weights, data=cohort)
ggsurvplot(myfit, conf.int=TRUE)


cohort <- cohort %>% mutate(predrug_af=as.factor(predrug_af))
myfit <- survfit(Surv(hf_censtime_yrs,hf_censvar) ~ studydrug + predrug_af,data=cohort)
cox_unadj <- coxph(Surv(hf_censtime_yrs,hf_censvar) ~ studydrug + predrug_af,data=cohort)

dummy <- list(predrug_af = levels(cohort$predrug_af), studydrug=levels(cohort$studydrug))
cox_surv <- survfit(cox_unadj, newdata = dummy)

plot(myfit, col = c(2, 4), lty = 1,
     xlab = 'years', ylab = 'survival')
lines(cox_surv, col = c(2, 4), lty = 2)
legend('topright',
       c('observed', 'expected', levels(undiag$hba1c_cat)),
       lty = c(1, 2, 4, 4), col = c(1, 1, 2, 4))
dev.off()




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



## Narrow HF outcome

f_adjusted_narrow_hf <- as.formula(paste0("Surv(narrow_hf_censtime_yrs, narrow_hf_censvar) ~  studydrug  +  ", return_cov_set("adjust")))

counts <- cohort %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(narrow_hf_censvar),
            total_time=sum(narrow_hf_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

adjoverlap <- coxph(f_adjusted_narrow_hf, cohort, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="narrow HF")
#Gives log likelihood error but coefficients seem fine - none are infinite



all_counts <- all_counts %>% bind_rows(counts)

all_hrs <- all_hrs %>% bind_rows(adjoverlap)



## Second line therapy after MFN only
### A) Intention to treat
### B) Per-protocol

second_line_cohort <- cohort %>% filter(drugline_all==2 & multi_drug_start==0 & MFN==1) #NB: drugline same as drugorder in this dataset as only have 1st instance of each drug class
#64,304

## need to remove drugline and ncurrtx from propensity score and adjustment

ps.formula.trial <- formula(paste("studydrug ~ malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + initiation_year + qdiabeteshf_5yr_score"))

overlap <- SumStat(ps.formula=ps.formula.trial, data = as.data.frame(second_line_cohort), weight = "overlap")

second_line_cohort$overlap_weights <- overlap$ps.weights$overlap

## A) ITT

f_adjusted_trial_itt <- as.formula(Surv(hf_itt_censtime_yrs, hf_itt_censvar) ~  studydrug  +  malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + INS + initiation_year + qdiabeteshf_5yr_score)


counts <- second_line_cohort %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(hf_itt_censvar),
            total_time=sum(hf_itt_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

adjoverlap <- coxph(f_adjusted_trial_itt, second_line_cohort, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="trial_itt")



all_counts <- all_counts %>% bind_rows(counts)

all_hrs <- all_hrs %>% bind_rows(adjoverlap)


## B) PP

f_adjusted_trial_pp <- as.formula(Surv(hf_pp_censtime_yrs, hf_pp_censvar) ~  studydrug  +  malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score)


counts <- second_line_cohort %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(hf_pp_censvar),
            total_time=sum(hf_pp_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))


adjoverlap <- coxph(f_adjusted_trial_pp, second_line_cohort, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="trial_pp")



all_counts <- all_counts %>% bind_rows(counts)

all_hrs <- all_hrs %>% bind_rows(adjoverlap)



## Death as competing risk

cohort <- cohort %>% 
  mutate(hf_censvar.cr=as.factor(ifelse(hf_censvar==0 & !is.na(death_date) & death_date==hf_censdate, 2, hf_censvar)))

counts <- cohort %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(hf_censvar.cr==1),
            total_time=sum(hf_censtime_yrs),
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
            event_count=sum(hf_censvar.cr==2),
            total_time=sum(hf_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

# Adjusted with overlap weighting
## Very slow to run

f_adjusted_crr <- as.formula(paste0("Surv(hf_censtime_yrs, hf_censvar.cr) ~  studydrug  +  ", return_cov_set("adjust")))

adjoverlap <- crr(f_adjusted_crr, data = cohort, weights=overlap) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="death competing")



all_counts <- all_counts %>% bind_rows(counts)

all_hrs <- all_hrs %>% bind_rows(adjoverlap)



## Subset on insulin

on_ins <- cohort %>% filter(INS==1)

return_cov_set("weight")
# remove insulin

ps.formula.insulin <- formula(paste("studydrug ~ malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + initiation_year + qdiabeteshf_5yr_score"))

return_cov_set("adjust")
# remove insulin

f_adjusted_insulin <- as.formula(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug  +  malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + initiation_year + qdiabeteshf_5yr_score)


overlap <- SumStat(ps.formula=ps.formula.insulin, data = as.data.frame(on_ins), weight = c("IPW", "overlap"))

on_ins$overlap_weights <- overlap$ps.weights$overlap

counts <- on_ins %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(hf_censvar),
            total_time=sum(hf_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

adjoverlap <- coxph(f_adjusted_insulin, on_ins, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="on insulin")

##LLs converge but looks OK



all_counts <- all_counts %>% bind_rows(counts)

all_hrs <- all_hrs %>% bind_rows(adjoverlap)



## By sex

### Males only

males <- cohort %>% filter(malesex==1)

return_cov_set("weight")
# remove sex

ps.formula.sex <- formula(paste("studydrug ~ dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score"))

return_cov_set("adjust")
# remove sex

f_adjusted_sex <- as.formula(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score)


overlap <- SumStat(ps.formula=ps.formula.sex, data = as.data.frame(males), weight = c("IPW", "overlap"))

males$overlap_weights <- overlap$ps.weights$overlap

counts <- males %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(hf_censvar),
            total_time=sum(hf_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

adjoverlap <- coxph(f_adjusted_sex, males, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="males")



all_counts <- all_counts %>% bind_rows(counts)

all_hrs <- all_hrs %>% bind_rows(adjoverlap)


### Females only

females <- cohort %>% filter(malesex==0)

overlap <- SumStat(ps.formula=ps.formula.sex, data = as.data.frame(females), weight = c("IPW", "overlap"))

females$overlap_weights <- overlap$ps.weights$overlap

counts <- females %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(hf_censvar),
            total_time=sum(hf_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

adjoverlap <- coxph(f_adjusted_sex, females, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="females")



all_counts <- all_counts %>% bind_rows(counts)

all_hrs <- all_hrs %>% bind_rows(adjoverlap)



# All methods

hrs_to_plot <- cbind(all_counts, all_hrs) %>%
  rbind(data.frame(DPP4SU_count=NA, DPP4SU_ir=NA, SGLT2_count=NA, SGLT2_ir=NA, analysis="trial", estimate=0.63, conf.low=0.50, conf.high=0.80))

setwd("/slade/CPRD_data/Katie SGLT2/Processed data/")

save(hrs_to_plot, file="hr_data.Rda")

load("hr_data.Rda")


# Format for plot

text <- hrs_to_plot %>%
  mutate(hr_text=paste0(round_pad(estimate, 2), " (", round_pad(conf.low, 2), "-", round_pad(conf.high, 2), ")")) %>%
  select(analysis, DPP4SU_count, DPP4SU_ir, SGLT2_count, SGLT2_ir, mean=estimate, lower=conf.low, upper=conf.high, hr_text)


tiff("/slade/CPRD_data/Katie SGLT2/Plots/HR_methods.tiff", width=18, height=15, units = "in", res=600) 

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
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "red", col="red")
  )
)


text %>% forestplot(labeltext = list(analysis_text=list("Primary analysis: adjusted\nwith overlap weighting", "SU-only comparator group", "DPP4i-only comparator group", "Unadjusted with overlap weighting", "Adjusted with no weighting", "Adjusted with IPTW", expression("More specific HF definition"^1), expression("Second-line intention-to-treat"^2), expression("Second-line per-protocol"^2), expression("Death as a competing risk"^3), "Insulin-treated subgroup only", "Males only", "Females only", "Trial meta-analysis"), SGLT2_count=text$SGLT2_count, SGLT2_ir=text$SGLT2_ir, DPP4SU_count=text$DPP4SU_count, DPP4SU_ir=text$DPP4SU_ir, hr_text=text$hr_text),
           xlog = T,
           ci.vertices = TRUE,
           xlab = "Hazard ratio ",
           fn.ci_norm = fpDrawCircleCI,
           xticks = c(log(0.4), log(0.5), log(0.6), log(0.7), log(0.8), log(0.9), log(1.0)),
           boxsize = 0.25,
           align = c("l", "l", "l", "l", "l"),
           graph.pos=6,
           lwd.zero = 1.2,
           shapes_gp = styles,
           col = fpColors(zero = c("grey50"), text=c('black','black','black', 'black','black','black','black','black','black','black','black','black','black', 'black', 'red')),
           lwd.ci=1,
           ci.vertices.height = 0.2,
           colgap = unit(6, "mm"),
           txt_gp = fpTxtGp(xlab=gpar(cex=1.8, fontface="bold"), ticks=gpar(cex=1.7), label=gpar(cex=1.7))) %>%
  fp_add_header(analysis_text = "", SGLT2_count = "SGLT2i\nEvents", SGLT2_ir = "\nIR", DPP4SU_count = "Comparator\nEvents", DPP4SU_ir = "\nIR", hr_text = "Hazard ratio\n(95% CI)")


dev.off()


############################################################################################

## Also check by ethnicity

table(cohort$ethnicity_5cat, useNA="always")

# 0       1       2      3      4      <NA> 
# 126191  24395   9865   2833   1843   3914 

white <- cohort %>% filter(ethnicity_5cat==0)
overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(white), weight = "overlap")
white$overlap_weights <- overlap$ps.weights$overlap
coxph(f_adjusted, white, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high)

#White
#0.716    0.636     0.806
  
sasian <- cohort %>% filter(ethnicity_5cat==1)
overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(sasian), weight = "overlap")
sasian$overlap_weights <- overlap$ps.weights$overlap
coxph(f_adjusted, sasian, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) 

#South Asian
#0.633    0.442     0.905


black <- cohort %>% filter(ethnicity_5cat==2)
overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(black), weight = "overlap")
black$overlap_weights <- overlap$ps.weights$overlap
coxph(f_adjusted, black, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) 

#Black
#0.618    0.345      1.11
  

other <- cohort %>% filter(ethnicity_5cat==3)
overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(other), weight = "overlap")
other$overlap_weights <- overlap$ps.weights$overlap
coxph(f_adjusted, other, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) 

#Other
# 0.344    0.139     0.849


mixed <- cohort %>% filter(ethnicity_5cat==4)

#have to remove ethnicity from PS score and adjusting as is constant in this group (only)

ps.formula.ethnicity <- formula("studydrug ~ malesex + dstartdate_age + dstartdate_dm_dur_all + 
    imd2015_10 + qrisk2_smoking_cat + 
    hypertension + predrug_af + hosp_admission_prev_year_count + 
    prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + 
    INS + initiation_year + qdiabeteshf_5yr_score")

f_adjusted_ethnicity <- as.formula("Surv(hf_censtime_yrs, hf_censvar) ~ studydrug + malesex + rcs(dstartdate_age, 5) + rcs(dstartdate_dm_dur_all, 5) + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi, 5) + rcs(prehba1c2yrs, 5) + rcs(presbp, 5) + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score")

overlap <- SumStat(ps.formula=ps.formula.ethnicity, data = as.data.frame(mixed), weight = "overlap")
mixed$overlap_weights <- overlap$ps.weights$overlap
coxph(f_adjusted_ethnicity, mixed, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) 

#Mixed
#0.0496   0.0240     0.103



# Osteoporosis

## Try without weights first

f_adjusted <- as.formula(paste0("Surv(hf_censtime_yrs, hf_censvar) ~ studydrug  +  ", return_cov_set("adjust")))

model <- cph(f_adjusted, cohort)

plot(anova(model))

vars <- return_cov_set("adjust")

vars <- paste0(vars, " + predrug_osteoporosis")


f_adjusted <- as.formula(paste0("Surv(hf_censtime_yrs, hf_censvar) ~ studydrug  +  ", vars))

model <- cph(f_adjusted, cohort)

plot(anova(model))
# not that significant?

model <- coxph(f_adjusted, cohort)
# is significantly associated - similar to insulin:  HR: 1.44 (1.17, 1.77)
#Not as much as AF though = 3.04 (2.46, 3.75




## By ethnicity

cohort <- cohort %>%
  mutate(ethnicity_decoded=factor(case_when(ethnicity_qrisk2_decoded=="missing" ~"Missing",
                                            ethnicity_qrisk2_decoded=="White" ~"White",
                                            ethnicity_qrisk2_decoded=="Indian" | ethnicity_qrisk2_decoded=="Pakistani" |  ethnicity_qrisk2_decoded=="Bangladeshi" | ethnicity_qrisk2_decoded=="Other Asian" ~"Asian",
                                            ethnicity_qrisk2_decoded=="Black Caribbean" | ethnicity_qrisk2_decoded=="Black African" ~"Black",
                                            ethnicity_qrisk2_decoded=="Chinese" | ethnicity_qrisk2_decoded=="Other" ~"Other"), levels=c("White", "Asian", "Black", "Other", "Missing")))


white <- cohort %>% filter(ethnicity_decoded=="White")

return_cov_set("weight")
# remove ethnicity

ps.formula.ethnicity <- formula(paste("studydrug ~ dstartdate_age + dstartdate_dm_dur_all + malesex + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score"))

return_cov_set("adjust")
# remove sex

f_adjusted_ethnicity <- as.formula(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + malesex + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score)

overlap <- SumStat(ps.formula=ps.formula.ethnicity, data = as.data.frame(white), weight = c("IPW", "overlap"))

white$overlap_weights <- overlap$ps.weights$overlap

all_eth_count <- white %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(hf_censvar),
            total_time=sum(hf_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

all_eth_adjoverlap <- coxph(f_adjusted_ethnicity, white, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="white")


asian <- cohort %>% filter(ethnicity_decoded=="Asian")

overlap <- SumStat(ps.formula=ps.formula.ethnicity, data = as.data.frame(asian), weight = c("IPW", "overlap"))

asian$overlap_weights <- overlap$ps.weights$overlap

eth_count <- asian %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(hf_censvar),
            total_time=sum(hf_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

eth_adjoverlap <- coxph(f_adjusted_ethnicity, asian, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="asian")


all_eth_count <- all_eth_count %>% rbind(eth_count)
all_eth_adjoverlap <- all_eth_adjoverlap %>% rbind(eth_adjoverlap)


black <- cohort %>% filter(ethnicity_decoded=="Black")

overlap <- SumStat(ps.formula=ps.formula.ethnicity, data = as.data.frame(black), weight = c("IPW", "overlap"))

black$overlap_weights <- overlap$ps.weights$overlap

eth_count <- black %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(hf_censvar),
            total_time=sum(hf_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

eth_adjoverlap <- coxph(f_adjusted_ethnicity, black, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="black")

all_eth_count <- all_eth_count %>% rbind(eth_count)
all_eth_adjoverlap <- all_eth_adjoverlap %>% rbind(eth_adjoverlap)


other <- cohort %>% filter(ethnicity_decoded=="Other")

overlap <- SumStat(ps.formula=ps.formula.ethnicity, data = as.data.frame(other), weight = c("IPW", "overlap"))

other$overlap_weights <- overlap$ps.weights$overlap

eth_count <- other %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(hf_censvar),
            total_time=sum(hf_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

eth_adjoverlap <- coxph(f_adjusted_ethnicity, other, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="other")

all_eth_count <- all_eth_count %>% rbind(eth_count)
all_eth_adjoverlap <- all_eth_adjoverlap %>% rbind(eth_adjoverlap)


missing <- cohort %>% filter(ethnicity_decoded=="Missing")

overlap <- SumStat(ps.formula=ps.formula.ethnicity, data = as.data.frame(missing), weight = c("IPW", "overlap"))

missing$overlap_weights <- overlap$ps.weights$overlap

eth_count <- missing %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(hf_censvar),
            total_time=sum(hf_censtime_yrs),
            count=paste0(event_count, "/", total_count),
            ir=round_pad((event_count*1000)/total_time, 2)) %>%
  select(-c(total_count, event_count, total_time)) %>%
  pivot_wider(names_from=studydrug,
              names_glue="{studydrug}_{.value}",
              values_from=c(count, ir))

eth_adjoverlap <- coxph(f_adjusted_ethnicity, missing, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="missing")

all_eth_count <- all_eth_count %>% rbind(eth_count)
all_eth_adjoverlap <- all_eth_adjoverlap %>% rbind(eth_adjoverlap)




