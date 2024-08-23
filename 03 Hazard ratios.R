
# HRs for SGLT2 vs DPP4SU for HF

# All sensitivity analyses: narrow HF, DPP4 or SU, death as competing risk, no weighting, IPTW

# Add insulin subgroup

# And with no DPP4i / no SU
 
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

setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions")
source("full_covariate_set.R")
# now using functions to define covariates adjusted and weighted for


## Overlap and IPTW weighting

ps.formula <- formula(paste0("studydrug ~ ", return_cov_set("weight")))

overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(cohort), weight = c("IPW", "overlap"))

cohort$overlap_weights <- overlap$ps.weights$overlap
cohort$iptw_weights <- overlap$ps.weights$IPW

f_adjusted <- as.formula(paste0("Surv(hf_censtime_yrs, hf_censvar) ~ studydrug  +  ", return_cov_set("adjust")))

f_unadjusted <- as.formula(paste0("Surv(hf_censtime_yrs, hf_censvar) ~ studydrug"))


## Different weighting

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
  bind_rows(counts) %>%
  bind_rows(counts) %>%
  bind_rows(counts)

all_hrs <- adjoverlap %>%
  bind_rows(unadjoverlap) %>%
  bind_rows(adjusted) %>%
  bind_rows(adjiptw)
 


## Narrow HF outcome

f_adjusted <- as.formula(paste0("Surv(narrow_hf_censtime_yrs, narrow_hf_censvar) ~  studydrug  +  ", return_cov_set("adjust")))

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

adjoverlap <- coxph(f_adjusted, cohort, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="narrow HF")
#Gives log likelihood error but coefficients seem fine 0 none are infinite

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

f_adjusted <- as.formula(paste0("Surv(hf_can_overlap_censtime_yrs, hf_can_overlap_censvar) ~  studydrug  +  ", return_cov_set("adjust")))

counts <- sensitivity_cohort %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(hf_can_overlap_censvar),
            total_time=sum(hf_can_overlap_censtime_yrs),
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

f_adjusted <- as.formula(paste0("Surv(hf_censtime_yrs, hf_censvar.cr) ~  studydrug  +  ", return_cov_set("adjust")))

adjoverlap <- crr(f_adjusted, data = cohort, weights=overlap) %>%
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

ps.formula <- formula(paste("studydrug ~ malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + initiation_year + qdiabeteshf_5yr_score"))

return_cov_set("adjust")
# remove insulin

f_adjusted <- as.formula(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug  +  malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + initiation_year + qdiabeteshf_5yr_score)


overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(on_ins), weight = c("IPW", "overlap"))

on_ins$overlap_weights <- overlap$ps.weights$overlap

counts <- on_ins %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(hf_can_overlap_censvar),
            total_time=sum(hf_can_overlap_censtime_yrs),
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


tiff("/slade/CPRD_data/Katie SGLT2/Plots/HR_methods.tiff", width=18, height=9, units = "in", res=400) 

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
    gpar(fill = "red", col="red")
  )
)


text %>% forestplot(labeltext = list(analysis_text=list("Primary analysis: adjusted\nwith overlap weighting", "Unadjusted with overlap weighting", "Adjusted with no weighting", "Adjusted with IPTW", expression("More specific HF definition"^1), expression("Patients on either DPP4i or SU"^2), expression("Death as a competing risk"^3), "Insulin-treated subgroup only", "Trial meta-analysis"), SGLT2_count=text$SGLT2_count, SGLT2_ir=text$SGLT2_ir, DPP4SU_count=text$DPP4SU_count, DPP4SU_ir=text$DPP4SU_ir, hr_text=text$hr_text),
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
           col = fpColors(zero = c("grey50"), text=c('black','black','black','black','black','black','black','black', 'black', 'red')),
           lwd.ci=1,
           ci.vertices.height = 0.2,
           colgap = unit(6, "mm"),
           txt_gp = fpTxtGp(xlab=gpar(cex=1.8, fontface="bold"), ticks=gpar(cex=1.7), label=gpar(cex=1.7))) %>%
  fp_add_header(analysis_text = "", SGLT2_count = "SGLT2i\nEvents", SGLT2_ir = "\nIR", DPP4SU_count = "Control\nEvents", DPP4SU_ir = "\nIR", hr_text = "Hazard ratio\n(95% CI)")


dev.off()



##########################################################################################################################

# SU-only and DPP4-only groups


# No DPP4i

no_dpp4i <- cohort %>% filter(drugclass!="DPP4")

ps.formula <- formula(paste("studydrug ~ ", return_cov_set("weight")))

overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(no_dpp4i), weight = "overlap")

no_dpp4i$overlap_weights <- overlap$ps.weights$overlap

f_adjusted <- as.formula(paste0("Surv(hf_censtime_yrs, hf_censvar) ~  studydrug  + ", return_cov_set("adjust")))

counts <- no_dpp4i %>%
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

counts
#DPP4SU_count SGLT2_count DPP4SU_ir SGLT2_ir
#826/42965    706/57368   7.69      5.60    

adjoverlap <- coxph(f_adjusted, no_dpp4i, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="overlap weight")

adjoverlap
#0.691    0.597     0.799



# No SU

no_su <- cohort %>% filter(drugclass!="SU")

ps.formula <- formula(paste("studydrug ~ ", return_cov_set("weights")))

overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(no_su), weight = "overlap")

no_su$overlap_weights <- overlap$ps.weights$overlap

f_adjusted <- as.formula(paste0("Surv(hf_censtime_yrs, hf_censvar) ~  studydrug  + ", return_cov_set("adjust")))

counts <- no_su %>%
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

counts
#DPP4SU_count SGLT2_count DPP4SU_ir SGLT2_ir
#1226/68708   706/57368   7.84      5.60 

adjoverlap <- coxph(f_adjusted, no_su, weights=overlap_weights) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(analysis="overlap weight")

adjoverlap
#0.701 0.623 0.788
