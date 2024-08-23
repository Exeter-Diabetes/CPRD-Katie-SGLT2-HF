
# Uncalibrated and calibrated QRISK2

# HR plot by CVD risk (calibrated QRISK2)


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

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# 1 Cohort selection and variable setup

setwd("/slade/CPRD_data/Katie SGLT2/Processed data/")
load("treatment_outcome_cohort_jun24.rda")
#169,041

table(cohort$studydrug)
# DPP4SU 111673   
# SGLT2 57368  


## B Make variables for survival analysis of all endpoints (see survival_variables function for details)

setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions/")
source("survival_variables.R")

cohort <- add_surv_vars(cohort, main_only=FALSE)


## Make sure don't accidentally use HF variables
cohort <- cohort %>% select(-c(hf_censtime_yrs, narrow_hf_censtime_yrs, qdiabeteshf_5yr_score))
         

############################################################################################

# Look at angina and TIA in addition

cohort <- cohort %>%
  
  mutate(full_cvd_outcome=pmin(postdrug_first_myocardialinfarction,
                               postdrug_first_angina,
                               postdrug_first_stroke,
                               postdrug_first_tia,
                               cv_death_date_any_cause,
                               na.rm=TRUE),
         
         full_cvd_censdate=pmin(full_cvd_outcome, cens_itt, na.rm=TRUE),
         full_cvd_censvar=ifelse(!is.na(full_cvd_outcome) & full_cvd_censdate==full_cvd_outcome, 1, 0),
         full_cvd_censtime_yrs=as.numeric(difftime(full_cvd_censdate, dstartdate, unit="days"))/365.25)
         
prop.table(table(cohort$mace_censvar))        
#2.73% with outcome

prop.table(table(cohort$full_cvd_censvar)) 
#3.96%



############################################################################################

# Calibrate QRISK2 - in whole cohort

dpp4su <- cohort %>%
  filter(studydrug=="DPP4SU") %>%
  mutate(qrisk2_decile=ntile(qrisk2_5yr_score, 10))

cal_females <- dpp4su %>% filter(malesex==0)
cal_males <- dpp4su %>% filter(malesex==1)

female_recal_mod <- coxph(as.formula("Surv(mace_censtime_yrs, mace_censvar)~offset(qrisk2_lin_predictor)"), data=cal_females)
female_qrisk2_surv <- summary(survfit(female_recal_mod),time=5)$surv

male_recal_mod <- coxph(as.formula("Surv(mace_censtime_yrs, mace_censvar)~offset(qrisk2_lin_predictor)"), data=cal_males)
male_qrisk2_surv <- summary(survfit(male_recal_mod),time=5)$surv

cohort <- cohort %>%
  group_by(malesex) %>%
  mutate(centred_qrisk2_lin_predictor=qrisk2_lin_predictor-mean(qrisk2_lin_predictor)) %>%
  ungroup() %>%
  mutate(qrisk2_5yr_score_cal=ifelse(malesex==1, (1-(male_qrisk2_surv^exp(centred_qrisk2_lin_predictor)))*100, (1-(female_qrisk2_surv^exp(centred_qrisk2_lin_predictor)))*100))

print(paste0("Mean linear predictor male: ", dpp4su %>% filter(malesex==1) %>% summarise(mean_qrisk2_lin_predictor=mean(qrisk2_lin_predictor))))
print(paste0("Mean linear predictor female: ", dpp4su %>% filter(malesex==0) %>% summarise(mean_qrisk2_lin_predictor=mean(qrisk2_lin_predictor))))
print(paste0("New male surv: ", male_qrisk2_surv))
print(paste0("New female surv: ", female_qrisk2_surv))

# [1] "Mean linear predictor male: 2.28977762600287"
# [1] "Mean linear predictor female: 2.62290329798239"
# [1] "New male surv: 0.940433624930564"
# [1] "New female surv: 0.956672342740104"


# Plot
dpp4su <- cohort %>%
  filter(studydrug=="DPP4SU") %>%
  mutate(qrisk2_decile=ntile(qrisk2_5yr_score, 10))

predicted <- dpp4su %>%
  group_by(qrisk2_decile) %>%
  summarise(median_qrisk2_pred=median(qrisk2_5yr_score)/100,
            median_qrisk2_pred_cal=median(qrisk2_5yr_score_cal)/100)

observed_mace <- survfit(Surv(mace_censtime_yrs, mace_censvar) ~ qrisk2_decile, data=dpp4su) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed=1-estimate,
         lower_ci=1-conf.low,
         upper_ci=1-conf.high)

observed_maceplus <- survfit(Surv(full_cvd_censtime_yrs, full_cvd_censvar) ~ qrisk2_decile, data=dpp4su) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_maceplus=1-estimate,
         lower_ci_maceplus=1-conf.low,
         upper_ci_maceplus=1-conf.high)


events_table <- t(dpp4su %>%
                    filter(mace_censvar==1) %>%
                    group_by(qrisk2_decile) %>%
                    summarise(events=n()))


obs_v_pred <- cbind(predicted, observed_mace, observed_maceplus) %>%
  select(qrisk2_decile, median_qrisk2_pred, median_qrisk2_pred_cal, observed, lower_ci, upper_ci, observed_maceplus, lower_ci_maceplus, upper_ci_maceplus)

dodge <- position_dodge(width=0.3)

tiff("../../Plots/qrisk2_with_full_cvd.tiff", width=6.5, height=5.5, units = "in", res=800)

ggplot(data=obs_v_pred, aes(x=qrisk2_decile)) +
  geom_point(aes(y = observed*100, color="observed", shape="observed"), size=3, position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100, ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = observed_maceplus*100, color="observed2", shape="observed2"), size=3, position=dodge) +
  geom_errorbar(aes(ymax=upper_ci_maceplus*100, ymin=lower_ci_maceplus*100),width=0.25,size=1, position=dodge, colour="grey") +
  geom_point(aes(y = median_qrisk2_pred*100, color="median_qrisk2_pred", shape="median_qrisk2_pred"), position=dodge, stroke=1.5, size=4) +
  geom_point(aes(y = median_qrisk2_pred_cal*100, color="median_qrisk2_pred_cal", shape="median_qrisk2_pred_cal"), position=dodge, stroke=1.5, size=4) +
  theme_bw() +
  xlab("Predicted 5-year absolute CVD risk decile") + ylab("5-year cardiovascular disease risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,40,by=5)), limits=c(0,41)) +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"),
        axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold"),
        plot.margin = margin(10,10,10,10),
        axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        legend.position=c(0.5, 0.9),
        legend.text=element_text(size=14)) +
  scale_color_manual(values = c(observed2="grey", observed = "black", median_qrisk2_pred = "blue", median_qrisk2_pred_cal = "red"), labels = c(observed2 = "Observed 3P-MACE + angina + TIA (with 95% CI)", observed = "Observed 3P-MACE (with 95% CI)", median_qrisk2_pred = "Median predicted risk per decile before calibration", median_qrisk2_pred_cal = "Median predicted risk per decile after calibration\nto 3P-MACE outcome"), name="") +
  scale_shape_manual(values = c(observed2=16, observed = 16, median_qrisk2_pred = 4, median_qrisk2_pred_cal = 4), labels = c(observed2 = "Observed 3P-MACE + angina + TIA (with 95% CI)", observed = "Observed 3P-MACE (with 95% CI)", median_qrisk2_pred = "Median predicted risk per decile before calibration", median_qrisk2_pred_cal = "Median predicted risk per decile after calibration\nto 3P-MACE outcome"), name="") #+
#ggtitle("5-year QDiabetes-Heart Failure vs heart failure incidence")

dev.off()
# 
# 
# 
# ## C-score
# 
# ### Before
# dpp4su <- dpp4su %>% mutate(qrisk2_survival=(100-qrisk2_5yr_score)/100)
# 
# surv_mod <- coxph(Surv(mace_censtime_yrs, mace_censvar)~qrisk2_survival, data=dpp4su, method="breslow")
# cstat <- round(summary(surv_mod)$concordance[1], 3)
# cstat_lower <- round(summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2]), 3)
# cstat_upper <- round(summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2]), 3)
# #0.68 (0.67, 0.69)
# 
# 
# ### After
# dpp4su <- dpp4su %>% mutate(qrisk2_survival=(100-qrisk2_5yr_score_cal)/100)
# 
# surv_mod <- coxph(Surv(mace_censtime_yrs, mace_censvar)~qrisk2_survival, data=dpp4su, method="breslow")
# cstat <- round(summary(surv_mod)$concordance[1], 3)
# cstat_lower <- round(summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2]), 3)
# cstat_upper <- round(summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2]), 3)
# #0.68 (0.67, 0.69)
# 
# 
# 
# ## Brier score
# 
# library(riskRegression)
# 
# ### Before
# dpp4su <- dpp4su %>% mutate(qrisk2_survival=(100-qrisk2_5yr_score)/100)
# raw_mod <- coxph(Surv(mace_censtime_yrs, mace_censvar) ~ qrisk2_survival, 
#                  data=dpp4su, x=T)
# 
# 
# ## At time 5 years
# score_raw <- 
#   Score(object = list(raw_mod),
#         formula = Surv(mace_censtime_yrs, mace_censvar) ~ 1, # null model
#         data = dpp4su,
#         metrics="brier",
#         times=5) 
# 
# score_raw$Brier$score$Brier[2]  #[1] is NULL model
# #0.058
# score_raw$Brier$score$lower[2]
# #0.055
# score_raw$Brier$score$upper[2]
# #0.060
# 
# 
# ### After
# dpp4su <- dpp4su %>% mutate(qrisk2_survival=(100-qrisk2_5yr_score_cal)/100)
# raw_mod <- coxph(Surv(mace_censtime_yrs, mace_censvar) ~ qrisk2_survival, 
#                  data=dpp4su, x=T)
# 
# 
# ## At time 5 years
# score_raw <- 
#   Score(object = list(raw_mod),
#         formula = Surv(mace_censtime_yrs, mace_censvar) ~ 1, # null model
#         data = dpp4su,
#         metrics="brier",
#         times=5) 
# 
# score_raw$Brier$score$Brier[2]  #[1] is NULL model
# #0.058
# score_raw$Brier$score$lower[2]
# #0.055
# score_raw$Brier$score$upper[2]
# #0.060
# 


## Both look unchanged?? (Brier, very slightly)


############################################################################################

# 2 Calculate hazard ratios -with all sensitivity analyses

setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions")
source("full_covariate_set.R")
# now using functions to define covariates adjusted and weighted for

return_cov_set("weight")
## replace qdhf with calibrated qrisk2
weight <- "malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + INS + initiation_year + qrisk2_5yr_score_cal"

return_cov_set("adjust")
## replace qdhf with calibrated qrisk2
adjust <- "malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year + qrisk2_5yr_score_cal"



## Overlap and IPTW weighting

ps.formula <- formula(paste0("studydrug ~ ", weight))

overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(cohort), weight = c("IPW", "overlap"))

cohort$overlap_weights <- overlap$ps.weights$overlap
cohort$iptw_weights <- overlap$ps.weights$IPW

f_adjusted <- as.formula(paste0("Surv(mace_censtime_yrs, mace_censvar) ~ studydrug  +  ", adjust))

f_unadjusted <- as.formula(paste0("Surv(mace_censtime_yrs, mace_censvar) ~ studydrug"))


## Different weighting

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



## MACE + angina + TIA

f_adjusted <- as.formula(paste0("Surv(full_cvd_censtime_yrs, full_cvd_censvar) ~  studydrug  +  ", adjust))

counts <- cohort %>%
  group_by(studydrug) %>%
  summarise(total_count=n(),
            event_count=sum(full_cvd_censvar),
            total_time=sum(full_cvd_censtime_yrs),
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
  mutate(analysis="broad mace")
#Gives log likelihood error but coefficients seem fine - none are infinite

all_counts <- all_counts %>% bind_rows(counts)

all_hrs <- all_hrs %>% bind_rows(adjoverlap)



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



## Subset on insulin

on_ins <- cohort %>% filter(INS==1)

return_cov_set("weight")
# remove insulin

ps.formula <- formula(paste("studydrug ~ malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + initiation_year + qrisk2_5yr_score_cal"))

return_cov_set("adjust")
# remove insulin

f_adjusted <- as.formula(Surv(mace_censtime_yrs, mace_censvar) ~  studydrug  +  malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + initiation_year + qrisk2_5yr_score_cal)


overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(on_ins), weight = c("IPW", "overlap"))

on_ins$overlap_weights <- overlap$ps.weights$overlap

counts <- on_ins %>%
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

setwd("/slade/CPRD_data/Katie SGLT2/Processed data/")

save(hrs_to_plot, file="hr_data_cvd.Rda")

load("hr_data_cvd.Rda")


# Format for plot

text <- hrs_to_plot %>%
  mutate(hr_text=paste0(round_pad(estimate, 2), " (", round_pad(conf.low, 2), "-", round_pad(conf.high, 2), ")")) %>%
  select(analysis, DPP4SU_count, DPP4SU_ir, SGLT2_count, SGLT2_ir, mean=estimate, lower=conf.low, upper=conf.high, hr_text)


tiff("/slade/CPRD_data/Katie SGLT2/Plots/HR_methods_CVD.tiff", width=18, height=9.2, units = "in", res=400) 

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
    gpar(fill = "red", col="red")
  )
)


text %>% forestplot(labeltext = list(analysis_text=list("Primary analysis: adjusted\nwith overlap weighting", "Unadjusted with overlap weighting", "Adjusted with no weighting", "Adjusted with IPTW", expression("3P-MACE + angina + TIA as outcome"^1), expression("More specific 3P-MACE definition"^2), expression("Patients on either DPP4i or SU"^3), expression("Death as a competing risk"^4), "Insulin-treated subgroup only", "Trial meta-analysis"), SGLT2_count=text$SGLT2_count, SGLT2_ir=text$SGLT2_ir, DPP4SU_count=text$DPP4SU_count, DPP4SU_ir=text$DPP4SU_ir, hr_text=text$hr_text),
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
                    col = fpColors(zero = c("grey50"), text=c('black', 'black', 'black', 'black','black','black','black','black','black','black', 'red')),
                    lwd.ci=1,
                    ci.vertices.height = 0.2,
                    colgap = unit(6, "mm"),
                    txt_gp = fpTxtGp(xlab=gpar(cex=1.8, fontface="bold"), ticks=gpar(cex=1.7), label=gpar(cex=1.7))) %>%
  fp_add_header(analysis_text = "", SGLT2_count = "SGLT2i\nEvents", SGLT2_ir = "\nIR", DPP4SU_count = "Control\nEvents", DPP4SU_ir = "\nIR", hr_text = "Hazard ratio\n(95% CI)")


dev.off()



############################################################################################


# HR by CVD risk

## Keep variables of interest only so don't accidentally use HF ones
cohort <- cohort %>%
  select(patid, studydrug, malesex, dstartdate_age, dstartdate_dm_dur_all, ethnicity_qrisk2_decoded, imd2015_10, qrisk2_smoking_cat, hypertension, predrug_af, hosp_admission_prev_year_count, prebmi, prehba1c2yrs, presbp, drugline_all, ncurrtx_cat, INS, initiation_year, qrisk2_5yr_score_cal, mace_censtime_yrs, mace_censvar)


# Overall HR:

setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions")
source("full_covariate_set.R")

return_cov_set("weight")
"malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score"
# replace qdiabeteshf_5yr_score with qrisk2_5yr_score_cal

ps.formula <- formula("studydrug ~ malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + INS + initiation_year + qrisk2_5yr_score_cal")

overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(cohort), weight = c("IPW", "overlap"))

cohort$overlap_weights <- overlap$ps.weights$overlap

return_cov_set("adjust")
"malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score"

f_adjusted <- as.formula("Surv(mace_censtime_yrs, mace_censvar) ~ studydrug  +  malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year + qrisk2_5yr_score_cal")

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
#remove qdiabeteshf_5yr_score and add qrisk2 (cal) as interaction term


model <- cph(Surv(mace_censtime_yrs, mace_censvar) ~  studydrug*rcs(qrisk2_5yr_score_cal,5) + malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year, data=cohort, x=T, y=T)

anova(model)
# no interaction - good (p=0.13)

c1 <- quantile(cohort$qrisk2_5yr_score_cal, .01, na.rm=TRUE)
c99 <- quantile(cohort$qrisk2_5yr_score_cal, .99, na.rm=TRUE)


contrast_spline <- contrast(model, list(studydrug = "SGLT2", qrisk2_5yr_score_cal = seq(c1, c99, by=0.05)), list(studydrug = "DPP4SU", qrisk2_5yr_score_cal = seq(c1, c99, by=0.05)))

contrast_spline_df <- as.data.frame(contrast_spline[c('qrisk2_5yr_score_cal','Contrast','Lower','Upper')])



# Add histogram
marginal_distribution <- function(x,var) {
  ggplot(x, aes_string(x = var)) +
    geom_histogram(bins = 64, alpha = 0.4, position = "identity") +
    guides(fill = FALSE) +
    theme(legend.title = element_blank(), panel.background = element_rect( fill = "white",color = "grey50")) +
    scale_x_continuous(breaks = seq(0,25,2.5)) +
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


spline_plot <- ggplot(data=contrast_spline_df, aes(x=qrisk2_5yr_score_cal, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=qrisk2_5yr_score_cal, y=exp(Contrast)), size=1) +
  xlab(expression(paste("5-year absolute cardiovascular disease risk (%)"))) +
  ylab("Cardiovascular disease HR for SGLT2i vs control") +
  ylim(0.5, 2.5) +
  coord_trans(y = "log10") +
  scale_x_continuous(breaks = seq(0,25,2.5)) +
  geom_ribbon(data=contrast_spline_df, aes(x=qrisk2_5yr_score_cal, ymin=exp(Lower), ymax=exp(Upper)), alpha=0.5) +
  
  geom_hline(yintercept = 1, linetype = "dashed")  +
  
  geom_hline(yintercept = 0.94, linetype = "dashed", size=0.7, color="red")  +
  geom_hline(yintercept = 0.83, linetype = "dashed", size=0.7, color="red")  +
  geom_hline(yintercept = 1.07, linetype = "dashed", size=0.7, color="red")  +
  
  geom_segment(aes(x = 1.9, xend = 2.8, y = 1.9, yend = 1.9), linetype = "dashed", linewidth=0.7, color="red") +
  
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
  annotate(geom="text", x=3, y=2.4, label="Overall HR in study cohort (95% CI): 0.84 (0.77-0.91)", color="black", hjust=0, size = 6) +
  annotate(geom="text", x=3, y=2.2, label="p-value for CVD risk * drug arm interaction on CVD outcome: 0.13", color="black", hjust=0, size = 6) +
  annotate(geom="text", x=3, y=1.9, label="Trial meta-analysis HR (95% CI): 0.94 (0.83-1.07)", color="red", hjust=0, size = 6) 




hist.dta <- cohort %>% filter(qrisk2_5yr_score_cal>=c1 &  qrisk2_5yr_score_cal <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "qrisk2_5yr_score_cal")



setwd("/slade/CPRD_data/Katie SGLT2/Processed data/Plots/")
tiff("HR_by_CVD_risk.tiff", width=9, height=8, units = "in", res=800) 

plot_grid(spline_plot, x_hist, ncol = 1,align = 'v',
          rel_heights = c(1,0.4), rel_widths = c(1,1))

dev.off()




