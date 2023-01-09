
# Exploring QRISK2 for predicting MACE incidence

## 1 Uncalibrated 5-year QRISK2 vs MACE incidence in both arms (report C-stat + number of events)

## 2 As per 1, but individual components of MACE

## 3 (Uncalibrated) 5-year QRISK2 vs hazard ratio for MACE (is HR constant by baseline QRISK2? Report ANOVA)

## 4 Calibrate 5-year QRISK2 on 20% sample of control arm (DPP4SU) - report C-stat before and after and plot
## Recalculate 5-year QRISK2 for control and study arm (report C-stat before and after and plot)

## 5 Plot unadjusted predicted benefit vs actual benefit by benefit decile/categories
## Plot adjusted version of above (adjust for...)

### 6 Number needed to treat per benefit decile/category for predicted and actual benefit

# Uses Rda dataset from Slade + cohort_definition_function as starting point

# Not including GLP1 for now


############################################################################################

# Setup
library(tidyverse)
library(gtsummary)
library(survival)
library(survminer)
library(broom)
library(gridExtra)
library(patchwork)
library(rms)
library(cowplot)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# Cohort selection and variable setup

## See cohort_definition_function scripts for details

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Raw data/")
load("20221212_t2d_1stinstance.Rda")
load("20221212_t2d_all_drug_periods.Rda")

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Scripts/")
source("cohort_definition_function.R")

cohort <- define_cohort(t2d_1stinstance, t2d_all_drug_periods)

write.table(count(cohort, studydrug), quote=F, col.names=F, row.names=F, sep=" ")
# DPP4SU 90853
# GLP1 12604
# SGLT2 48279

cohort <- cohort %>%
  filter(studydrug!="GLP1") %>%
  select(patid, malesex, ethnicity_5cat_decoded, imd2015_10, regstartdate, gp_record_end, death_date, drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, prehba1c, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preckdstage, qrisk2_smoking_cat, postdrug_first_myocardialinfarction, postdrug_first_primary_incident_mi, postdrug_first_stroke, postdrug_first_primary_incident_stroke, postdrug_first_heartfailure, postdrug_first_primary_hhf, postdrug_first_all_cause_hosp, next_glp1_start, next_sglt2_start, next_tzd_start, last_sglt2_stop, cv_death_date_primary_cause, cv_death_date_any_cause, hf_death_date_primary_cause, hf_death_date_any_cause, qrisk2_lin_predictor, qrisk2_5yr_score, qrisk2_10yr_score, qdiabeteshf_lin_predictor, qdiabeteshf_5yr_score, contains("statins"))

rm(list=setdiff(ls(), "cohort"))


############################################################################################

# Make outcome and survival analysis variables for main analysis
## Broad MACE (hospitalisation/death any cause or in GP records)

cohort <- cohort %>%
  
  mutate(
    
    # broad_mace = GP (MI/stroke) + expanded HES codelist (MI/stroke; any cause) + CV death (any cause)
    postdrug_broad_mace=pmin(postdrug_first_myocardialinfarction, postdrug_first_stroke, cv_death_date_any_cause, na.rm=TRUE),
    
  )

# Broad MACE

# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Start SGLT2 (not SGLT2 arm), GLP1, or TZD
## End of GP records
## Broad MACE outcome

cohort <- cohort %>%
  mutate(five_years_post_dstart=dstartdate+(365.25*5),
         
         mace_broad_censdate=if_else(studydrug=="SGLT2",
                                     pmin(five_years_post_dstart,
                                          death_date,
                                          next_glp1_start,
                                          next_tzd_start,
                                          gp_record_end,
                                          postdrug_broad_mace, na.rm=TRUE),
                                     if_else(studydrug=="DPP4SU",
                                             pmin(five_years_post_dstart,
                                                  death_date,
                                                  next_sglt2_start,
                                                  next_glp1_start,
                                                  next_tzd_start,
                                                  gp_record_end,
                                                  postdrug_broad_mace, na.rm=TRUE),
                                             as.Date(NA))),
         
         mace_broad_censvar=ifelse(!is.na(postdrug_broad_mace) & mace_broad_censdate==postdrug_broad_mace, 1, 0),
         
         mace_broad_censtime_yrs=as.numeric(difftime(mace_broad_censdate, dstartdate, unit="days"))/365.25)


############################################################################################

# 1 How well does uncalibrated QRISK predict 5-year MACE incidence?

## Use ~deciles of QRISK2, but convert to actual values
quantile(cohort$qrisk2_5yr_score, probs = seq(.1, .9, by = .1))
# 3, 5, 6.5, 8, 9.5, 11.5, 13.5, 16, 20

cohort <- cohort %>%
  mutate(qrisk2_cat=case_when(
    qrisk2_5yr_score<=3 ~ 1,
    qrisk2_5yr_score>3 & qrisk2_5yr_score<=5 ~ 2,
    qrisk2_5yr_score>5 & qrisk2_5yr_score<=6.5 ~ 3,
    qrisk2_5yr_score>6.5 & qrisk2_5yr_score<=8 ~ 4,
    qrisk2_5yr_score>8 & qrisk2_5yr_score<=9.5 ~ 5,
    qrisk2_5yr_score>9.5 & qrisk2_5yr_score<=11.5 ~ 6,
    qrisk2_5yr_score>11.5 & qrisk2_5yr_score<=13.5 ~ 7,
    qrisk2_5yr_score>13.5 & qrisk2_5yr_score<=16 ~ 8,
    qrisk2_5yr_score>16 & qrisk2_5yr_score<=20 ~ 9,
    qrisk2_5yr_score>20 ~ 10
  ))

table(cohort$qrisk2_cat, useNA="always")
# No NAs and roughly equal counts



# Get mean predicted probabilities for each QRISK2 category by studydrug
predicted <- cohort %>%
  group_by(qrisk2_cat, studydrug) %>%
  summarise(mean_pred=mean(qrisk2_5yr_score)/100)


# Find actual observed probabilities by QRISK2 category and studydrug
observed_dpp4su <- survfit(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_dpp4su=1-estimate,
         lower_ci_dpp4su=1-conf.low,
         upper_ci_dpp4su=1-conf.high) %>%
  select(observed_dpp4su, lower_ci_dpp4su, upper_ci_dpp4su, strata)

dpp4su_events <- cohort %>%
  filter(studydrug=="DPP4SU" & mace_broad_censvar==1) %>%
  group_by(qrisk2_cat) %>%
  summarise(DPP4SU=n())

observed_sglt2 <- survfit(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_sglt2=1-estimate,
         lower_ci_sglt2=1-conf.low,
         upper_ci_sglt2=1-conf.high) %>%
  select(observed_sglt2, lower_ci_sglt2, upper_ci_sglt2, strata)

sglt2_events <- cohort %>%
  filter(studydrug=="SGLT2" & mace_broad_censvar==1) %>%
  group_by(qrisk2_cat) %>%
  summarise(SGLT2=n())


obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug=="DPP4SU")), observed_dpp4su),
  cbind((predicted %>% filter(studydrug=="SGLT2")), observed_sglt2)
) %>%
  mutate(observed=coalesce(observed_dpp4su, observed_sglt2),
         lower_ci=coalesce(lower_ci_dpp4su, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_dpp4su, upper_ci_sglt2))

events_table <- data.frame(t(dpp4su_events %>%
  inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="qrisk2_cat")


dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(qrisk2_cat==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_pred=NA, qrisk2_cat=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=qrisk2_cat, group=studydrug, color=studydrug)) + 
  geom_point(aes(y = mean_pred*100), position=dodge) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1,width=0.25,size=1, position=dodge) +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,25,by=5)), limits=c(0,26)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Uncalibrated QRISK2 vs MACE incidence (5 year)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))


# Way underestimates


############################################################################################

# 2 Look at contributions of different MACE components

cohort <- cohort %>%
  mutate(which_mace=ifelse(!is.na(postdrug_broad_mace) & !is.na(postdrug_first_myocardialinfarction) & !is.na(postdrug_first_stroke) & !is.na(cv_death_date_primary_cause) & postdrug_broad_mace==postdrug_first_myocardialinfarction & postdrug_broad_mace==postdrug_first_stroke & postdrug_broad_mace==cv_death_date_primary_cause, "all three",
                           ifelse(!is.na(postdrug_broad_mace) & !is.na(postdrug_first_myocardialinfarction) & !is.na(postdrug_first_stroke) & postdrug_broad_mace==postdrug_first_myocardialinfarction & postdrug_broad_mace==postdrug_first_stroke, "mi & stroke",
                                  ifelse(!is.na(postdrug_broad_mace) & !is.na(postdrug_first_myocardialinfarction) & !is.na(cv_death_date_primary_cause) & postdrug_broad_mace==postdrug_first_myocardialinfarction & postdrug_broad_mace==cv_death_date_primary_cause, "mi & death",
                                         ifelse(!is.na(postdrug_broad_mace) & !is.na(postdrug_first_stroke) & !is.na(cv_death_date_primary_cause) & postdrug_broad_mace==postdrug_first_stroke & postdrug_broad_mace==cv_death_date_primary_cause, "stroke & death",
                                                ifelse(!is.na(postdrug_broad_mace) & !is.na(postdrug_first_myocardialinfarction) & postdrug_broad_mace==postdrug_first_myocardialinfarction, "mi",
                                                       ifelse(!is.na(postdrug_broad_mace) & !is.na(postdrug_first_stroke) & postdrug_broad_mace==postdrug_first_stroke, "stroke",
                                                              ifelse(!is.na(postdrug_broad_mace) & !is.na(cv_death_date_primary_cause) & postdrug_broad_mace==cv_death_date_primary_cause, "cv death", NA))))))))


prop.table(table(cohort$which_mace))
# 13% CV death, 44% MI, 42% stroke

prop.table(table(cohort$studydrug, cohort$which_mace), margin=1)
# 13% CV death, 44% MI, 42% stroke


# Make new survival variables; don't censor at postdrug_broad_mace
cohort <- cohort %>%
  mutate(mi_censdate=if_else(studydrug=="SGLT2",
                                     pmin(five_years_post_dstart,
                                          death_date,
                                          next_glp1_start,
                                          next_tzd_start,
                                          gp_record_end,
                                          postdrug_first_myocardialinfarction, na.rm=TRUE),
                                     if_else(studydrug=="DPP4SU",
                                             pmin(five_years_post_dstart,
                                                  death_date,
                                                  next_sglt2_start,
                                                  next_glp1_start,
                                                  next_tzd_start,
                                                  gp_record_end,
                                                  postdrug_first_myocardialinfarction, na.rm=TRUE),
                                             as.Date(NA))),
         
         mi_censvar=ifelse(!is.na(postdrug_first_myocardialinfarction) & mi_censdate==postdrug_first_myocardialinfarction, 1, 0),
         
         mi_censtime_yrs=as.numeric(difftime(mi_censdate, dstartdate, unit="days"))/365.25,
         
         
         cv_death_censdate=if_else(studydrug=="SGLT2",
                             pmin(five_years_post_dstart,
                                  death_date,
                                  next_glp1_start,
                                  next_tzd_start,
                                  gp_record_end,
                                  cv_death_date_any_cause, na.rm=TRUE),
                             if_else(studydrug=="DPP4SU",
                                     pmin(five_years_post_dstart,
                                          death_date,
                                          next_sglt2_start,
                                          next_glp1_start,
                                          next_tzd_start,
                                          gp_record_end,
                                          cv_death_date_any_cause, na.rm=TRUE),
                                     as.Date(NA))),
         
         cv_death_censvar=ifelse(!is.na(cv_death_date_any_cause) & cv_death_censdate==cv_death_date_any_cause, 1, 0),
         
         cv_death_censtime_yrs=as.numeric(difftime(cv_death_censdate, dstartdate, unit="days"))/365.25,
         
         
         stroke_censdate=if_else(studydrug=="SGLT2",
                             pmin(five_years_post_dstart,
                                  death_date,
                                  next_glp1_start,
                                  next_tzd_start,
                                  gp_record_end,
                                  postdrug_first_stroke, na.rm=TRUE),
                             if_else(studydrug=="DPP4SU",
                                     pmin(five_years_post_dstart,
                                          death_date,
                                          next_sglt2_start,
                                          next_glp1_start,
                                          next_tzd_start,
                                          gp_record_end,
                                          postdrug_first_stroke, na.rm=TRUE),
                                     as.Date(NA))),
         
         stroke_censvar=ifelse(!is.na(postdrug_first_stroke) & stroke_censdate==postdrug_first_stroke, 1, 0),
         
         stroke_censtime_yrs=as.numeric(difftime(stroke_censdate, dstartdate, unit="days"))/365.25)


# By arm

## Find actual observed probabilities by QRISK2 category

### MI

mi_dpp4su <- survfit(Surv(mi_censtime_yrs, mi_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_dpp4su=1-estimate,
         lower_ci_dpp4su=1-conf.low,
         upper_ci_dpp4su=1-conf.high) %>%
  select(observed_dpp4su, lower_ci_dpp4su, upper_ci_dpp4su, strata)

mi_sglt2 <- survfit(Surv(mi_censtime_yrs, mi_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_sglt2=1-estimate,
         lower_ci_sglt2=1-conf.low,
         upper_ci_sglt2=1-conf.high) %>%
  select(observed_sglt2, lower_ci_sglt2, upper_ci_sglt2, strata)

obs_v_pred <- rbind(mi_dpp4su, mi_sglt2) %>%
  mutate(observed=coalesce(observed_dpp4su, observed_sglt2),
         lower_ci=coalesce(lower_ci_dpp4su, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_dpp4su, upper_ci_sglt2),
         qrisk2_cat=case_when(
           strata=="qrisk2_cat=1" ~ 1,
           strata=="qrisk2_cat=2" ~ 2,
           strata=="qrisk2_cat=3" ~ 3,
           strata=="qrisk2_cat=4" ~ 4,
           strata=="qrisk2_cat=5" ~ 5,
           strata=="qrisk2_cat=6" ~ 6,
           strata=="qrisk2_cat=7" ~ 7,
           strata=="qrisk2_cat=8" ~ 8,
           strata=="qrisk2_cat=9" ~ 9,
           strata=="qrisk2_cat=10" ~ 10
         ),
         studydrug=ifelse(is.na(observed_sglt2), "DPP4SU", "SGLT2"))

dpp4su_events <- cohort %>%
  filter(studydrug=="DPP4SU" & mi_censvar==1) %>%
  group_by(qrisk2_cat) %>%
  summarise(DPP4SU=n())

sglt2_events <- cohort %>%
  filter(studydrug=="SGLT2" & mi_censvar==1) %>%
  group_by(qrisk2_cat) %>%
  summarise(SGLT2=n())

events_table <- data.frame(t(dpp4su_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="qrisk2_cat")


dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(qrisk2_cat==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_pred=NA, qrisk2_cat=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=qrisk2_cat, group=studydrug, color=studydrug)) + 
  geom_point(aes(y = mean_pred*100), position=dodge) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1,width=0.25,size=1, position=dodge) +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,10,by=2)), limits=c(0,11)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("MI incidence")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))



### Stroke

stroke_dpp4su <- survfit(Surv(stroke_censtime_yrs, stroke_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_dpp4su=1-estimate,
         lower_ci_dpp4su=1-conf.low,
         upper_ci_dpp4su=1-conf.high) %>%
  select(observed_dpp4su, lower_ci_dpp4su, upper_ci_dpp4su, strata)

stroke_sglt2 <- survfit(Surv(stroke_censtime_yrs, stroke_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_sglt2=1-estimate,
         lower_ci_sglt2=1-conf.low,
         upper_ci_sglt2=1-conf.high) %>%
  select(observed_sglt2, lower_ci_sglt2, upper_ci_sglt2, strata)

obs_v_pred <- rbind(stroke_dpp4su, stroke_sglt2) %>%
  mutate(observed=coalesce(observed_dpp4su, observed_sglt2),
         lower_ci=coalesce(lower_ci_dpp4su, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_dpp4su, upper_ci_sglt2),
         qrisk2_cat=case_when(
           strata=="qrisk2_cat=1" ~ 1,
           strata=="qrisk2_cat=2" ~ 2,
           strata=="qrisk2_cat=3" ~ 3,
           strata=="qrisk2_cat=4" ~ 4,
           strata=="qrisk2_cat=5" ~ 5,
           strata=="qrisk2_cat=6" ~ 6,
           strata=="qrisk2_cat=7" ~ 7,
           strata=="qrisk2_cat=8" ~ 8,
           strata=="qrisk2_cat=9" ~ 9,
           strata=="qrisk2_cat=10" ~ 10
         ),
         studydrug=ifelse(is.na(observed_sglt2), "DPP4SU", "SGLT2"))

dpp4su_events <- cohort %>%
  filter(studydrug=="DPP4SU" & stroke_censvar==1) %>%
  group_by(qrisk2_cat) %>%
  summarise(DPP4SU=n())

sglt2_events <- cohort %>%
  filter(studydrug=="SGLT2" & stroke_censvar==1) %>%
  group_by(qrisk2_cat) %>%
  summarise(SGLT2=n())

events_table <- data.frame(t(dpp4su_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="qrisk2_cat")


dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(qrisk2_cat==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_pred=NA, qrisk2_cat=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=qrisk2_cat, group=studydrug, color=studydrug)) + 
  geom_point(aes(y = mean_pred*100), position=dodge) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1,width=0.25,size=1, position=dodge) +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,10,by=2)), limits=c(0,11)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Stroke incidence")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))


### CV death

cv_death_dpp4su <- survfit(Surv(cv_death_censtime_yrs, cv_death_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_dpp4su=1-estimate,
         lower_ci_dpp4su=1-conf.low,
         upper_ci_dpp4su=1-conf.high) %>%
  select(observed_dpp4su, lower_ci_dpp4su, upper_ci_dpp4su, strata)

cv_death_sglt2 <- survfit(Surv(cv_death_censtime_yrs, cv_death_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_sglt2=1-estimate,
         lower_ci_sglt2=1-conf.low,
         upper_ci_sglt2=1-conf.high) %>%
  select(observed_sglt2, lower_ci_sglt2, upper_ci_sglt2, strata)

obs_v_pred <- rbind(cv_death_dpp4su, cv_death_sglt2) %>%
  mutate(observed=coalesce(observed_dpp4su, observed_sglt2),
         lower_ci=coalesce(lower_ci_dpp4su, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_dpp4su, upper_ci_sglt2),
         qrisk2_cat=case_when(
           strata=="qrisk2_cat=1" ~ 1,
           strata=="qrisk2_cat=2" ~ 2,
           strata=="qrisk2_cat=3" ~ 3,
           strata=="qrisk2_cat=4" ~ 4,
           strata=="qrisk2_cat=5" ~ 5,
           strata=="qrisk2_cat=6" ~ 6,
           strata=="qrisk2_cat=7" ~ 7,
           strata=="qrisk2_cat=8" ~ 8,
           strata=="qrisk2_cat=9" ~ 9,
           strata=="qrisk2_cat=10" ~ 10
         ),
         studydrug=ifelse(is.na(observed_sglt2), "DPP4SU", "SGLT2"))

dpp4su_events <- cohort %>%
  filter(studydrug=="DPP4SU" & cv_death_censvar==1) %>%
  group_by(qrisk2_cat) %>%
  summarise(DPP4SU=n())

sglt2_events <- cohort %>%
  filter(studydrug=="SGLT2" & cv_death_censvar==1) %>%
  group_by(qrisk2_cat) %>%
  summarise(SGLT2=n())

events_table <- data.frame(t(dpp4su_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="qrisk2_cat")


dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(qrisk2_cat==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_pred=NA, qrisk2_cat=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=qrisk2_cat, group=studydrug, color=studydrug)) + 
  geom_point(aes(y = mean_pred*100), position=dodge) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1,width=0.25,size=1, position=dodge) +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,10,by=2)), limits=c(0,11)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("CV death incidence")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))



############################################################################################

# 3 Is HR constant with baseline QRISK2?
ddist <- datadist(cohort); options(datadist='ddist')


# Unadjusted (QRISK2 5 year score only)

m2 <- cph(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ studydrug*rcs(qrisk2_5yr_score,5), data = cohort,x=T,y=T)
anova(m2)

describe(cohort$qrisk2_5yr_score)
quantile(cohort$qrisk2_5yr_score, c(.01, .99), na.rm=TRUE)
c1 <- quantile(cohort$qrisk2_5yr_score, .01, na.rm=TRUE)
c99 <- quantile(cohort$qrisk2_5yr_score, .99, na.rm=TRUE)

contrast_spline.1 <- contrast(m2,list(studydrug = "SGLT2", qrisk2_5yr_score = seq(c1,c99,by=0.05)),list(studydrug = "DPP4SU", qrisk2_5yr_score = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('qrisk2_5yr_score','Contrast','Lower','Upper')])

#plot and save
contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=qrisk2_5yr_score, y=exp(Contrast))) + 
  geom_line(data=contrast_spline_df,aes(x=qrisk2_5yr_score, y=exp(Contrast)), size=1) + 
  xlab(expression(paste("QRISK2 5 yr score"))) + 
  ylab("HR") +
  scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=qrisk2_5yr_score,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  geom_hline(yintercept = 0.94, linetype = "twodash", color="red", size=1)  +
  geom_hline(yintercept = 0.83, linetype = "twodash", color="red")  +
  geom_hline(yintercept = 1.07, linetype = "twodash", color="red")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +  
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("")

# define a marginal histogram
marginal_distribution <- function(x,var) {
  ggplot(x, aes_string(x = var)) +
    geom_histogram(bins = 64, alpha = 0.4, position = "identity") +
    guides(fill = FALSE) +
    #theme_void() +
    theme(legend.title = element_blank(), panel.background = element_rect( fill = "white",color = "grey50")) +  
    scale_x_continuous(breaks = seq(0,50,5)) +
    xlab(expression(paste("QRISK2 5 yr score"))) + 
    #theme(plot.margin = margin()) +
    theme(text = element_text(size = 14),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) + 
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(color="grey50"))
  
}


hist.dta <- cohort %>% filter(qrisk2_5yr_score>=c1 &  qrisk2_5yr_score <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "qrisk2_5yr_score")


# Arranging the plot using cowplot
plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
          rel_heights = c(1,0.4), rel_widths = c(1,1))




# Plot with 'deciles' of QRISK2 shown
#plot and save
ggplot(data=contrast_spline_df,aes(x=qrisk2_5yr_score, y=exp(Contrast))) + 
  geom_line(data=contrast_spline_df,aes(x=qrisk2_5yr_score, y=exp(Contrast)), size=1) + 
  xlab(expression(paste("QRISK2 5 yr score"))) + 
  ylab("HR") +
  scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=qrisk2_5yr_score,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  geom_hline(yintercept = 0.94, linetype = "twodash", color="red", size=1)  +
  geom_hline(yintercept = 0.83, linetype = "twodash", color="red")  +
  geom_hline(yintercept = 1.07, linetype = "twodash", color="red")  +
  geom_vline(xintercept = 3, size=1, colour="grey80")  +
  geom_vline(xintercept = 5, size=1, colour="grey80")  +
  geom_vline(xintercept = 6.5, size=1, colour="grey80")  +
  geom_vline(xintercept = 8, size=1, colour="grey80")  +
  geom_vline(xintercept = 9.5, size=1, colour="grey80")  +
  geom_vline(xintercept = 11.5, size=1, colour="grey80")  +
  geom_vline(xintercept = 13.5, size=1, colour="grey80")  +
  geom_vline(xintercept = 16, size=1, colour="grey80")  +
  geom_vline(xintercept = 20, size=1, colour="grey80")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +  
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("")



# Adjusted plot
m3 <- cph(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ studydrug*rcs(qrisk2_5yr_score,5) + dstartdate_age  + malesex + dstartdate_dm_dur_all + imd2015_10 + drugline_all + ncurrtx,
          data = cohort,x=T,y=T)
anova(m3)
#no interaction

new_contrast_spline.1 <- contrast(m3,list(studydrug = "SGLT2", qrisk2_5yr_score = seq(c1,c99,by=0.05)),list(studydrug = "DPP4SU", qrisk2_5yr_score = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
new_contrast_spline_df <- as.data.frame(new_contrast_spline.1[c('qrisk2_5yr_score','Contrast','Lower','Upper')])

#plot and save
contrast_spline_plot_1 <- ggplot(data=new_contrast_spline_df,aes(x=qrisk2_5yr_score, y=exp(Contrast))) + 
  geom_line(data=new_contrast_spline_df,aes(x=qrisk2_5yr_score, y=exp(Contrast)), size=1) + 
  xlab(expression(paste("QRISK2 5 yr score"))) + 
  ylab("HR") +
  scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=new_contrast_spline_df,aes(x=qrisk2_5yr_score,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  geom_hline(yintercept = 0.94, linetype = "twodash", color="red", size=1)  +
  geom_hline(yintercept = 0.83, linetype = "twodash", color="red")  +
  geom_hline(yintercept = 1.07, linetype = "twodash", color="red")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +  
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("")

plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
          rel_heights = c(1,0.4), rel_widths = c(1,1))




## Statins at baseline - in last 91 days
cohort <- cohort %>%
  mutate(baseline_statins=ifelse(!is.na(predrug_latest_statins) & dstartdate-predrug_latest_statins<=91, 1, 0))

cohort %>%
  group_by(qrisk2_cat, studydrug) %>%
  summarise(total=n(),
            baseline_statins_count=sum(baseline_statins, na.rm=TRUE)) %>%
  mutate(baseline_statins_perc=round(baseline_statins_count*100/total, 0))



## Statins in followup
cohort <- cohort %>%
  mutate(statins_in_followup=ifelse(!is.na(postdrug_first_statins) & postdrug_first_statins<mace_broad_censdate, 1, 0))

cohort %>%
  group_by(qrisk2_cat, studydrug) %>%
  summarise(total=n(),
            statins_in_followup_count=sum(statins_in_followup, na.rm=TRUE)) %>%
  mutate(statins_in_followup_perc=round(statins_in_followup_count*100/total, 0))



## HbA1c
cohort %>% group_by(qrisk2_cat, studydrug) %>% summarise(hba1c=median(prehba1c, na.rm=TRUE))


############################################################################################

# 4 Does calibrating help?

# Re-estimate baseline hazard on 20% random sample of control (DPP4SU) arm

## Assign random 20% of DPP4SU arm as calibration cohort and remove from main cohort
set.seed(123)

cal_cohort <- cohort %>%
  filter(studydrug=="DPP4SU") %>%
  slice_sample(prop=0.2)
#18,179

noncal_cohort <- cohort %>%
  anti_join(cal_cohort, by=c("patid", "dstartdate", "studydrug"))
table(noncal_cohort$studydrug)
#SGLT2: 48,279 
#DPP4SU: 72,683  


## Re-estimate baseline hazard for females
### Original: 0.994671821594238
cal_females <- cal_cohort %>%
  filter(malesex==0)

recal_mod <- coxph(Surv(mace_broad_censtime_yrs,mace_broad_censvar)~offset(qrisk2_lin_predictor), data=cal_females)
female_surv <- summary(survfit(recal_mod),time=5)$surv
sprintf("%.15f",female_surv)
# 0.955939959552193


## Re-estimate baseline hazard for males
### Original: 0.989570081233978
cal_males <- cal_cohort %>%
  filter(malesex==1)

recal_mod <- coxph(Surv(mace_broad_censtime_yrs,mace_broad_censvar)~offset(qrisk2_lin_predictor), data=cal_males)
male_surv <- summary(survfit(recal_mod),time=5)$surv
sprintf("%.15f",male_surv)
# 0.939068982235158


## Recalculate QRISK2 and compare to uncalibrated and observed results in calibration cohort
cal_cohort <- cal_cohort %>%
  group_by(malesex) %>%
  mutate(centred_lin_predictor=qrisk2_lin_predictor-mean(qrisk2_lin_predictor)) %>%
  ungroup() %>%
  mutate(qrisk2_5yr_score_cal=ifelse(malesex==1, (1-(male_surv^exp(centred_lin_predictor)))*100, (1-(female_surv^exp(centred_lin_predictor)))*100))

predicted_cal_cohort <- cal_cohort %>%
  group_by(qrisk2_cat) %>%
  summarise(mean_pred=mean(qrisk2_5yr_score)/100,
            mean_pred_cal=mean(qrisk2_5yr_score_cal)/100)

observed_cal_cohort <- survfit(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ qrisk2_cat, data=cal_cohort) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed=1-estimate,
         lower_ci=1-conf.low,
         upper_ci=1-conf.high) %>%
  select(observed, lower_ci, upper_ci, strata)

events_cal_cohort <- cal_cohort %>%
  filter(mace_broad_censvar==1) %>%
  group_by(qrisk2_cat) %>%
  summarise(events=n())


obs_v_pred <- cbind(predicted_cal_cohort, observed_cal_cohort)

events_table <- data.frame(t(events_cal_cohort)) %>%
  rownames_to_column() %>%
  filter(rowname!="qrisk2_cat")


dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(qrisk2_cat==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_pred=NA, mean_pred_cal=NA, qrisk2_cat=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=qrisk2_cat)) + 
  geom_point(aes(y = mean_pred*100), position=dodge, shape=4, color="black", size=2, stroke=2) +
  geom_point(aes(y = mean_pred_cal*100), position=dodge, shape=4, color="blue", size=2, stroke=2) +
  geom_point(aes(y = observed*100), position=dodge, color="darkgreen", alpha=0.5, size=3) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge, color="darkgreen", alpha=0.5) +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,25,by=5)), limits=c(0,26)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("(Un)calibrated QRISK2 vs MACE in calibration cohort (n=18,179)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(6,1))


## C-stats
cal_cohort <- cal_cohort %>%
  mutate(qrisk2_survival=(100-qrisk2_5yr_score)/100,
         qrisk2_survival_cal=(100-qrisk2_5yr_score_cal)/100)

### Uncalibrated
surv_mod <- coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar)~qrisk2_survival,method="breslow",data=cal_cohort)
summary(surv_mod)$concordance[1]
summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2])
summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2])
# 0.6762914 (0.6549061-0.6976768)

### Calibrated
surv_mod <- coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar)~qrisk2_survival_cal,method="breslow",data=cal_cohort)
summary(surv_mod)$concordance[1]
summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2])
summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2])
# 0.676328 (0.6549465-0.6977095)




## Recalculate QRISK2 in rest of cohort and compare to uncalibrated and observed results
noncal_cohort <- noncal_cohort %>%
  group_by(malesex) %>%
  mutate(centred_lin_predictor=qrisk2_lin_predictor-mean(qrisk2_lin_predictor)) %>%
  ungroup() %>%
  mutate(qrisk2_5yr_score_cal=ifelse(malesex==1, (1-(male_surv^exp(centred_lin_predictor)))*100, (1-(female_surv^exp(centred_lin_predictor)))*100))


### DPP4SU
predicted_noncal_dpp4su <- noncal_cohort %>%
  filter(studydrug=="DPP4SU") %>%
  group_by(qrisk2_cat) %>%
  summarise(mean_pred=mean(qrisk2_5yr_score)/100,
            mean_pred_cal=mean(qrisk2_5yr_score_cal)/100)

observed_noncal_dpp4su <- survfit(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ qrisk2_cat, data=noncal_cohort[noncal_cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed=1-estimate,
         lower_ci=1-conf.low,
         upper_ci=1-conf.high) %>%
  select(observed, lower_ci, upper_ci, strata)

events_noncal_dpp4su <- noncal_cohort %>%
  filter(studydrug=="DPP4SU" & mace_broad_censvar==1) %>%
  group_by(qrisk2_cat) %>%
  summarise(events=n())


obs_v_pred <- cbind(predicted_noncal_dpp4su, observed_noncal_dpp4su)

events_table <- data.frame(t(events_noncal_dpp4su)) %>%
  rownames_to_column() %>%
  filter(rowname!="qrisk2_cat")


dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(qrisk2_cat==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_pred=NA, mean_pred_cal=NA, qrisk2_cat=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=qrisk2_cat)) + 
  geom_point(aes(y = mean_pred*100), position=dodge, shape=4, color="black", size=2, stroke=2) +
  geom_point(aes(y = mean_pred_cal*100), position=dodge, shape=4, color="blue", size=2, stroke=2) +
  geom_point(aes(y = observed*100), position=dodge, color="darkgreen", alpha=0.5, size=3) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge, color="darkgreen", alpha=0.5) +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,25,by=5)), limits=c(0,26)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("QRISK2 vs MACE in DPP4SU non-calibration cohort (n=72,683)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(6,1))


## C-stats
noncal_cohort <- noncal_cohort %>%
  mutate(qrisk2_survival=(100-qrisk2_5yr_score)/100,
         qrisk2_survival_cal=(100-qrisk2_5yr_score_cal)/100)

### Uncalibrated
surv_mod <- coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar)~qrisk2_survival,method="breslow",data=noncal_cohort[noncal_cohort$studydrug=="DPP4SU",])
summary(surv_mod)$concordance[1]
summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2])
summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2])
# 0.6756263 (0.6650895-0.6861631)

### Calibrated
surv_mod <- coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar)~qrisk2_survival_cal,method="breslow",data=noncal_cohort[noncal_cohort$studydrug=="DPP4SU",])
summary(surv_mod)$concordance[1]
summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2])
summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2])
# 0.6755896 (0.6650576-0.6861215)



### SGLT2
predicted_noncal_sglt2 <- noncal_cohort %>%
  filter(studydrug=="SGLT2") %>%
  group_by(qrisk2_cat) %>%
  summarise(mean_pred=mean(qrisk2_5yr_score)/100,
            mean_pred_cal=mean(qrisk2_5yr_score_cal)/100)

observed_noncal_sglt2 <- survfit(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ qrisk2_cat, data=noncal_cohort[noncal_cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed=1-estimate,
         lower_ci=1-conf.low,
         upper_ci=1-conf.high) %>%
  select(observed, lower_ci, upper_ci, strata)

events_noncal_sglt2 <- noncal_cohort %>%
  filter(studydrug=="SGLT2" & mace_broad_censvar==1) %>%
  group_by(qrisk2_cat) %>%
  summarise(events=n())


obs_v_pred <- cbind(predicted_noncal_sglt2, observed_noncal_sglt2)

events_table <- data.frame(t(events_noncal_sglt2)) %>%
  rownames_to_column() %>%
  filter(rowname!="qrisk2_cat")


dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(qrisk2_cat==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_pred=NA, mean_pred_cal=NA, qrisk2_cat=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=qrisk2_cat)) + 
  geom_point(aes(y = mean_pred*100), position=dodge, shape=4, color="black", size=2, stroke=2) +
  geom_point(aes(y = mean_pred_cal*100), position=dodge, shape=4, color="blue", size=2, stroke=2) +
  geom_point(aes(y = observed*100), position=dodge, color="darkgreen", alpha=0.5, size=3) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge, color="darkgreen", alpha=0.5) +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,25,by=5)), limits=c(0,26)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("QRISK2 vs MACE in SGLT2 arm (n=48,279)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(6,1))


## C-stats

### Uncalibrated
surv_mod <- coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar)~qrisk2_survival,method="breslow",data=noncal_cohort[noncal_cohort$studydrug=="SGLT2",])
summary(surv_mod)$concordance[1]
summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2])
summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2])
# 0.6729413 (0.6561626-0.68972)

### Calibrated
surv_mod <- coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar)~qrisk2_survival_cal,method="breslow",data=noncal_cohort[noncal_cohort$studydrug=="SGLT2",])
summary(surv_mod)$concordance[1]
summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2])
summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2])
# 0.6730351 (0.6562606-0.6898097)

############################################################################################

# 5 Look at survival benefit from SGLT2s when add hazard ratio from trials meta-analysis
# https://jamanetwork.com/journals/jamacardiology/fullarticle/2771459


# Add HR from trials
noncal_cohort <- noncal_cohort %>%
  mutate(qrisk2_survival_cal_sglt2=qrisk2_survival_cal^0.94,
         qrisk2_5yr_score_cal_sglt2=100-(qrisk2_survival_cal_sglt2*100))


# Does QRISK2+HR perform better than QRISK2 in SGLT2 arm?
predicted_sglt2 <- noncal_cohort %>%
  filter(studydrug=="SGLT2") %>%
  group_by(qrisk2_cat) %>%
  summarise(mean_pred_cal=mean(qrisk2_5yr_score_cal)/100,
            mean_pred_cal_sglt2=mean(qrisk2_5yr_score_cal_sglt2)/100)

observed_sglt2 <- survfit(Surv(mace_broad_censtime_yrs, mace_broad_censvar) ~ qrisk2_cat, data=noncal_cohort[noncal_cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed=1-estimate,
         lower_ci=1-conf.low,
         upper_ci=1-conf.high) %>%
  select(observed, lower_ci, upper_ci, strata)

events_sglt2 <- noncal_cohort %>%
  filter(studydrug=="SGLT2" & mace_broad_censvar==1) %>%
  group_by(qrisk2_cat) %>%
  summarise(events=n())


obs_v_pred <- cbind(predicted_sglt2, observed_sglt2)

events_table <- data.frame(t(events_sglt2)) %>%
  rownames_to_column() %>%
  filter(rowname!="qrisk2_cat")


dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(qrisk2_cat==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_pred=NA, mean_pred_cal=NA, qrisk2_cat=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=qrisk2_cat)) + 
  geom_point(aes(y = mean_pred_cal_sglt2*100), position=dodge, shape=4, color="red", size=2, stroke=2) +
  geom_point(aes(y = mean_pred_cal*100), position=dodge, shape=4, color="blue", size=2, stroke=2) +
  geom_point(aes(y = observed*100), position=dodge, color="darkgreen", alpha=0.5, size=3) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge, color="darkgreen", alpha=0.5) +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,25,by=5)), limits=c(0,26)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("QRISK2 vs MACE in SGLT2 arm (n=48,279)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(6,1))


## C-stats

### QRISK2 alone
surv_mod <- coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar)~qrisk2_survival_cal,method="breslow",data=noncal_cohort[noncal_cohort$studydrug=="SGLT2",])
summary(surv_mod)$concordance[1]
summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2])
summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2])
# 0.6730351 (0.6562606-0.6898097)

### QRISK2 + HR
surv_mod <- coxph(Surv(mace_broad_censtime_yrs, mace_broad_censvar)~qrisk2_survival_cal_sglt2,method="breslow",data=noncal_cohort[noncal_cohort$studydrug=="SGLT2",])
summary(surv_mod)$concordance[1]
summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2])
summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2])
# 0.6730351 (0.6562606-0.6898097)
  
  
