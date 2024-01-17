
# Exploring QRISK2 and QDiabetes-Heart Failure (QDHF) for predicting heart failure (HF) incidence

## 1 Uncalibrated 5-year risk scores vs HF incidence in both arms (report C-stat + number of events)

## 2 (Uncalibrated) 5-year QRISK2 / QDHF vs hazard ratio for HF (is HR constant by baseline QRISK2/QDHF? Report ANOVA)

## 3 Calibrate 5-year risk scores on 20% sample of control arm (DPP4SU) - report C-stat before and after and plot
## Recalculate 5-year QRISK2 for control and study arm (report C-stat before and after and plot)

## 4 Plot unadjusted predicted benefit vs actual benefit by benefit decile/categories
## Plot adjusted version of above (adjust for...)

### 5 Number needed to treat per benefit decile/category for predicted and actual benefit

# Not including GLP1 for now

############################################################################################

# Setup
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(patchwork)
library(rms)
library(cowplot)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# 1 Cohort selection and variable setup

## A Cohort selection (see cohort_definition function for details)

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Raw data/")
load("20230308_t2d_1stinstance.Rda")
load("20230308_t2d_all_drug_periods.Rda")

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Scripts/Functions")
source("cohort_definition.R")

cohort <- define_cohort(t2d_1stinstance, t2d_all_drug_periods)

cohort <- cohort %>% filter(studydrug!="GLP1")

table(cohort$studydrug)
# DPP4SU 95362
# SGLT2 50234


## B Make variables for survival analysis of all endpoints (see survival_variables function for details)

source("survival_variables.R")

cohort <- add_surv_vars(cohort, main_only=TRUE)


## C Just keep variables of interest

cohort <- cohort %>%
  
  select(patid, malesex, ethnicity_5cat_decoded, imd2015_10, regstartdate, gp_record_end, death_date, drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, prehba1c, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preegfr, preckdstage, qrisk2_smoking_cat, contains("cens"), qrisk2_lin_predictor, qrisk2_5yr_score, qdiabeteshf_lin_predictor, qdiabeteshf_5yr_score, starts_with("ckdpc"), last_sglt2_stop, contains("statins"))

rm(list=setdiff(ls(), "cohort"))


############################################################################################

# 1 How well do uncalibrated QRISK2 and QDHF predict 5-year HF incidence?


## A Compare QRISK2 and QDHF distributions

scaled_qrisk2 <- data.frame(scale(cohort %>% select(qrisk2_5yr_score)))
scaled_qdhf <- data.frame(scale(cohort %>% select(qdiabeteshf_5yr_score)))

test <- cbind(scaled_qrisk2, scaled_qdhf) %>%
  pivot_longer(cols=c(qrisk2_5yr_score, qdiabeteshf_5yr_score))

ggplot(test, aes(x=value, fill=name)) +
  geom_density(alpha=.3) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))

ggplot(cohort, aes(x=qrisk2_5yr_score)) + 
  geom_point(aes(y=qdiabeteshf_5yr_score))



## B Use ~QRISK2 deciles - define these

quantile(cohort$qrisk2_5yr_score, probs = seq(.1, .9, by = .1))
# 3.3, 5, 6.5, 8, 9.6, 11.4, 13.5, 16.2, 20.4

cohort <- cohort %>%
  mutate(qrisk2_cat=case_when(
    qrisk2_5yr_score<=3.3 ~ 1,
    qrisk2_5yr_score>3.3 & qrisk2_5yr_score<=5 ~ 2,
    qrisk2_5yr_score>5 & qrisk2_5yr_score<=6.5 ~ 3,
    qrisk2_5yr_score>6.5 & qrisk2_5yr_score<=8 ~ 4,
    qrisk2_5yr_score>8 & qrisk2_5yr_score<=9.6 ~ 5,
    qrisk2_5yr_score>9.6 & qrisk2_5yr_score<=11.4 ~ 6,
    qrisk2_5yr_score>11.4 & qrisk2_5yr_score<=13.5 ~ 7,
    qrisk2_5yr_score>13.5 & qrisk2_5yr_score<=16.2 ~ 8,
    qrisk2_5yr_score>16.2 & qrisk2_5yr_score<=20.4 ~ 9,
    qrisk2_5yr_score>20.4 ~ 10
  ))

table(cohort$qrisk2_cat, useNA="always")
# No NAs and roughly equal counts



## C Compare predicted and observed results

# Get mean predicted probabilities from QRISK2 and QDHF for each QRISK2 category by studydrug
predicted <- cohort %>%
  group_by(qrisk2_cat, studydrug) %>%
  summarise(mean_qrisk2_pred=mean(qrisk2_5yr_score)/100,
            mean_qdhf_pred=mean(qdiabeteshf_5yr_score)/100)


# Find actual observed probabilities by QRISK2 category and studydrug
observed_dpp4su <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_dpp4su=1-estimate,
         lower_ci_dpp4su=1-conf.low,
         upper_ci_dpp4su=1-conf.high) %>%
  select(observed_dpp4su, lower_ci_dpp4su, upper_ci_dpp4su, strata)

dpp4su_events <- cohort %>%
  filter(studydrug=="DPP4SU" & hf_censvar==1) %>%
  group_by(qrisk2_cat) %>%
  summarise(DPP4SU=n())

observed_sglt2 <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ qrisk2_cat, data=cohort[cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_sglt2=1-estimate,
         lower_ci_sglt2=1-conf.low,
         upper_ci_sglt2=1-conf.high) %>%
  select(observed_sglt2, lower_ci_sglt2, upper_ci_sglt2, strata)

sglt2_events <- cohort %>%
  filter(studydrug=="SGLT2" & hf_censvar==1) %>%
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
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_qrisk2_pred=NA, mean_qdhf_pred=NA, qrisk2_cat=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=qrisk2_cat, group=studydrug, color=studydrug, fill=studydrug)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_qrisk2_pred*100), position=dodge, shape=19, size=2) +
  geom_point(aes(y = mean_qdhf_pred*100), position=dodge, shape=22, size=2, color="black") +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,25,by=5)), limits=c(0,27)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Uncalibrated QRISK2/QDHF vs HF incidence (5 year)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))

# QRISK2 way overestimates but might expect this as predicting MACE = much more common
# QDHF looks good


## D Look at C-stats

cohort <- cohort %>%
  mutate(qrisk2_survival=(100-qrisk2_5yr_score)/100,
         qdiabeteshf_survival=(100-qdiabeteshf_5yr_score)/100)
  
### QRISK2
surv_mod <- coxph(Surv(hf_censtime_yrs, hf_censvar)~qrisk2_survival,data=cohort,method="breslow")
summary(surv_mod)$concordance[1]
summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2])
summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2])
# 0.7139811 (0.7044466-0.7235156)

### QDHF
surv_mod <- coxph(Surv(hf_censtime_yrs, hf_censvar)~qdiabeteshf_survival,data=cohort,method="breslow")
summary(surv_mod)$concordance[1]
summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2])
summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2])
# 0.7244267 (0.7147242-0.7341292)


############################################################################################

# # 2 Is HR constant with baseline QRISK2 and baseline QDHF?
# 
# ddist <- datadist(cohort); options(datadist='ddist')
# 
# 
# # Unadjusted QRISK2
# 
# m <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug*rcs(qrisk2_5yr_score,5), data = cohort,x=T,y=T)
# anova(m)
# 
# describe(cohort$qrisk2_5yr_score)
# quantile(cohort$qrisk2_5yr_score, c(.01, .99), na.rm=TRUE)
# c1 <- quantile(cohort$qrisk2_5yr_score, .01, na.rm=TRUE)
# c99 <- quantile(cohort$qrisk2_5yr_score, .99, na.rm=TRUE)
# 
# contrast_spline.1 <- contrast(m,list(studydrug = "SGLT2", qrisk2_5yr_score = seq(c1,c99,by=0.05)),list(studydrug = "DPP4SU", qrisk2_5yr_score = seq(c1,c99,by=0.05)))
# # save the contrast calculations in a dataframe
# contrast_spline_df <- as.data.frame(contrast_spline.1[c('qrisk2_5yr_score','Contrast','Lower','Upper')])
# 
# contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=qrisk2_5yr_score, y=exp(Contrast))) +
#   geom_line(data=contrast_spline_df,aes(x=qrisk2_5yr_score, y=exp(Contrast)), size=1) +
#   xlab(expression(paste("QRISK2 5 yr score"))) +
#   ylab("HR") +
#   scale_x_continuous(breaks = seq(0,50,5)) +
#   #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
#   geom_ribbon(data=contrast_spline_df,aes(x=qrisk2_5yr_score,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
#   geom_hline(yintercept = 1, linetype = "dashed")  +
#   geom_hline(yintercept = 0.63, linetype = "twodash", color="red", size=1)  +
#   geom_hline(yintercept = 0.50, linetype = "twodash", color="red")  +
#   geom_hline(yintercept = 0.80, linetype = "twodash", color="red")  +
#   theme(legend.position=c(0.8, 0.1)) +
#   theme(legend.title = element_blank()) +
#   theme_bw() +
#   theme(text = element_text(size = 14),
#         axis.line = element_line(colour =  "grey50" ),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   theme(legend.text = element_text(colour="black", size=rel(1))) +
#   ggtitle("")
# 
# # define a marginal histogram
# marginal_distribution <- function(x,var) {
#   ggplot(x, aes_string(x = var)) +
#     geom_histogram(bins = 64, alpha = 0.4, position = "identity") +
#     guides(fill = FALSE) +
#     #theme_void() +
#     theme(legend.title = element_blank(), panel.background = element_rect( fill = "white",color = "grey50")) +
#     scale_x_continuous(breaks = seq(0,50,5)) +
#     xlab(expression(paste("QRISK2 5 yr score"))) +
#     #theme(plot.margin = margin()) +
#     theme(text = element_text(size = 14),
#           axis.ticks.y = element_blank(),
#           axis.text.y = element_blank(),
#           axis.title.y = element_blank()) +
#     theme(axis.ticks.y = element_blank(),
#           axis.text.y = element_blank(),
#           axis.title.y = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank(),
#           panel.background = element_blank(),
#           axis.line.x = element_line(color="grey50"))
# 
# }
# 
# 
# hist.dta <- cohort %>% filter(qrisk2_5yr_score>=c1 &  qrisk2_5yr_score <= c99)
# hist.dta$dummy <- 1
# x_hist <- marginal_distribution(hist.dta, "qrisk2_5yr_score")
# 
# 
# # Arranging the plot using cowplot
# plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
#           rel_heights = c(1,0.4), rel_widths = c(1,1))
# 
# 
# 
# 
# # Plot with 'deciles' of QRISK2 shown
# #plot and save
# ggplot(data=contrast_spline_df,aes(x=qrisk2_5yr_score, y=exp(Contrast))) +
#   geom_line(data=contrast_spline_df,aes(x=qrisk2_5yr_score, y=exp(Contrast)), size=1) +
#   xlab(expression(paste("QRISK2 5 yr score"))) +
#   ylab("HR") +
#   scale_x_continuous(breaks = seq(0,50,5)) +
#   #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
#   geom_ribbon(data=contrast_spline_df,aes(x=qrisk2_5yr_score,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
#   geom_hline(yintercept = 1, linetype = "dashed")  +
#   geom_hline(yintercept = 0.94, linetype = "twodash", color="red", size=1)  +
#   geom_hline(yintercept = 0.83, linetype = "twodash", color="red")  +
#   geom_hline(yintercept = 1.07, linetype = "twodash", color="red")  +
#   geom_vline(xintercept = 3, size=1, colour="grey80")  +
#   geom_vline(xintercept = 5, size=1, colour="grey80")  +
#   geom_vline(xintercept = 6.5, size=1, colour="grey80")  +
#   geom_vline(xintercept = 8, size=1, colour="grey80")  +
#   geom_vline(xintercept = 9.5, size=1, colour="grey80")  +
#   geom_vline(xintercept = 11.5, size=1, colour="grey80")  +
#   geom_vline(xintercept = 13.5, size=1, colour="grey80")  +
#   geom_vline(xintercept = 16, size=1, colour="grey80")  +
#   geom_vline(xintercept = 20, size=1, colour="grey80")  +
#   theme(legend.position=c(0.8, 0.1)) +
#   theme(legend.title = element_blank()) +
#   theme_bw() +
#   theme(text = element_text(size = 14),
#         axis.line = element_line(colour =  "grey50" ),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   theme(legend.text = element_text(colour="black", size=rel(1))) +
#   ggtitle("")
# 
# 
# 
# # Adjusted plot
# m3 <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug*rcs(qrisk2_5yr_score,5) + dstartdate_age  + malesex + dstartdate_dm_dur_all + imd2015_10 + drugline_all + ncurrtx,
#           data = cohort,x=T,y=T)
# anova(m3)
# #no interaction
# 
# new_contrast_spline.1 <- contrast(m3,list(studydrug = "SGLT2", qrisk2_5yr_score = seq(c1,c99,by=0.05)),list(studydrug = "DPP4SU", qrisk2_5yr_score = seq(c1,c99,by=0.05)))
# # save the contrast calculations in a dataframe
# new_contrast_spline_df <- as.data.frame(new_contrast_spline.1[c('qrisk2_5yr_score','Contrast','Lower','Upper')])
# 
# #plot and save
# contrast_spline_plot_1 <- ggplot(data=new_contrast_spline_df,aes(x=qrisk2_5yr_score, y=exp(Contrast))) +
#   geom_line(data=new_contrast_spline_df,aes(x=qrisk2_5yr_score, y=exp(Contrast)), size=1) +
#   xlab(expression(paste("QRISK2 5 yr score"))) +
#   ylab("HR") +
#   scale_x_continuous(breaks = seq(0,50,5)) +
#   #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
#   geom_ribbon(data=new_contrast_spline_df,aes(x=qrisk2_5yr_score,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
#   geom_hline(yintercept = 1, linetype = "dashed")  +
#   geom_hline(yintercept = 0.63, linetype = "twodash", color="red", size=1)  +
#   geom_hline(yintercept = 0.50, linetype = "twodash", color="red")  +
#   geom_hline(yintercept = 0.80, linetype = "twodash", color="red")  +
#   theme(legend.position=c(0.8, 0.1)) +
#   theme(legend.title = element_blank()) +
#   theme_bw() +
#   theme(text = element_text(size = 14),
#         axis.line = element_line(colour =  "grey50" ),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   theme(legend.text = element_text(colour="black", size=rel(1))) +
#   ggtitle("")
# 
# plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
#           rel_heights = c(1,0.4), rel_widths = c(1,1))


############################################################################################

# 3 Does calibrating help?

# Assign random 20% of DPP4SU arm as calibration cohort and remove from main cohort
set.seed(123)

cal_cohort <- cohort %>%
  filter(studydrug=="DPP4SU") %>%
  slice_sample(prop=0.2)
#19,072

noncal_cohort <- cohort %>%
  anti_join(cal_cohort, by=c("patid", "dstartdate", "studydrug"))
table(noncal_cohort$studydrug)
#SGLT2: 50,234
#DPP4SU: 76,290


# QRISK2

## Females
### Original: 0.994671821594238
cal_females <- cal_cohort %>% filter(malesex==0)
recal_mod <- coxph(Surv(hf_censtime_yrs,hf_censvar)~offset(qrisk2_lin_predictor), data=cal_females)
female_qrisk2_surv <- summary(survfit(recal_mod),time=5)$surv
sprintf("%.15f", female_qrisk2_surv)
# 0.961228613501420 - unsurprisingly much lower

## Males
### Original: 0.989570081233978
cal_males <- cal_cohort %>% filter(malesex==1)
recal_mod <- coxph(Surv(hf_censtime_yrs,hf_censvar)~offset(qrisk2_lin_predictor), data=cal_males)
male_qrisk2_surv <- summary(survfit(recal_mod),time=5)$surv
sprintf("%.15f", male_qrisk2_surv)
# 0.961097423064088 - unsurprisingly much lower


# QDiabetes-HF

## Females
### Original: 0.985216200351715
cal_females <- cal_cohort %>% filter(malesex==0)
recal_mod <- coxph(Surv(hf_censtime_yrs,hf_censvar)~offset(qdiabeteshf_lin_predictor), data=cal_females)
female_qdhf_surv <- summary(survfit(recal_mod),time=5)$surv
sprintf("%.15f", female_qdhf_surv)
# 0.965202142600253

## Males
### Original: 0.981611728668213
cal_males <- cal_cohort %>% filter(malesex==1)
recal_mod <- coxph(Surv(hf_censtime_yrs,hf_censvar)~offset(qdiabeteshf_lin_predictor), data=cal_males)
male_qdhf_surv <- summary(survfit(recal_mod),time=5)$surv
sprintf("%.15f", male_qdhf_surv)
# 0.965209400040619



# Recalculate scores in rest of cohort

noncal_cohort <- noncal_cohort %>%
  group_by(malesex) %>%
  mutate(centred_qrisk2_lin_predictor=qrisk2_lin_predictor-mean(qrisk2_lin_predictor),
         centred_qdhf_lin_predictor=qdiabeteshf_lin_predictor-mean(qdiabeteshf_lin_predictor)) %>%
  ungroup() %>%
  mutate(qrisk2_5yr_score_cal=ifelse(malesex==1, (1-(male_qrisk2_surv^exp(centred_qrisk2_lin_predictor)))*100, (1-(female_qrisk2_surv^exp(centred_qrisk2_lin_predictor)))*100),
         qdiabeteshf_5yr_score_cal=ifelse(malesex==1, (1-(male_qdhf_surv^exp(centred_qdhf_lin_predictor)))*100, (1-(female_qdhf_surv^exp(centred_qdhf_lin_predictor)))*100))


# Plot calibrated + uncalibrated + observed scores

# Get mean predicted probabilities for calibrated and uncalibrated
predicted <- noncal_cohort %>%
  group_by(qrisk2_cat, studydrug) %>%
  summarise(mean_qrisk2_pred=mean(qrisk2_5yr_score)/100,
            mean_qrisk2_pred_cal=mean(qrisk2_5yr_score_cal)/100,
            mean_qdhf_pred=mean(qdiabeteshf_5yr_score)/100,
            mean_qdhf_pred_cal=mean(qdiabeteshf_5yr_score_cal)/100)


# Find actual observed probabilities by QRISK2 category and studydrug
observed_dpp4su <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ qrisk2_cat, data=noncal_cohort[noncal_cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_dpp4su=1-estimate,
         lower_ci_dpp4su=1-conf.low,
         upper_ci_dpp4su=1-conf.high) %>%
  select(observed_dpp4su, lower_ci_dpp4su, upper_ci_dpp4su, strata)

dpp4su_events <- cohort %>%
  filter(studydrug=="DPP4SU" & hf_censvar==1) %>%
  group_by(qrisk2_cat) %>%
  summarise(DPP4SU=n())

observed_sglt2 <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ qrisk2_cat, data=noncal_cohort[noncal_cohort$studydrug=="SGLT2",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed_sglt2=1-estimate,
         lower_ci_sglt2=1-conf.low,
         upper_ci_sglt2=1-conf.high) %>%
  select(observed_sglt2, lower_ci_sglt2, upper_ci_sglt2, strata)

sglt2_events <- cohort %>%
  filter(studydrug=="SGLT2" & hf_censvar==1) %>%
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
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_qrisk2_pred=NA, mean_qrisk2_pred_cal=NA, mean_qdhf_pred=NA, mean_qdhf_pred_cal=NA, qrisk2_cat=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=qrisk2_cat, group=studydrug, color=studydrug)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_qrisk2_pred*100), position=dodge, shape=17, size=2) +
  geom_point(aes(y = mean_qrisk2_pred_cal*100), position=dodge, shape=24, size=2, color="black") +
  geom_point(aes(y = mean_qdhf_pred*100), position=dodge, shape=15, size=2) +
  geom_point(aes(y = mean_qdhf_pred_cal*100), position=dodge, shape=22, size=2, color="black") +
  theme_bw() +
  xlab("QRISK2 category") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,30,by=5)), limits=c(0,31)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Uncalibrated/calibrated QRISK2/QDHF vs HF incidence (5 year)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))


############################################################################################

# 4 Look at survival benefit from SGLT2s when add hazard ratio from trials meta-analysis
# https://jamanetwork.com/journals/jamacardiology/fullarticle/2771459

# Add HR from trials
noncal_cohort <- noncal_cohort %>%
  mutate(qrisk2_survival_cal=(100-qrisk2_5yr_score_cal)/100,
         qrisk2_survival_cal_sglt2=qrisk2_survival_cal^0.63,
         qrisk2_5yr_score_cal_sglt2=100-(qrisk2_survival_cal_sglt2*100),
         
         qdhf_survival_cal=(100-qdiabeteshf_5yr_score_cal)/100,
         qdhf_survival_cal_sglt2=qdhf_survival_cal^0.63,
         qdiabeteshf_5yr_score_cal_sglt2=100-(qdhf_survival_cal_sglt2*100))


# Distribution of predicted benefits
noncal_cohort <- noncal_cohort %>%
  mutate(qrisk2_sglt2_benefit=qrisk2_survival_cal_sglt2-qrisk2_survival_cal,
         qdhf_sglt2_benefit=qdhf_survival_cal_sglt2-qdhf_survival_cal)

describe(noncal_cohort$qrisk2_sglt2_benefit)
#mean=1.8%

describe(noncal_cohort$qdhf_sglt2_benefit)
#mean=1.7%

test <- noncal_cohort %>% select(qrisk2_sglt2_benefit, qdhf_sglt2_benefit) %>% pivot_longer(cols=c(qrisk2_sglt2_benefit, qdhf_sglt2_benefit))

ggplot(test, aes(x=value*100, color=name, fill=name)) + 
  geom_histogram(binwidth=0.05, alpha=0.5, position="identity") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))





# Quartile plots

# Cut into quartiles of absolute risk
noncal_cohort$hf.risk.quartile <- ntile(noncal_cohort$qdiabeteshf_5yr_score_cal, 4)

## Average benefit per quartile
noncal_cohort %>%
  group_by(hf.risk.quartile) %>%
  summarise(pred.mean=mean(qdiabeteshf_5yr_score_cal,na.rm=T),
            pred.median=median(qdiabeteshf_5yr_score_cal,na.rm=T),
            sd=sd(qdiabeteshf_5yr_score_cal,na.rm=T),
            l_iqr=quantile(qdiabeteshf_5yr_score_cal,na.rm=T,probs=0.25),
            u_iqr=quantile(qdiabeteshf_5yr_score_cal,na.rm=T,probs=0.75),
            min=min(qdiabeteshf_5yr_score_cal,na.rm=T),
            max=max(qdiabeteshf_5yr_score_cal,na.rm=T))


### Quartile 4
noncal_cohort$studydrug=factor(noncal_cohort$studydrug, levels=c("SGLT2","DPP4SU"))
sfit <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, data = noncal_cohort, subset=hf.risk.quartile==4)
plot4 <- ggsurvplot(sfit,
                    fun = function(x) {100 - x*100},             # Cumulative probability plot
                    censor = F,
                    size= 1.5,
                    legend.title = "" ,
                    legend.labs = c("SGLT2","DPP4SU"),
                    #legend.title = "",
                    risk.table = F,       # show risk table.
                    conf.int = T,         # show confidence intervals for
                    # point estimates of survival curves.
                    #xlim = c(0,0.06),         # present narrower X axis, but not affect
                    # survival estimates.
                    ylim = c(0,21),
                    xlab = "Years",   # customize X axis label.
                    ylab = "Heart failure incidence (%)",   # customize X axis label.
                    break.time.by = 1,     # break X axis in time intervals every two years.
                    ggtheme = theme_classic(), # customize plot and risk table with a theme.
                    risk.table.y.text.col = TRUE, # colour risk table text annotations.
                    #fontsize = 10, # in legend of risk table
                    font.legend=16,
                    font.x = c(18),font.y = c(18),font.tickslab = c(18),
                    axes.offset = FALSE, # start the axis at the origin
                    #palette = c("#998ec3","#f1a340"),
                    #linetype = c("strata"),
                    #legend = "none",
                    legend = "right", 
                    #title = "Predicted 5yr MACE risk 12.7%",
                    #tables.theme = theme_cleantable() # clean theme for tables
                    
) 
plot4




# Compare predicted and actual benefits by decile (predicted = mean of decile, as previous)

## Unadjusted

## QRISK2

noncal_cohort$sglt2_benefit_decile <- ntile(noncal_cohort$qrisk2_sglt2_benefit, 10)
noncal_cohort <- noncal_cohort %>% mutate(sglt2_benefit_decile=as.factor(sglt2_benefit_decile))

pred <- noncal_cohort %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(mean_qrisk2_predicted_benefit=mean(qrisk2_sglt2_benefit, na.rm=T))


noncal_cohort <- noncal_cohort %>% mutate(studydrug=as.vector(studydrug))
ddist <- datadist(noncal_cohort)

model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug*sglt2_benefit_decile, data=noncal_cohort, x=TRUE, y=TRUE, surv=TRUE)

survival_est <- survest(model, newdata=expand.grid(studydrug=c("SGLT2","DPP4SU"), sglt2_benefit_decile=c(1:10)), times=5)


obs <- data.frame(surv=unlist(survival_est$surv), studydrug=rep(c("SGLT2","DPP4SU"),5), sglt2_benefit_decile=rep(1:10, rep_len(2, 10)))

obs <- obs %>%
  pivot_wider(id_cols=sglt2_benefit_decile, values_from=surv, names_from=studydrug) %>%
  mutate(surv_diff=SGLT2-DPP4SU) %>%
  select(sglt2_benefit_decile, surv_diff)


se <- data.frame(se=unlist(survival_est$std.err), studydrug=rep(c("SGLT2","DPP4SU"),5), sglt2_benefit_decile=rep(1:10, rep_len(2, 10)))

se <- se %>%
  pivot_wider(id_cols=sglt2_benefit_decile, values_from=se, names_from=studydrug) %>%
  mutate(se=sqrt((SGLT2^2)+(DPP4SU^2))) %>%
  select(sglt2_benefit_decile, se)

obs <- obs %>%
  inner_join(se, by="sglt2_benefit_decile") %>%
  mutate(lower_ci=surv_diff-(1.96*se),
         upper_ci=surv_diff+(1.96*se))

dpp4su_events <- noncal_cohort %>%
  filter(studydrug=="DPP4SU" & mace_censvar==1) %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(DPP4SU=n())

sglt2_events <- noncal_cohort %>%
  filter(studydrug=="SGLT2" & mace_censvar==1) %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(SGLT2=n())

events_table <- data.frame(t(dpp4su_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="sglt2_benefit_decile")

obs <- obs %>% mutate(sglt2_benefit_decile=as.factor(sglt2_benefit_decile))
obs_v_pred <- pred %>% inner_join(obs, by="sglt2_benefit_decile")


empty_tick <- obs_v_pred %>%
  filter(sglt2_benefit_decile==1) %>%
  mutate(mean_predicted_benefit=NA, surv_diff=NA, lower_ci=NA, upper_ci=NA, sglt2_benefit_decile=as.factor(0))

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_qrisk2_predicted_benefit*100)) + 
  geom_point(aes(y = surv_diff*100)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1) +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  xlab("Predicted SGLT2 benefit (mean by predicted decile)") + ylab("Observed SGLT2 benefit (by predicted decile)") +
  scale_x_continuous(breaks=c(seq(0,5,by=1)), limits=c(0,5)) +
  scale_y_continuous(breaks=c(seq(-2,6,by=1)), limits=c(-2,6)) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5)))

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))

#NNT
obs_v_pred <- obs_v_pred %>%
  mutate(nnt_predicted=1/(mean_qdhf_predicted_benefit*100),
         nnt_observed=1/(surv_diff*100))
obs_v_pred



## QDHF

noncal_cohort$sglt2_benefit_decile <- ntile(noncal_cohort$qdhf_sglt2_benefit, 10)
noncal_cohort <- noncal_cohort %>% mutate(sglt2_benefit_decile=as.factor(sglt2_benefit_decile))

pred <- noncal_cohort %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(mean_qdhf_predicted_benefit=mean(qdhf_sglt2_benefit, na.rm=T))


noncal_cohort <- noncal_cohort %>% mutate(studydrug=as.vector(studydrug))
ddist <- datadist(noncal_cohort)

model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug*sglt2_benefit_decile, data=noncal_cohort, x=TRUE, y=TRUE, surv=TRUE)

survival_est <- survest(model, newdata=expand.grid(studydrug=c("SGLT2","DPP4SU"), sglt2_benefit_decile=c(1:10)), times=5)


obs <- data.frame(surv=unlist(survival_est$surv), studydrug=rep(c("SGLT2","DPP4SU"),5), sglt2_benefit_decile=rep(1:10, rep_len(2, 10)))

obs <- obs %>%
  pivot_wider(id_cols=sglt2_benefit_decile, values_from=surv, names_from=studydrug) %>%
  mutate(surv_diff=SGLT2-DPP4SU) %>%
  select(sglt2_benefit_decile, surv_diff)


se <- data.frame(se=unlist(survival_est$std.err), studydrug=rep(c("SGLT2","DPP4SU"),5), sglt2_benefit_decile=rep(1:10, rep_len(2, 10)))

se <- se %>%
  pivot_wider(id_cols=sglt2_benefit_decile, values_from=se, names_from=studydrug) %>%
  mutate(se=sqrt((SGLT2^2)+(DPP4SU^2))) %>%
  select(sglt2_benefit_decile, se)

obs <- obs %>%
  inner_join(se, by="sglt2_benefit_decile") %>%
  mutate(lower_ci=surv_diff-(1.96*se),
         upper_ci=surv_diff+(1.96*se))

dpp4su_events <- noncal_cohort %>%
  filter(studydrug=="DPP4SU" & mace_censvar==1) %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(DPP4SU=n())

sglt2_events <- noncal_cohort %>%
  filter(studydrug=="SGLT2" & mace_censvar==1) %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(SGLT2=n())

events_table <- data.frame(t(dpp4su_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="sglt2_benefit_decile")

obs <- obs %>% mutate(sglt2_benefit_decile=as.factor(sglt2_benefit_decile))
obs_v_pred <- pred %>% inner_join(obs, by="sglt2_benefit_decile")


empty_tick <- obs_v_pred %>%
  filter(sglt2_benefit_decile==1) %>%
  mutate(mean_predicted_benefit=NA, surv_diff=NA, lower_ci=NA, upper_ci=NA, sglt2_benefit_decile=as.factor(0))

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_qdhf_predicted_benefit*100)) + 
  geom_point(aes(y = surv_diff*100)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1) +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  xlab("Predicted SGLT2 benefit (mean by predicted decile)") + ylab("Observed SGLT2 benefit (by predicted decile)") +
  scale_x_continuous(breaks=c(seq(0,6,by=1)), limits=c(0,6)) +
  scale_y_continuous(breaks=c(seq(-2,8,by=1)), limits=c(-2,8)) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5)))

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))

#NNT
obs_v_pred <- obs_v_pred %>%
  mutate(nnt_predicted=1/(mean_qdhf_predicted_benefit*100),
         nnt_observed=1/(surv_diff*100))
obs_v_pred




## Observed adjusted for age, sex, duration, IMD, drugline and ncurrtx

## QRISK2

memory.limit(size=40000)
gc()

noncal_cohort$sglt2_benefit_decile <- ntile(noncal_cohort$qrisk2_sglt2_benefit, 10)
noncal_cohort <- noncal_cohort %>% mutate(sglt2_benefit_decile=as.factor(sglt2_benefit_decile))

# model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + drugline_all + ncurrtx, data=noncal_cohort, x=TRUE, y=TRUE, surv=TRUE)
# 
# obs_SGLT2 <- noncal_cohort %>%
#   select(patid, sglt2_benefit_decile, dstartdate_age, malesex, dstartdate_dm_dur_all, imd2015_10, drugline_all, ncurrtx) %>%
#   mutate(studydrug="SGLT2",
#          rowno=row_number())
# 
# observed_sglt2 <- survfit(model, newdata=as.data.frame(obs_SGLT2)) %>%
#   tidy() %>%
#   filter(time==5) %>%
#   pivot_longer(cols=-c(time, n.risk, n.event, n.censor), names_to = c(".value", "group"), names_pattern = "(.*)\\.(.*)") %>%
#   select(group, estimate_sglt2=estimate) %>%
#   mutate(group=as.numeric(group)) %>%
#   inner_join(obs_SGLT2, by=c("group"="rowno"))
# 
# save(observed_sglt2, file="../../Raw data/HF_adjust_SGLT2.Rda")
# 
# obs_DPP4SU <- noncal_cohort %>%
#   select(patid, sglt2_benefit_decile, dstartdate_age, malesex, dstartdate_dm_dur_all, imd2015_10, drugline_all, ncurrtx) %>%
#   mutate(studydrug="DPP4SU",
#          rowno=row_number())
# 
# observed_dpp4su <- survfit(model, newdata=as.data.frame(obs_DPP4SU)) %>%
#   tidy() %>%
#   filter(time==5) %>%
#   pivot_longer(cols=-c(time, n.risk, n.event, n.censor), names_to = c(".value", "group"), names_pattern = "(.*)\\.(.*)") %>%
#   select(group, estimate_dpp4su=estimate) %>%
#   mutate(group=as.numeric(group)) %>%
#   inner_join(obs_DPP4SU, by=c("group"="rowno"))
# 
# save(observed_dpp4su, file="../../Raw data/HF_adjust_DPP4SU.Rda")

load("../../Raw data/HF_adjust_SGLT2.Rda")
load("../../Raw data/HF_adjust_DPP4SU.Rda")

observed <- observed_sglt2 %>%
  select(group, estimate_sglt2) %>%
  inner_join(observed_dpp4su, by="group") %>%
  mutate(survdiff=estimate_sglt2-estimate_dpp4su) %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(mean_benefit=mean(survdiff),
            median_benefit=median(survdiff),
            lq_benefit=quantile(survdiff, prob=c(.25)),
            uq_benefit=quantile(survdiff, prob=c(.75)))

### Predicted same as previous
pred <- noncal_cohort %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(median_qrisk2_predicted_benefit=median(qrisk2_sglt2_benefit, na.rm=T))


obs_v_pred <- pred %>% inner_join(observed, by="sglt2_benefit_decile")


empty_tick <- obs_v_pred %>%
  filter(sglt2_benefit_decile==1) %>%
  mutate(mean_benefit=NA, median_benefit=NA, lq_benefit=NA, uq_benefit=NA, mean_qrisk2_predicted_benefit=NA, sglt2_benefit_decile=as.factor(0))

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=median_qrisk2_predicted_benefit*100)) + 
  geom_point(aes(y = median_benefit*100)) +
  geom_errorbar(aes(ymax=uq_benefit*100,ymin=lq_benefit*100),alpha=1) +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  xlab("Predicted SGLT2 benefit (mean by predicted decile)") + ylab("Observed SGLT2 benefit (by predicted decile)") +
  scale_x_continuous(breaks=c(seq(0,6,by=1)), limits=c(0,6)) +
  scale_y_continuous(breaks=c(seq(-2,6,by=1)), limits=c(-2,8)) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5)))

p1





noncal_cohort$sglt2_benefit_decile <- ntile(noncal_cohort$qrisk2_sglt2_benefit, 10)
noncal_cohort <- noncal_cohort %>% mutate(sglt2_benefit_decile=as.factor(sglt2_benefit_decile))

pred <- noncal_cohort %>%
  group_by(sglt2_benefit_decile) %>%
  summarise(mean_qrisk2_predicted_benefit=mean(qrisk2_sglt2_benefit, na.rm=T))


noncal_cohort <- noncal_cohort %>% mutate(studydrug=as.vector(studydrug))
ddist <- datadist(noncal_cohort)

model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + drugline_all + ncurrtx + qrisk2_5yr_score, data=noncal_cohort, x=TRUE, y=TRUE, surv=TRUE)

obs_SGLT2 <- noncal_cohort %>%
    select(patid, sglt2_benefit_decile, dstartdate_age, malesex, dstartdate_dm_dur_all, imd2015_10, drugline_all, ncurrtx, qrisk2_5yr_score) %>%
    mutate(studydrug="SGLT2",
           rowno=row_number())


memory.limit(size=40000)
gc()
survival_est2 <- survest(model, newdata=as.data.frame(obs_SGLT2), times=5)







