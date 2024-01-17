
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
load("20231121_t2d_1stinstance.Rda")
load("20231121_t2d_all_drug_periods.Rda")

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/Scripts/Functions")
source("cohort_definition.R")

cohort <- define_cohort(t2d_1stinstance, t2d_all_drug_periods)

cohort <- cohort %>% filter(drugclass!="GLP1") %>% mutate(studydrug=factor(ifelse(drugclass=="SGLT2", "SGLT2", "DPP4SU"), levels=c("DPP4SU", "SGLT2")))
#cohort <- cohort %>% filter(drugclass!="GLP1" & drugclass!="DPP4") %>% mutate(studydrug=factor(drugclass, levels=c("SU", "SGLT2")))
#cohort <- cohort %>% filter(drugclass!="GLP1" & drugclass!="SGLT2") %>% mutate(studydrug=factor(drugclass, levels=c("SU", "DPP4")))

table(cohort$studydrug)
# DPP4SU 98397 
# SGLT2 45781

# SU 36378
# SGLT2 45781

# SU 36378
# DPP4 58964 


## B Make variables for survival analysis of all endpoints (see survival_variables function for details)

source("survival_variables.R")

cohort <- add_surv_vars(cohort, main_only=FALSE)


## C Just keep variables of interest

cohort <- cohort %>%
  
  select(patid, studydrug, drugclass, drugsubstances, dstartdate, dstopdate, dstartdate_age, malesex, ethnicity_qrisk2_decoded, imd2015_10, qrisk2_smoking_cat, drugline_all, ncurrtx, initiation_year, drugorder, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_dm_dur_all, prehba1c, prebmi, pretotalcholesterol, prehdl, preegfr, preckdstage, prealt, presbp, contains("last_6_months"), hosp_admission_prev_year, hosp_admission_prev_year_count, predrug_cld, predrug_hypertension, predrug_af, predrug_copd, predrug_haem_cancer, predrug_solid_cancer, predrug_diabeticnephropathy, predrug_neuropathy, predrug_retinopathy, predrug_dementia, predrug_otherneuroconditions, qdiabeteshf_5yr_score, qdiabeteshf_lin_predictor, contains("cens"), regstartdate, gp_record_end, death_date)

rm(list=setdiff(ls(), "cohort"))

############################################################################################


#Define deciles of predicted probabilities
cohort <- cohort %>% mutate(qdhf_decile=ntile(qdiabeteshf_5yr_score, 10))

# Get mean predicted probabilities from QDHF for each QDHF cat by studydrug
predicted <- cohort %>%
  filter(studydrug=="DPP4SU") %>%
  group_by(qdhf_decile) %>%
  summarise(mean_qdhf_pred=mean(qdiabeteshf_5yr_score)/100)


# Find actual observed probabilities by QRISK2 category and studydrug
observed <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ qdhf_decile, data=cohort[cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed=1-estimate,
         lower_ci=1-conf.low,
         upper_ci=1-conf.high) %>%
  select(observed, lower_ci, upper_ci, strata)

events <- cohort %>%
  filter(studydrug=="DPP4SU" & hf_censvar==1) %>%
  group_by(qdhf_decile) %>%
  summarise(events=n())

obs_v_pred <- cbind(predicted, observed)

events_table <- data.frame(t(events)) %>%
  rownames_to_column() %>%
  filter(rowname!="qdhf_decile")


dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(qdhf_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_qdhf_pred=NA, qdhf_decile=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=qdhf_decile)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_qdhf_pred*100), position=dodge, shape=17, size=2) +
  theme_bw() +
  xlab("QDiabetes-Heart Failure decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,25,by=5)), limits=c(0,27)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Uncalibrated QDHF vs HF incidence (5 year)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))


# QDHF looks good


## D Look at C-stats

cohort <- cohort %>%
  mutate(qdiabeteshf_survival=(100-qdiabeteshf_5yr_score)/100)
  
### QDHF
surv_mod <- coxph(Surv(hf_censtime_yrs, hf_censvar)~qdiabeteshf_survival,data=cohort,method="breslow")
summary(surv_mod)$concordance[1]
summary(surv_mod)$concordance[1]-(1.96*summary(surv_mod)$concordance[2])
summary(surv_mod)$concordance[1]+(1.96*summary(surv_mod)$concordance[2])
# 0.7244267 (0.7147242-0.7341292)

#Brier score
library()





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
  mutate(centred_qdhf_lin_predictor=qdiabeteshf_lin_predictor-mean(qdiabeteshf_lin_predictor)) %>%
  ungroup() %>%
  mutate(qdiabeteshf_5yr_score_cal=ifelse(malesex==1, (1-(male_qdhf_surv^exp(centred_qdhf_lin_predictor)))*100, (1-(female_qdhf_surv^exp(centred_qdhf_lin_predictor)))*100))


# Plot calibrated + uncalibrated + observed scores

# Get mean predicted probabilities for calibrated and uncalibrated
predicted <- noncal_cohort %>%
  filter(studydrug=="DPP4SU") %>%
  group_by(qdhf_decile) %>%
  summarise(mean_qdhf_pred=mean(qdiabeteshf_5yr_score)/100,
            mean_qdhf_pred_cal=mean(qdiabeteshf_5yr_score_cal)/100)


# Find actual observed probabilities by QRISK2 category and studydrug
observed <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ qdhf_decile, data=noncal_cohort[noncal_cohort$studydrug=="DPP4SU",]) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed=1-estimate,
         lower_ci=1-conf.low,
         upper_ci=1-conf.high) %>%
  select(observed, lower_ci, upper_ci, strata)

events <- cohort %>%
  filter(studydrug=="DPP4SU" & hf_censvar==1) %>%
  group_by(qdhf_decile) %>%
  summarise(events=n())


obs_v_pred <- cbind(predicted, observed) 
events_table <- data.frame(t(events)) %>%
  rownames_to_column() %>%
  filter(rowname!="qdhf_decile")


dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(qdhf_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_qdhf_pred=NA, mean_qdhf_pred_cal=NA, qdhf_decile=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=qdhf_decile)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_qdhf_pred_cal*100), position=dodge, shape=17, size=2) +
  theme_bw() +
  xlab("QDiabetes-Heart Failure decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,25,by=5)), limits=c(0,27)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Calibrated QDHF vs HF incidence (5 year)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))



############################################################################################

# Add HR from trials
cohort <- noncal_cohort

cohort <- cohort %>%
  mutate(qdhf_survival_cal=(100-qdiabeteshf_5yr_score_cal)/100,
         qdhf_survival_cal_sglt2=qdhf_survival_cal^0.63,
         qdiabeteshf_5yr_score_cal_sglt2=100-(qdhf_survival_cal_sglt2*100))


# Distribution of predicted benefits
cohort <- cohort %>%
  mutate(qdhf_sglt2_benefit=qdhf_survival_cal_sglt2-qdhf_survival_cal)

summary(cohort$qdhf_sglt2_benefit)
#mean=1.7%


# Cut into quartiles of predicted benefit
cohort$sglt2.benefit.quartile <- ntile(cohort$qdhf_sglt2_benefit, 4)

## Average benefit per quartile
cohort %>%
  group_by(sglt2.benefit.quartile) %>%
  summarise(pred.mean=mean(qdhf_sglt2_benefit,na.rm=T),
            pred.median=median(qdhf_sglt2_benefit,na.rm=T),
            sd=sd(qdhf_sglt2_benefit,na.rm=T),
            l_iqr=quantile(qdhf_sglt2_benefit,na.rm=T,probs=0.25),
            u_iqr=quantile(qdhf_sglt2_benefit,na.rm=T,probs=0.75),
            min=min(qdhf_sglt2_benefit,na.rm=T),
            max=max(qdhf_sglt2_benefit,na.rm=T),
            NNT=1/(mean(qdhf_sglt2_benefit,na.rm=T)))



## Purple and orange plots

pdf.options(reset = TRUE, onefile = FALSE)


#sfit <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug,data = cohort,subset=sglt2.benefit.quartile==1)
#pdf("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/quartile1.pdf", width=4, height=5)
#sfit <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug,data = cohort,subset=sglt2.benefit.quartile==2)
#pdf("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/quartile2.pdf", width=4, height=5)
#sfit <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug,data = cohort,subset=sglt2.benefit.quartile==3)
#pdf("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/quartile3.pdf", width=4, height=5)
sfit <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug,data = cohort,subset=sglt2.benefit.quartile==4)
pdf("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2023/1 SGLT2 CVD project/quartile4.pdf", width=4, height=5) 

plot1 <- ggsurvplot(sfit,
                    fun = function(x) {100 - x*100},             # Cumulative probability plot
                    censor = F,
                    size= 1.5,
                    #legend.title = "" ,
                    #legend.labs = c("DPP4/SU", "SGLT2"),
                    #legend.title = "",
                    risk.table = F,       # show risk table.
                    conf.int = T,         # show confidence intervals for
                    # point estimates of survival curves.
                    #xlim = c(0,0.06),         # present narrower X axis, but not affect
                    # survival estimates.
                    ylim = c(0,12),
                    xlab = "Years",   # customize X axis label.
                    ylab = "Heart failure incidence (%)",   # customize X axis label.
                    break.time.by = 1,     # break X axis in time intervals every two years.
                    ggtheme = theme_classic(), # customize plot and risk table with a theme.
                    risk.table.y.text.col = TRUE, # colour risk table text annotations.
                    #fontsize = 3, # in legend of risk table
                    font.x = c(18),font.y = c(18),font.tickslab = c(18),
                    axes.offset = FALSE, # start the axis at the origin
                    palette = c("#998ec3","#f1a340"),
                    #linetype = c("strata"),
                    legend = "none",
                    #title = "Predicted 5yr SGLT2 benefit 0.2%",
                    #tables.theme = theme_cleantable() # clean theme for tables
                    
)
plot1$plot + theme(plot.margin = margin(20,20,20,20))
dev.off()



### NNT

## Need SGLT2 benefit for each quartile; NNT=100/this

### Quartile 1
x <- summary(survfit(Surv(hf_censtime_yrs,hf_censvar) ~ studydrug, data=cohort, subset=sglt2.benefit.quartile==1), times=5)
## Difference between survival
(x$surv[2]-x$surv[1])*100
#0.14%
(x$lower[2]-x$lower[1])*100
#-0.017%
(x$upper[2]-x$upper[1])*100
#0.30%

#NNT=100/0.14=>700

NNT=100/0.1
=>700

100/0.1


### Quartile 2
x <- summary(survfit(Surv(hf_censtime_yrs,hf_censvar) ~ studydrug, data=cohort, subset=sglt2.benefit.quartile==2), times=5)
## Difference between survival
(x$surv[2]-x$surv[1])*100
#0.33%
(x$lower[2]-x$lower[1])*100
#0.178%
(x$upper[2]-x$upper[1])*100
#0.488%

#NNT=100/0.33=300


### Quartile 3
x <- summary(survfit(Surv(hf_censtime_yrs,hf_censvar) ~ studydrug, data=cohort, subset=sglt2.benefit.quartile==3), times=5)
## Difference between survival
(x$surv[2]-x$surv[1])*100
#0.858%
(x$lower[2]-x$lower[1])*100
#0.737%
(x$upper[2]-x$upper[1])*100
#0.979%

#NNT=100/0.858=116


### Quartile 4
x <- summary(survfit(Surv(hf_censtime_yrs,hf_censvar) ~ studydrug, data=cohort, subset=sglt2.benefit.quartile==4), times=5)
## Difference between survival
(x$surv[2]-x$surv[1])*100
#2.68%
(x$lower[2]-x$lower[1])*100
#2.35%
(x$upper[2]-x$upper[1])*100
#3.02%

#NNT=100/2.68=37





model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug*sglt2.benefit.quartile, data=cohort, x=TRUE, y=TRUE, surv=TRUE)

survival_est <- survest(model, newdata=expand.grid(studydrug=c("SGLT2","DPP4SU"), sglt2.benefit.quartile=c(1:4)), times=5)

obs <- data.frame(surv=unlist(survival_est$surv), studydrug=rep(c("SGLT2","DPP4SU"),4), sglt2.benefit.quartile=rep(1:4, rep_len(2, 4)))

obs <- obs %>%
  pivot_wider(id_cols=sglt2.benefit.quartile, values_from=surv, names_from=studydrug) %>%
  mutate(surv_diff=(SGLT2-DPP4SU)*100) %>%
  select(sglt2.benefit.quartile, surv_diff)




# quartile characteristics
library(gtsummary)

cohort %>%
  mutate(imd2015_10=as.factor(imd2015_10),
         age_70_or_over=ifelse(dstartdate_age>=70, 1, 0),
         smoker=ifelse(qrisk2_smoking_cat<2, qrisk2_smoking_cat, 2)) %>%
  select(sglt2.benefit.quartile, malesex, ethnicity_5cat_decoded, imd2015_10, dstartdate_age, age_70_or_over, dstartdate_dm_dur_all, prehba1c, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preegfr, preckdstage, smoker) %>%
  tbl_summary(by=sglt2.benefit.quartile, missing="no")


cohort %>%
  mutate(imd2015_10=as.factor(imd2015_10),
         age_70_or_over=ifelse(dstartdate_age>=70, 1, 0),
         smoker=ifelse(qrisk2_smoking_cat<2, qrisk2_smoking_cat, 2)) %>%
  select(sglt2.benefit.quartile, malesex, ethnicity_5cat_decoded, imd2015_10, dstartdate_age, age_70_or_over, dstartdate_dm_dur_all, prehba1c, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preegfr, preckdstage, smoker) %>%
  tbl_summary(missing="no")



