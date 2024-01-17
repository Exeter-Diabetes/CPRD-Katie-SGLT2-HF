
## 1 Setup and calibration

## 2 Is QDiabetes constant by HF risk?

## 3 Distribution of predicted risk

## 4 Quartiles


############################################################################################

# Setup
library(tidyverse)
library(survival)
library(rms)
library(cowplot)
library(survminer)
#library(broom)
#library(patchwork)



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
  
  select(patid, malesex, ethnicity_5cat_decoded, imd2015_10, regstartdate, gp_record_end, death_date, drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, prehba1c, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preegfr, preckdstage, qrisk2_smoking_cat, contains("cens"), qrisk2_lin_predictor, qrisk2_5yr_score, qdiabeteshf_lin_predictor, qdiabeteshf_5yr_score, starts_with("ckdpc"), last_sglt2_stop)

rm(list=setdiff(ls(), "cohort"))


## D Calibrate score

full_cohort <- cohort

source("calibrate_risk_scores.R")

cohort <- calibrate_risk_score(cohort, risk_score="qdiabeteshf", outcome="hf")

table(cohort$studydrug)
# DPP4SU 95362
# SGLT2 50234


############################################################################################

# 2 Is HR constant with baseline QDiabetes-HF?

ddist <- datadist(cohort); options(datadist='ddist')


# Adjusted plot
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug*rcs(qdiabeteshf_5yr_score_cal,5) + dstartdate_age  + malesex + dstartdate_dm_dur_all + imd2015_10 + drugline_all + ncurrtx,
                      data = cohort,x=T,y=T)
anova(model)
#no interaction

c1 <- quantile(cohort$qdiabeteshf_5yr_score_cal, .01, na.rm=TRUE)
c99 <- quantile(cohort$qdiabeteshf_5yr_score_cal, .99, na.rm=TRUE)

contrast_spline.1 <- contrast(model, list(studydrug="SGLT2", qdiabeteshf_5yr_score_cal=seq(c1,c99,by=0.05)), list(studydrug="DPP4SU", qdiabeteshf_5yr_score_cal=seq(c1, c99, by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('qdiabeteshf_5yr_score_cal','Contrast','Lower','Upper')])

#plot and save
contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=qdiabeteshf_5yr_score_cal, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=qdiabeteshf_5yr_score_cal, y=exp(Contrast)), size=1) +
  xlab(expression(paste("QDiabetes-Heart Failure 5 year score"))) +
  ylab("Hazard Ratio") +
  scale_x_continuous(breaks = seq(0,25,5)) +
  geom_ribbon(data=contrast_spline_df,aes(x=qdiabeteshf_5yr_score_cal,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  geom_hline(yintercept = 0.63, color="red", size=0.5)  +
  geom_hline(yintercept = 0.5, linetype = "twodash", color="red", size=0.5)  +
  geom_hline(yintercept = 0.8, linetype = "twodash", color="red", size=0.5)  +
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
    geom_histogram(bins = 100, alpha = 0.4, position = "identity") +
    guides(fill = FALSE) +
    theme(legend.title = element_blank(), panel.background = element_rect( fill = "white",color = "grey50")) +
    scale_x_continuous(breaks = seq(0,50,5)) +
    scale_y_continuous(breaks = seq(0,5000,2000)) + #- for small graph only
    xlab(expression(paste("QDiabetes-Heart Failure 5 year score"))) +
    ylab(expression(paste("Frequency"))) +
    theme(text = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(color="grey50"),
          axis.line.y = element_line(color="grey50"))
  
}


hist.dta <- cohort %>% filter(qdiabeteshf_5yr_score_cal>=c1 &  qdiabeteshf_5yr_score_cal <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "qdiabeteshf_5yr_score_cal")


# Arranging the plot using cowplot
plot_grid(x_hist, contrast_spline_plot_1, ncol = 1,align = 'hv',
          rel_heights = c(1,1), rel_widths = c(1,1))
#export as portrait pdf 7 x 6

# Arranging the plot using cowplot
plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'v',
          rel_heights = c(1,0.3), rel_widths = c(1,1))
#export as landscape pdf 7 x 5.5


# Hazard ratio by QDHF

cohort <- cohort %>%
  mutate(qdhf_cat=ifelse(qdiabeteshf_5yr_score_cal<5, 1,
                         ifelse(qdiabeteshf_5yr_score_cal<10, 2,
                                ifelse(qdiabeteshf_5yr_score_cal<15, 3, 4))))

table(cohort$qdhf_cat)

summary(coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qdiabeteshf_5yr_score_cal + drugline_all + ncurrtx, data=cohort, subset=qdhf_cat==1))

summary(coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qdiabeteshf_5yr_score_cal + drugline_all + ncurrtx, data=cohort, subset=qdhf_cat==2))

summary(coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qdiabeteshf_5yr_score_cal + drugline_all + ncurrtx, data=cohort, subset=qdhf_cat==3))

summary(coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qdiabeteshf_5yr_score_cal + drugline_all + ncurrtx, data=cohort, subset=qdhf_cat==4))

summary(coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qdiabeteshf_5yr_score_cal + drugline_all + ncurrtx, data=cohort))


## Forest plot
library(broom)
library(forestplot)

trial <- data.frame(estimate=0.63, conf.low=0.5, conf.high=0.8, cat="trial")

whole <- coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qdiabeteshf_5yr_score_cal + drugline_all + ncurrtx, data=cohort) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(cat="whole")


qdhf1 <- coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qdiabeteshf_5yr_score_cal + drugline_all + ncurrtx, data=cohort, subset=qdhf_cat==1) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(cat="qdhf==1")

qdhf2 <- coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qdiabeteshf_5yr_score_cal + drugline_all + ncurrtx, data=cohort, subset=qdhf_cat==2) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(cat="qdhf==2")

qdhf3 <- coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qdiabeteshf_5yr_score_cal + drugline_all + ncurrtx, data=cohort, subset=qdhf_cat==3) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(cat="qdhf==3")

qdhf4 <- coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qdiabeteshf_5yr_score_cal + drugline_all + ncurrtx, data=cohort, subset=qdhf_cat==4) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(cat="qdhf==4")

data <- rbind (trial, whole, empty=data.frame(estimate=NA, conf.low=NA, conf.high=NA, cat=NA), qdhf1, qdhf2, qdhf3, qdhf4)




## With interaction term

cohort <- cohort %>% mutate(qdhf_cat=as.factor(qdhf_cat))

cohort$qdhf_cat <- relevel(cohort$qdhf_cat, ref=1)

qdhf1 <- coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug*qdhf_cat + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qdiabeteshf_5yr_score_cal + drugline_all + ncurrtx, data=cohort) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(cat="qdhf==1")

cohort$qdhf_cat <- relevel(cohort$qdhf_cat, ref=2)

qdhf2 <- coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug*qdhf_cat + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qdiabeteshf_5yr_score_cal + drugline_all + ncurrtx, data=cohort) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(cat="qdhf==2")

cohort$qdhf_cat <- relevel(cohort$qdhf_cat, ref=3)

qdhf3 <- coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug*qdhf_cat + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qdiabeteshf_5yr_score_cal + drugline_all + ncurrtx, data=cohort) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(cat="qdhf==3")

cohort$qdhf_cat <- relevel(cohort$qdhf_cat, ref=4)

qdhf4 <- coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug*qdhf_cat + dstartdate_age + malesex + dstartdate_dm_dur_all + imd2015_10 + qdiabeteshf_5yr_score_cal + drugline_all + ncurrtx, data=cohort) %>%
  tidy(conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term=="studydrugSGLT2") %>%
  select(estimate, conf.low, conf.high) %>%
  mutate(cat="qdhf==4")


data <- rbind(trial, whole, empty=data.frame(estimate=NA, conf.low=NA, conf.high=NA, cat=NA), qdhf1, qdhf2, qdhf3, qdhf4)





styles <- fpShapesGp(
  lines = list(
    gpar(col = "#C60046"),
    gpar(col = "black"),
    gpar(col = "white"),
    gpar(col = "blue"),
    gpar(col = "blue"),
    gpar(col = "blue"),
    gpar(col = "blue")
  ),
  box = list(
    gpar(fill = "#C60046", col="#C60046"),
    gpar(fill = "black", col="black"),
    gpar(fill = "white", col="white"),
    gpar(fill = "blue", col="blue"),
    gpar(fill = "blue", col="blue"),
    gpar(fill = "blue", col="blue"),
    gpar(fill = "blue", col="blue")
  ) 
)


forestplot(data$cat,
           mean = data$estimate,
           lower= data$conf.low,
           upper = data$conf.high,
           xticks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
           zero = 0,
           ci.vertices=TRUE,
           col=fpColors(box="blue", lines="blue", zero = "gray50"),
           #xlab="HbA1c response difference per SD change (mmol/mol)",cex=1,
           lwd.ci=4,
           new_page = TRUE,
           fn.ci_norm = fpDrawCircleCI,
           shapes_gp = styles,
           boxsize = .2,
           txt_gp = fpTxtGp(legend  = gpar(cex = 2.5),xlab  = gpar(cex = 2.5),ticks  = gpar(cex = 2.5)))








############################################################################################

# 3 Distribution of predicted risk is in histogram above

############################################################################################

# 4 Quartiles

# Add hazard ratio from trials meta-analysis
# https://jamanetwork.com/journals/jamacardiology/fullarticle/2771459

# Add HR from trials
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

sfit <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug,data = cohort,subset=sglt2.benefit.quartile==1)
pdf("../../../DUK/Presentation and Poster/quartile1.pdf", width=4, height=5)
sfit <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug,data = cohort,subset=sglt2.benefit.quartile==2)
pdf("../../../DUK/Presentation and Poster/quartile2.pdf", width=4, height=5)
sfit <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug,data = cohort,subset=sglt2.benefit.quartile==3)
pdf("../../../DUK/Presentation and Poster/quartile3.pdf", width=4, height=5)
sfit <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug,data = cohort,subset=sglt2.benefit.quartile==4)
pdf("../../../DUK/Presentation and Poster/quartile4.pdf", width=4, height=5) 

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








########################################################################################################################

#calibration plots

library(broom)
library(patchwork)

cohort$qdhf_decile <- ntile(cohort$qdiabeteshf_5yr_score, 10)

predicted <- cohort %>%
  group_by(qdhf_decile) %>%
  summarise(mean_pred=mean(qdiabeteshf_5yr_score)/100)


observed <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ qdhf_decile, data=cohort) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed=1-estimate,
         lower_ci=1-conf.low,
         upper_ci=1-conf.high) %>%
  select(observed, lower_ci, upper_ci, strata)

events <- cohort %>%
  filter(mace_censvar==1) %>%
  group_by(studydrug, qdhf_decile) %>%
  summarise(overall_events=n())


obs_v_pred <- cbind(predicted, observed)

events_table <- data.frame(t(events %>%
  pivot_wider(id_cols=qdhf_decile, names_from=studydrug, values_from=overall_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="qdhf_decile")


dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(qdhf_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_pred=NA, qdhf_decile=0)

p1 <- ggplot(data=bind_rows(empty_tick, obs_v_pred), aes(x=qdhf_decile)) +
  geom_point(aes(y = mean_pred*100), position=dodge, color="red", size=2) +
  geom_point(aes(y = observed*100), position=dodge, size=2) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1,width=0.25, size=0.8, position=dodge) +
  theme_bw() +
  xlab("QDiabetes-HF decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,20,by=5)), limits=c(0,26)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Uncalibrated QDiabetes-HF vs HF (5 year)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))






cohort$qdhf_decile <- ntile(cohort$qdiabeteshf_5yr_score_cal, 10)

predicted <- cohort %>%
  group_by(qdhf_decile) %>%
  summarise(mean_pred=mean(qdiabeteshf_5yr_score_cal)/100)


observed <- survfit(Surv(hf_censtime_yrs, hf_censvar) ~ qdhf_decile, data=cohort) %>%
  tidy() %>%
  group_by(strata) %>%
  filter(time==max(time)) %>%
  mutate(observed=1-estimate,
         lower_ci=1-conf.low,
         upper_ci=1-conf.high) %>%
  select(observed, lower_ci, upper_ci, strata)

events <- cohort %>%
  filter(mace_censvar==1) %>%
  group_by(studydrug, qdhf_decile) %>%
  summarise(overall_events=n())


obs_v_pred <- cbind(predicted, observed)

events_table <- data.frame(t(events %>%
                               pivot_wider(id_cols=qdhf_decile, names_from=studydrug, values_from=overall_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="qdhf_decile")


dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(qdhf_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_pred=NA, qdhf_decile=0)

p1 <- ggplot(data=bind_rows(empty_tick, obs_v_pred), aes(x=qdhf_decile)) +
  geom_point(aes(y = mean_pred*100), position=dodge, color="red", size=2) +
  geom_point(aes(y = observed*100), position=dodge, size=2) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1,width=0.25, size=0.8, position=dodge) +
  theme_bw() +
  xlab("QDiabetes-HF decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,20,by=5)), limits=c(0,26)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Calibrated QDiabetes-HF vs HF (5 year)")

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))
