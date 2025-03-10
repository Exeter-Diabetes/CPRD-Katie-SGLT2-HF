
# HRs for SGLT2 vs DPP4SU for HF
# Stratified by QDHF
# Now 3 splines as main analysis
 
############################################################################################

# Setup
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(forestplot)
library(rms)
library(cowplot)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# 1 Cohort selection (see script 00)

setwd("/slade/CPRD_data/Katie SGLT2/Processed data")
load("treatment_outcome_cohort_jun24.rda")

# Add survival variables
setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions")
source("survival_variables.R")

cohort <- add_surv_vars(cohort, main_only=T)

# now using functions to define covariates adjusted and weighted for
setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions")
source("full_covariate_set.R")


############################################################################################

# 2 Make spline of HR by HF risk

# Adjusted plot
## Adjusting for same as previous except not QDHF - put in as interaction term with studydrug

cohort <- cohort %>% mutate(studydrug=as.vector(studydrug))

# Remove constant variables otherwise datadist has issues
cohort <- cohort %>%
  select(patid, studydrug, malesex, dstartdate_age, dstartdate_dm_dur_all, ethnicity_qrisk2_decoded, imd2015_10, qrisk2_smoking_cat, hypertension, predrug_af, hosp_admission_prev_year_count, prebmi, prehba1c2yrs, presbp, drugline_all, ncurrtx_cat, INS, initiation_year, qdiabeteshf_5yr_score, hf_censtime_yrs, hf_censvar)

ddist <- datadist(cohort)
options(datadist='ddist')


return_cov_set("adjust")
#"malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score"
#remove qdiabeteshf_5yr_score and add as interaction term


model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug*rcs(qdiabeteshf_5yr_score,3) + malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year, data=cohort, x=T, y=T)

anova(model)
# no interaction - good (p=0.82)

c1 <- quantile(cohort$qdiabeteshf_5yr_score, .01, na.rm=TRUE)
c99 <- quantile(cohort$qdiabeteshf_5yr_score, .99, na.rm=TRUE)


contrast_spline <- contrast(model, list(studydrug = "SGLT2", qdiabeteshf_5yr_score = seq(c1, c99, by=0.05)), list(studydrug = "DPP4SU", qdiabeteshf_5yr_score = seq(c1, c99, by=0.05)))

contrast_spline_df <- as.data.frame(contrast_spline[c('qdiabeteshf_5yr_score','Contrast','Lower','Upper')])



# Add histogram
marginal_distribution <- function(x,var) {
  ggplot(x, aes_string(x = var)) +
    geom_histogram(bins = 64, alpha = 0.4, position = "identity") +
    guides(fill = FALSE) +
    theme(legend.title = element_blank(), panel.background = element_rect( fill = "white",color = "grey50")) +
    scale_x_continuous(breaks = seq(0,15,2.5)) +
    xlab(expression(paste("5-year absolute heart failure risk (%)"))) +
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


spline_plot <- ggplot(data=contrast_spline_df, aes(x=qdiabeteshf_5yr_score, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=qdiabeteshf_5yr_score, y=exp(Contrast)), size=1) +
  xlab(expression(paste("5-year absolute heart failure risk (%)"))) +
  ylab("Heart failure HR: SGLT2i vs comparator") +
  coord_trans(y = "log10") +
  scale_x_continuous(breaks = seq(0,15,2.5)) +
  geom_ribbon(data=contrast_spline_df, aes(x=qdiabeteshf_5yr_score, ymin=exp(Lower), ymax=exp(Upper)), alpha=0.5) +
  
  geom_hline(yintercept = 1, linetype = "dashed")  +
  
  geom_hline(yintercept = 0.63, linetype = "dashed", size=0.7, color="red")  +
  geom_hline(yintercept = 0.5, linetype = "dashed", size=0.7, color="red")  +
  geom_hline(yintercept = 0.8, linetype = "dashed", size=0.7, color="red")  +
  
  geom_segment(aes(x = 1.9, xend = 2.8, y = 1.6, yend = 1.6), linetype = "dashed", linewidth=0.7, color="red") +
  
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
  annotate(geom="text", x=3, y=2, label="Overall HR in study cohort (95% CI): 0.70 (0.63 to 0.78)", color="black", hjust=0, size = 6) +
  annotate(geom="text", x=3, y=1.8, label="p-value for HF risk * drug arm interaction on HF outcome: 0.82", color="black", hjust=0, size = 6) +
  annotate(geom="text", x=3, y=1.6, label="Trial meta-analysis HR (95% CI): 0.63 (0.50 to 0.80)", color="red", hjust=0, size = 6) 




hist.dta <- cohort %>% filter(qdiabeteshf_5yr_score>=c1 &  qdiabeteshf_5yr_score <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "qdiabeteshf_5yr_score")



setwd("/slade/CPRD_data/Katie SGLT2/Plots/")
tiff("HR_by_HF_risk.tiff", width=9, height=7, units = "in", res=800) 

plot_grid(spline_plot, x_hist, ncol = 1,align = 'v',
          rel_heights = c(1,0.5), rel_widths = c(1,1))

dev.off()






# Testing different numbers of splines

model_3 <- cph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug*rcs(qdiabeteshf_5yr_score,3) + malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year, data=cohort, x=T, y=T)

anova(model_3)
#p=0.82 as above

model_4 <- cph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug*rcs(qdiabeteshf_5yr_score,4) + malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year, data=cohort, x=T, y=T)

anova(model_4)
#p=0.93

model_5 <- cph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug*rcs(qdiabeteshf_5yr_score,5) + malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year, data=cohort, x=T, y=T)

anova(model_5)
#p=0.71



lrtest(model_4, model_3)
# L.R. Chisq       d.f.          P 
# 5.64705483 2.00000000 0.05939606 

lrtest(model_5, model_4)
# L.R. Chisq       d.f.          P 
# 8.17178651 2.00000000 0.01680812 



model_3 <- coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug*rcs(qdiabeteshf_5yr_score,3) + malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year, data=cohort, x=T, y=T)

model_4 <- coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug*rcs(qdiabeteshf_5yr_score,4) + malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year, data=cohort, x=T, y=T)

model_5 <- coxph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug*rcs(qdiabeteshf_5yr_score,5) + malesex + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year, data=cohort, x=T, y=T)


stats::anova(model_3, model_4)
# loglik  Chisq Df Pr(>|Chi|)  
# 1 -29966                       
# 2 -29963 5.6473  2    0.05939 .
#same as above

stats::anova(model_4, model_5)
# loglik  Chisq Df Pr(>|Chi|)  
# 1 -29963                       
# 2 -29959 8.1717  2    0.01681 *
#same as above




############################################################################################

# Stratified by sex

males <- cohort %>% filter(malesex==1) %>% select(-malesex)

ddist <- datadist(males)
options(datadist='ddist')

# sex removed
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug*rcs(qdiabeteshf_5yr_score,3) + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year, data=males, x=T, y=T)

anova(model)
# no interaction - good (p=0.80)



females <- cohort %>% filter(malesex==0) %>% select(-malesex)

ddist <- datadist(females)
options(datadist='ddist')

# sex removed
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~  studydrug*rcs(qdiabeteshf_5yr_score,3) + rcs(dstartdate_age,5) + rcs(dstartdate_dm_dur_all,5) + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + rcs(prebmi,5) + rcs(prehba1c2yrs,5) + rcs(presbp,5) + drugline_all + ncurrtx_cat + INS + initiation_year, data=females, x=T, y=T)

anova(model)
# no interaction - good (p=0.17)
