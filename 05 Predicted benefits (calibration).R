
# Predict benefits

# Plot calibration plot
 
############################################################################################

# Setup
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(rms)
library(ggthemes)
library(PSweight)
library(grid)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# 1 Cohort selection (see script 00)

setwd("/slade/CPRD_data/Katie SGLT2/Processed data/")
load("treatment_outcome_cohort_jun24.rda")

# Add survival variables
setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions/")
source("survival_variables.R")

cohort <- add_surv_vars(cohort, main_only=T)

# now using functions to define covariates adjusted and weighted for
setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions/")
source("full_covariate_set.R")


############################################################################################

# 2 Calculate SGLT2 benefits

cohort <- cohort %>%
  mutate(qdhf_survival=(100-qdiabeteshf_5yr_score)/100,
         qdhf_survival_sglt2=qdhf_survival^0.63,
         qdiabeteshf_5yr_score_sglt2=100-(qdhf_survival_sglt2*100),
         qdhf_sglt2_benefit=qdhf_survival_sglt2-qdhf_survival)


############################################################################################

# 3 Histogram of predicted benefits (not truncated)

summary(cohort$qdhf_sglt2_benefit*100)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.02814  0.55854  1.00996  1.36487  1.79742 14.09242 

# setwd("/slade/CPRD_data/Katie SGLT2/Plots/")
# tiff("histogram_benefits.tiff", width=9, height=5, units = "in", res=600)
# 
# ggplot(cohort, aes(x=qdhf_sglt2_benefit*100)) +
#   geom_histogram(aes(y = after_stat(count / sum(count))*100), binwidth=0.1, alpha=0.5, position="identity") +
#   ylab("Proportion of study population (%)") + xlab("Predicted 5-year SGLT2i benefit (%)") +
#   scale_x_continuous(limits=c(0,15),breaks=c(seq(0,15,by=1))) +
#   scale_y_continuous(limits=c(0,7),breaks=c(seq(0, 7, by=1))) +
#   theme_base() +
#   theme(plot.background = element_blank(),
#         panel.border=element_blank())
# 
# dev.off()



############################################################################################

# 4 Compare to adjusted observed by decile

# Cut into decile of predicted benefit
cohort$predicted_benefit_decile <- as.factor(ntile(cohort$qdhf_sglt2_benefit,10))

## Average benefit per decile
predicted <- cohort %>%
  group_by(predicted_benefit_decile) %>%
  summarise(pred.mean=mean(qdhf_sglt2_benefit,na.rm=T),
            pred.median=median(qdhf_sglt2_benefit,na.rm=T),
            sd=sd(qdhf_sglt2_benefit,na.rm=T),
            l_iqr=quantile(qdhf_sglt2_benefit,na.rm=T,probs=0.25),
            u_iqr=quantile(qdhf_sglt2_benefit,na.rm=T,probs=0.75),
            min=min(qdhf_sglt2_benefit,na.rm=T),
            max=max(qdhf_sglt2_benefit,na.rm=T))




## Observed - adjusted and weighted with overlap weights - using this now

ps.formula <- formula(paste0("studydrug ~ ", return_cov_set("weight")))

overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort), weight="overlap")
cohort$overlap_weights <- overlap$ps.weights$overlap

# Remove constant variables otherwise datadist has issues
cohort <- cohort %>% select(patid, studydrug, qdiabeteshf_5yr_score, dstartdate_age, ethnicity_qrisk2_decoded, malesex, dstartdate_dm_dur_all, imd2015_10, drugline_all, ncurrtx_cat, INS, initiation_year, prebmi, prehba1c2yrs, presbp, qrisk2_smoking_cat, hypertension, hosp_admission_prev_year_count, hf_censtime_yrs, hf_censvar, overlap_weights, qdhf_sglt2_benefit, predicted_benefit_decile, predrug_af)

ddist <- datadist(cohort)
options(datadist='ddist')

f_adjusted <- as.formula(paste0("Surv(hf_censtime_yrs, hf_censvar) ~ studydrug  +  ", return_cov_set("adjust")))

# This is without studydrug * risk interaction term - which would allow HR to vary with risk
model <- cph(f_adjusted, data=cohort, weights=overlap_weights, x=T, y=T, surv=TRUE)



cohort.t <- data.frame(cohort) %>% mutate(studydrug="DPP4SU")
DPP4SU.pred = survest(model,newdata=cohort.t,times=5,se.fit=F)
cohort.t <- data.frame(cohort) %>% mutate(studydrug="SGLT2")
SGLT2.pred = survest(model,newdata=cohort.t,times=5,se.fit=F)

cohort$DPP4SU.pred <- DPP4SU.pred$surv
cohort$SGLT2.pred <- SGLT2.pred$surv

obs <- cohort %>%
  mutate(obs=(SGLT2.pred-DPP4SU.pred)*100) %>%
  group_by(predicted_benefit_decile) %>%
  dplyr::summarise(
    obs.gp = median(obs,na.rm=T),
    l_iqr=quantile(obs,0.25,na.rm=T),
    u_iqr=quantile(obs,0.75,na.rm=T)
  )


predicted.p <- predicted %>% 
  mutate(pred=pred.median*100) %>%
  select(predicted_benefit_decile,pred)
plot.data <- merge(obs,predicted.p,by="predicted_benefit_decile") 


# ## Observed - adjusted and weighted with overlap weights, with studydrug*decile interaction - not using
# 
# ps.formula <- formula(paste0("studydrug ~ ", return_cov_set("weight")))
# 
# overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort), weight="overlap")
# cohort$overlap_weights <- overlap$ps.weights$overlap
# 
# # Remove constant variables otherwise datadist has issues
# cohort <- cohort %>% select(patid, studydrug, qdiabeteshf_5yr_score, dstartdate_age, ethnicity_qrisk2_decoded, malesex, dstartdate_dm_dur_all, imd2015_10, drugline_all, ncurrtx_cat, INS, initiation_year, prebmi, prehba1c2yrs, presbp, qrisk2_smoking_cat, hypertension, hosp_admission_prev_year_count, hf_censtime_yrs, hf_censvar, overlap_weights, predicted_benefit_decile, predrug_af)
# 
# ddist <- datadist(cohort)
# options(datadist='ddist')
# 
# f_adjusted <- as.formula(paste0("Surv(hf_censtime_yrs, hf_censvar) ~ studydrug*predicted_benefit_decile  +  ", return_cov_set("adjust")))
# 
# model <- cph(f_adjusted, data=cohort, weights=overlap_weights, x=T, y=T, surv=TRUE)
# 
# cohort.t <- data.frame(cohort) %>% mutate(studydrug="DPP4SU")
# DPP4SU.pred = survest(model,newdata=cohort.t,times=5,se.fit=F)
# cohort.t <- data.frame(cohort) %>% mutate(studydrug="SGLT2")
# SGLT2.pred = survest(model,newdata=cohort.t,times=5,se.fit=F)


# # from Thijs code
# SGLT2.pred2 = survfit(model,newdata=cohort.t) %>%
#   tidy() %>%
#   filter(time==5) %>%
#   pivot_longer(cols=-c(time, n.risk, n.event, n.censor), names_to = c(".value", "group"), names_pattern = "(.*)\\.(.*)") %>%
#   select(group, estimate_sglt2=estimate)
# 
# 
# test <- cbind(SGLT2.pred2, SGLT2.pred$surv)
# #gives exactly the same answer!!



# Overall estimates

summary((cohort %>% mutate(obs=(SGLT2.pred-DPP4SU.pred)*100))$obs)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.02795  0.50804  0.89648  1.30174  1.61432 13.13283 

summary(cohort$qdhf_sglt2_benefit*100)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.02814  0.55854  1.00996  1.36487  1.79742 14.09242 




cal_plot <- ggplot(data=plot.data, aes(x=pred,y=obs.gp)) +
  geom_point(alpha=1) + theme_bw() +
  geom_errorbar(aes(ymin=l_iqr, ymax=u_iqr), colour="black", width=.1) +
  ylab("Observed 5-year SGLT2i benefit (%)") + xlab("Predicted 5-year SGLT2i benefit (%)") +
  scale_x_continuous(limits=c(0,4.5),breaks=c(seq(0,4.5,by=1))) +
  scale_y_continuous(limits=c(0,5),breaks=c(seq(0,5,by=1))) +
  theme_base() + 
  theme_bw() +
  theme(text = element_text(size = 18),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        #plot.margin = unit(c(0.5, 0, 0, 0), "cm"),
        legend.position = "none") +
  geom_abline(intercept=0,slope=1, color="red", lwd=0.75) + ggtitle("") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey60") + geom_hline(yintercept=0, linetype="dashed", color = "grey60")


cal_plot <- cal_plot +
  annotate(geom="text", x=0.1, y=4.8, label="Overall median observed benefit (IQR): 0.90% (0.51-1.6%)\nOverall median predicted benefit: 1.0%", color="black", hjust=0, size = 5)


hist_plot <- ggplot(cohort, aes(x=qdhf_sglt2_benefit*100)) +
  geom_histogram(aes(y = after_stat(count / sum(count))*100), binwidth=0.1, alpha=0.5, position="identity") +
  guides(fill = FALSE) +
  theme(legend.title = element_blank(), panel.background = element_rect( fill = "white",color = "grey50")) +
  scale_x_continuous(limits=c(0,4.5),breaks=c(seq(0,4.5,by=1))) +
  xlab(expression(paste("Predicted 5-year SGLT2i benefit (%)"))) +
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




setwd("/slade/CPRD_data/Katie SGLT2/Plots/")
tiff("Figure_2.tiff", width=8, height=8, units = "in", res=600) 

plot_grid(cal_plot, hist_plot, ncol = 1,align = 'v',
          rel_heights = c(1,0.5), rel_widths = c(1,1))

dev.off()  








# ### with studydrug*qdhf interaction - not using
# 
# 
# ## Observed - adjusted and weighted with overlap weights, with studydrug*decile interaction
# 
# ps.formula <- formula(paste0("studydrug ~ ", return_cov_set("weight")))
# 
# overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort), weight="overlap")
# cohort$overlap_weights <- overlap$ps.weights$overlap
# 
# # Remove constant variables otherwise datadist has issues
# cohort <- cohort %>% select(patid, studydrug, qdiabeteshf_5yr_score, dstartdate_age, ethnicity_qrisk2_decoded, malesex, dstartdate_dm_dur_all, imd2015_10, drugline_all, ncurrtx_cat, INS, initiation_year, prebmi, prehba1c2yrs, presbp, qrisk2_smoking_cat, hypertension, hosp_admission_prev_year_count, hf_censtime_yrs, hf_censvar, overlap_weights, predicted_benefit_decile, predrug_af)
# 
# ddist <- datadist(cohort)
# options(datadist='ddist')
# 
# f_adjusted <- as.formula(paste0("Surv(hf_censtime_yrs, hf_censvar) ~ studydrug*rcs(qdiabeteshf_5yr_score,5)  +  ", return_cov_set("adjust")))
# 
# model <- cph(f_adjusted, data=cohort, weights=overlap_weights, x=T, y=T, surv=TRUE)
# 
# cohort.t <- data.frame(cohort) %>% mutate(studydrug="DPP4SU")
# DPP4SU.pred = survest(model,newdata=cohort.t,times=5,se.fit=F)
# cohort.t <- data.frame(cohort) %>% mutate(studydrug="SGLT2")
# SGLT2.pred = survest(model,newdata=cohort.t,times=5,se.fit=F)
# 
# cohort$DPP4SU.pred <- DPP4SU.pred$surv
# cohort$SGLT2.pred <- SGLT2.pred$surv
# 
# obs <- cohort %>%
#   mutate(obs=(SGLT2.pred-DPP4SU.pred)*100) %>%
#   group_by(predicted_benefit_decile) %>%
#   dplyr::summarise(
#     obs.gp = median(obs,na.rm=T),
#     l_iqr=quantile(obs,0.25,na.rm=T),
#     u_iqr=quantile(obs,0.75,na.rm=T)
#   )
# 
# 
# predicted.p <- predicted %>% 
#   mutate(pred=pred.median*100) %>%
#   select(predicted_benefit_decile,pred)
# plot.data <- merge(obs,predicted.p,by="predicted_benefit_decile") 
# 
# 
# setwd("/slade/CPRD_data/Katie SGLT2/Plots/")
# tiff("calibration_benefits_qdhf_interaction.tiff", width=7, height=6, units = "in", res=800) 
# 
# ggplot(data=plot.data,aes(x=pred,y=obs.gp)) +
#   geom_point(alpha=1) + theme_bw() +
#   geom_errorbar(aes(ymin=l_iqr, ymax=u_iqr), colour="black", width=.1) +
#   ylab("Observed 5-year SGLT2i benefit (%)") + xlab("Predicted 5-year SGLT2i benefit (%)") +
#   scale_x_continuous(limits=c(0,4.5),breaks=c(seq(0,4.5,by=1))) +
#   scale_y_continuous(limits=c(0,5.5),breaks=c(seq(0,5.5,by=1))) +
#   theme_base() + 
#   theme(plot.background = element_blank()) +
#   geom_abline(intercept=0,slope=1, color="red", lwd=0.75) + ggtitle("") +
#   geom_vline(xintercept=0, linetype="dashed", color = "grey60") + geom_hline(yintercept=0, linetype="dashed", color = "grey60")
# 
# dev.off()  




