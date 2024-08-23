
# Plot KM curves for control vs SGLT2

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

cohort <- add_surv_vars(cohort, main_only=TRUE)
         

## C Just keep variables of interest

cohort <- cohort %>%
  
  select(patid, malesex, ethnicity_qrisk2_decoded, imd2015_10, regstartdate, gp_record_end, death_date, drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx_cat, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, prehba1c2yrs, prebmi, prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, preegfr, preckdstage, contains("cens"), qdiabeteshf_5yr_score, qdiabeteshf_lin_predictor, initiation_year, hosp_admission_prev_year_count, hypertension, qrisk2_smoking_cat, drugorder, qrisk2_5yr_score, qrisk2_10yr_score, qrisk2_lin_predictor, predrug_af, tds_2011, uacr)

rm(list=setdiff(ls(), "cohort"))


# now using functions to define covariates adjusted and weighted for
setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions")
source("full_covariate_set.R")


############################################################################################

# 2 Calculate SGLT2 benefits

cohort <- cohort %>%
  mutate(qdhf_survival=(100-qdiabeteshf_5yr_score)/100,
         qdhf_survival_sglt2=qdhf_survival^0.63,
         qdiabeteshf_5yr_score_sglt2=100-(qdhf_survival_sglt2*100),
         qdhf_sglt2_benefit=qdhf_survival_sglt2-qdhf_survival)


############################################################################################

# 3 Km plots by absolute HF risk

## Split into groups based on absolute risk

cohort <- cohort %>% mutate(predicted_benefit_strata=cut(qdhf_sglt2_benefit*100, breaks = c(-Inf, 0.5, 1, 2, 3, Inf), labels = c(1, 2, 3, 4, 5)))
  

### Patient count by strata
table(cohort$predicted_benefit_strata)
# 1     2     3     4     5 
# 35989 47741 50004 20584 14723 
prop.table(table(cohort$predicted_benefit_strata))
# 1          2          3          4          5 
# 0.21290101 0.28242261 0.29580989 0.12176927 0.08709721 


### Median predicted by strata
cohort %>% group_by(predicted_benefit_strata) %>% summarise(pred.median=median(qdhf_sglt2_benefit,na.rm=T)*100)
# predicted_benefit_strata pred.median
# 1 1                              0.333
# 2 2                              0.729
# 3 3                              1.39 
# 4 4                              2.39 
# 5 5                              3.80 


### Observed benefit by strata

### Previously did unweighted/unadjusted: just
#model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug*predicted_benefit_strata, data=cohort, x=TRUE, y=TRUE, surv=TRUE)

### Now weighting within strata as per plots themselves
ps.formula <- formula(paste0("studydrug ~ ", return_cov_set("weight")))

cohort_1 <- cohort %>% filter(predicted_benefit_strata==1)
overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort_1), weight="overlap")
cohort_1$overlap_weights <- overlap$ps.weights$overlap

cohort_2 <- cohort %>% filter(predicted_benefit_strata==2)
overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort_2), weight="overlap")
cohort_2$overlap_weights <- overlap$ps.weights$overlap

cohort_3 <- cohort %>% filter(predicted_benefit_strata==3)
overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort_3), weight="overlap")
cohort_3$overlap_weights <- overlap$ps.weights$overlap

cohort_4 <- cohort %>% filter(predicted_benefit_strata==4)
overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort_4), weight="overlap")
cohort_4$overlap_weights <- overlap$ps.weights$overlap

cohort_5 <- cohort %>% filter(predicted_benefit_strata==5)
overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort_5), weight="overlap")
cohort_5$overlap_weights <- overlap$ps.weights$overlap


ddist <- datadist(cohort_1)
options(datadist='ddist')
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, weights=overlap_weights, data=cohort_1, x=TRUE, y=TRUE, surv=TRUE)

obs_SGLT2 <- survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv
obs_DPP4 <- survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv           
(obs_SGLT2-obs_DPP4)*100 
#0.2975442%              
              
ddist <- datadist(cohort_2)
options(datadist='ddist')
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, weights=overlap_weights, data=cohort_2, x=TRUE, y=TRUE, surv=TRUE)

obs_SGLT2 <- survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv
obs_DPP4 <- survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv           
(obs_SGLT2-obs_DPP4)*100 
#0.9258248%  

ddist <- datadist(cohort_3)
options(datadist='ddist')
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, weights=overlap_weights, data=cohort_3, x=TRUE, y=TRUE, surv=TRUE)

obs_SGLT2 <- survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv
obs_DPP4 <- survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv           
(obs_SGLT2-obs_DPP4)*100 
#1.298041%  

ddist <- datadist(cohort_4)
options(datadist='ddist')
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, weights=overlap_weights, data=cohort_4, x=TRUE, y=TRUE, surv=TRUE)

obs_SGLT2 <- survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv
obs_DPP4 <- survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv           
(obs_SGLT2-obs_DPP4)*100 
#1.703508%  

ddist <- datadist(cohort_5)
options(datadist='ddist')
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, weights=overlap_weights, data=cohort_5, x=TRUE, y=TRUE, surv=TRUE)

obs_SGLT2 <- survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv
obs_DPP4 <- survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv           
(obs_SGLT2-obs_DPP4)*100 
#3.627542%


## Km plots

cohort <- rbind(cohort_1, cohort_2, cohort_3, cohort_4, cohort_5)

survfit_list <- lapply(split(cohort, f=cohort$predicted_benefit_strata), function(x) survfit(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, data=x, weights=x$overlap_weights))

names(survfit_list[[1]][["strata"]]) <- c("DPP4i/SU   ", "SGLT2i")


p1 <- ggsurvplot_list(survfit_list,
                      data=cohort,
                      title=c("Predicted benefit <0.5%", "Predicted benefit 0.5-1%", "Predicted benefit 1-2%", "Predicted benefit 2-3%", "Predicted benefit >3%"),
                      size= 1.5,
                      fun = function(x) {100 - x*100},
                      conf.int = T,
                      ylim = c(0, 16),
                      xlab = "Years",
                      ylab = "Heart failure incidence (%)",
                      break.time.by = 1,
                      ggtheme = theme_classic(),
                      font.x = c(18),font.y = c(18),font.tickslab = c(16),
                      axes.offset = FALSE, # start the axis at the origin
                      palette = c("#782e6a","#E69F00"))

p1[[1]][["plot"]][["labels"]][["subtitle"]] <- "n=35,989 (21.3%)\nMedian 5-year predicted benefit: 0.33%\nObserved 5-year benefit: 0.30%"
p1[[2]][["plot"]][["labels"]][["subtitle"]] <- "n=47,741 (28.2%)\nMedian 5-year predicted benefit: 0.73%\nObserved 5-year benefit: 0.93%"
p1[[3]][["plot"]][["labels"]][["subtitle"]] <- "n=50,004 (29.6%)\nMedian 5-year predicted benefit: 1.4%\nObserved 5-year benefit: 1.3%"
p1[[4]][["plot"]][["labels"]][["subtitle"]] <- "n=20,584 (12.2%)\nMedian 5-year predicted benefit: 2.4%\nObserved 5-year benefit: 1.7%"
p1[[5]][["plot"]][["labels"]][["subtitle"]] <- "n=14,723 (8.7%)\nMedian 5-year predicted benefit: 3.8%\nObserved 5-year benefit: 3.6%"


setwd("/slade/CPRD_data/Katie SGLT2/Processed data/Plots/")
tiff("km_plots.tiff", width=15, height=11, units = "in", res=800)

grid.arrange(arrangeGrob(p1[[1]][["plot"]] + theme(legend.position=c(0.5, 0.7), legend.title=element_blank(), legend.text=element_text(size=16, face="bold"), plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)) + guides(colour = guide_legend(nrow = 1)),
                         p1[[2]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
                         p1[[3]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
                         p1[[4]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
                         p1[[5]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)), ncol=3, nrow=2, widths=c(1,1,1)))

#grid.text("Decile of predicted SGLT2i benefit:", x = unit(0.5, "npc"), y = unit(.99, "npc"), just=c("centre", "top"), gp = gpar(fontsize=22, fontface="bold"))

dev.off()



