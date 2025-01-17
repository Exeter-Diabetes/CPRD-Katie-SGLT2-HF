############################################################################################

# Setup
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(gridExtra)
library(cowplot)
library(grid)
library(PSweight)
library(rms)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# 1 Cohort selection and variable setup

setwd("/slade/CPRD_data/Katie SGLT2/Processed data")
load("treatment_outcome_cohort_jun24.rda")
#169,041

table(cohort$studydrug)
# DPP4SU 111673   
# SGLT2 57368  


## B Make variables for survival analysis of all endpoints (see survival_variables function for details)

setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions")
source("survival_variables.R")

cohort <- add_surv_vars(cohort, main_only=TRUE)

## C Just keep variables of interest

cohort <- cohort %>%
  
  select(patid, malesex, ethnicity_qrisk2_decoded, imd2015_10, regstartdate, gp_record_end, death_date, drugclass, studydrug, dstartdate, dstopdate, drugline_all, 
         drugsubstances, ncurrtx_cat, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, prehba1c2yrs, prebmi, prehdl, preldl, 
         pretriglyceride, pretotalcholesterol, prealt, presbp, preegfr, preckdstage, contains("cens"), qdiabeteshf_5yr_score, qdiabeteshf_lin_predictor, initiation_year, 
         hosp_admission_prev_year_count, hypertension, qrisk2_smoking_cat, drugorder, qrisk2_5yr_score, qrisk2_10yr_score, qrisk2_lin_predictor, predrug_af, tds_2011, uacr,
         predrug_hypertension, predbp, statins, preacr, preacr_from_separate)

rm(list=setdiff(ls(), "cohort"))

############################################################################################

# 2 Add SGLT2 benefits

cohort <- cohort %>%
  mutate(qdhf_survival=(100-qdiabeteshf_5yr_score)/100,
         qdhf_survival_sglt2=qdhf_survival^0.63,
         qdiabeteshf_5yr_score_sglt2=100-(qdhf_survival_sglt2*100),
         qdhf_sglt2_benefit=qdhf_survival_sglt2-qdhf_survival)


############################################################################################

# 3 Add in ADA/EASD and NICE categories

# ADA/EASD

## Obesity and smoking straightforward

## Hypertension
# 4 have missing DBP (0.002%) - definitely don't exclude as SBP is main measure
# 13785 have DBP>90, 38025 have SBP>140, 140808 have hypertension (from GP records - all time - use? Will no longer have high BP if treated, but are using ever)
# End up with 64% of cohort with hypertension, 69% through GP diagnosis code alone

## Dyslipidemia
# Definitions all use total cholesterol, many also use LDL/non-HDL cholesterol (we don't have latter and former is very missing - 37%), some also use HDL (minimum level)
# Have used total cholesterol and statin usage here: 56617 have chol >5, 163454 on statins, 20619 both
# Not sure if statin usage is valid as given to everyone..., but then most people won't have high cholesterol after starting statin
# 78% with dyslipidemia, 72% of these from statin usage alone
# 1488 have missing total cholesterol (<1%)

## Albuminuria
# Assume not present if UACR measurement missing


cohort <- cohort %>%
  mutate(obesity=ifelse(prebmi>=30, 1, 0),
         hypertension=ifelse(presbp>=140 | predrug_hypertension==1 | (!is.na(predbp) & predbp>90), 1, 0),
         smoking=ifelse(as.numeric(qrisk2_smoking_cat)>2, 1, 0),
         dyslipidemia=ifelse((!is.na(pretotalcholesterol) & pretotalcholesterol>5) | statins==1, 1, 0),
         
         uacr=ifelse(!is.na(preacr), preacr, preacr_from_separate),
         albuminuria=ifelse(!is.na(uacr) & uacr>3, 1, 0),
         
         ada_easd=ifelse(dstartdate_age>=55 & (obesity+hypertension+smoking+dyslipidemia+albuminuria)>=2, 1, 0))


# NICE
cohort <- cohort %>%
  mutate(nice_qrisk2_10=ifelse(qrisk2_10yr_score>10, 1, 0))


############################################################################################

# 4 Define subgroups and weight within subgroup

## A ADA/EASD vs QDHF

cohort %>% group_by(ada_easd) %>% dplyr::summarise(n=n(),
                                                   pc=n()/nrow(.),
                                                   hfb.mean=mean(qdhf_sglt2_benefit*100),
                                                   hfb=median(qdhf_sglt2_benefit*100),
                                                   hfb.l=min(qdhf_sglt2_benefit*100),
                                                   hfb.u=max(qdhf_sglt2_benefit*100),
                                                   hfb.l.iqr=quantile(qdhf_sglt2_benefit*100,0.25),
                                                   hfb.u.iqr=quantile(qdhf_sglt2_benefit*100,0.75))


source("full_covariate_set.R")
ps.formula <- formula(paste0("studydrug ~ ", return_cov_set("weight")))

#ADA
cohort_adaN <- cohort %>% filter(ada_easd==0) %>% mutate(subgp="cohort_adaN")
overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort_adaN), weight="overlap")
cohort_adaN$overlap_weights <- overlap$ps.weights$overlap

cohort_adaY <- cohort %>% filter(ada_easd==1) %>% mutate(subgp="cohort_adaY")
overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort_adaY), weight="overlap")
cohort_adaY$overlap_weights <- overlap$ps.weights$overlap


#ADA QDHF 75th percentile
ada.threshold <- 2.44 #hfb.u.iqr in ada/easd treat group

cohort_ada_qdhfN <- cohort %>% filter((100*qdhf_sglt2_benefit)<ada.threshold) %>% mutate(subgp="cohort_ada_qdhfN")
overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort_ada_qdhfN), weight="overlap")
cohort_ada_qdhfN$overlap_weights <- overlap$ps.weights$overlap

cohort_ada_qdhfY <- cohort %>% filter((100*qdhf_sglt2_benefit)>=ada.threshold) %>% mutate(subgp="cohort_ada_qdhfY")
overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort_ada_qdhfY), weight="overlap")
cohort_ada_qdhfY$overlap_weights <- overlap$ps.weights$overlap

describe(cohort_ada_qdhfN$qdhf_sglt2_benefit*100)
describe(cohort_ada_qdhfY$qdhf_sglt2_benefit*100)



#Observed benefits

##cohort_adaN
ddist <- datadist(cohort_adaN);options(datadist='ddist')
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, weights=overlap_weights, data=cohort_adaN, x=TRUE, y=TRUE, surv=TRUE)
obs_SGLT2 <- survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv
obs_DPP4 <- survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv           
(obs_SGLT2-obs_DPP4)*100 #0.79%
(1-survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv) * 100 #1.84%
(1-survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv) * 100 #2.63%   

##cohort_adaY
ddist <- datadist(cohort_adaY);options(datadist='ddist')
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, weights=overlap_weights, data=cohort_adaY, x=TRUE, y=TRUE, surv=TRUE)
obs_SGLT2 <- survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv
obs_DPP4 <- survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv           
(obs_SGLT2-obs_DPP4)*100 #1.81%
(1-survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv) * 100 #4.32%
(1-survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv) * 100 #6.12%   

##cohort_ada_qdhfN
ddist <- datadist(cohort_ada_qdhfN);options(datadist='ddist')
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, weights=overlap_weights, data=cohort_ada_qdhfN, x=TRUE, y=TRUE, surv=TRUE)
obs_SGLT2 <- survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv
obs_DPP4 <- survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv           
(obs_SGLT2-obs_DPP4)*100 #1.01%
(1-survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv) * 100 #2.25%
(1-survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv) * 100 #3.26%   

##cohort_ada_qdhfY
ddist <- datadist(cohort_ada_qdhfY);options(datadist='ddist')
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, weights=overlap_weights, data=cohort_ada_qdhfY, x=TRUE, y=TRUE, surv=TRUE)
obs_SGLT2 <- survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv
obs_DPP4 <- survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv           
(obs_SGLT2-obs_DPP4)*100  #3.11%
(1-survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv) * 100 #8.55%
(1-survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv) * 100 #11/67%


#Km plots

new_cohort <- rbind(cohort_adaN, cohort_adaY, cohort_ada_qdhfN, cohort_ada_qdhfY)

table(new_cohort$subgp)
# cohort_adaN      cohort_adaY cohort_ada_qdhfN cohort_ada_qdhfY 
# 84359            84682           145087            23954 

new_cohort <- new_cohort %>% mutate(subgp=factor(subgp, levels=c("cohort_adaN", "cohort_adaY", "cohort_ada_qdhfN", "cohort_ada_qdhfY")))

survfit_list <- lapply(split(new_cohort, f=new_cohort$subgp), function(x) survfit(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, data=x, weights=x$overlap_weights))

names(survfit_list[[1]][["strata"]]) <- c("DPP4i/SU   ", "SGLT2i")


p1 <- ggsurvplot_list(survfit_list,
                      data=new_cohort,
                      title=c(
                        paste0("ADA/EASD: SGLT2i not\nrecommended (n=84,359 [", round(100*(nrow(cohort_adaN)/nrow(cohort)),1),"%])"), 
                        paste0("ADA/EASD: SGLT2i\nrecommended (n=84,682 [", round(100*(nrow(cohort_adaY)/nrow(cohort)),1),"%])"), 
                        paste0("SABRE model: SGLT2i not\nrecommended (n=145,087 [", round(100*(nrow(cohort_ada_qdhfN)/nrow(cohort)),1),"%])"),
                        paste0("SABRE model: SGLT2i\nrecommended (n=23,954 [", round(100*(nrow(cohort_ada_qdhfY)/nrow(cohort)),1),"%])")),
                      size= 1.5,
                      fun = function(x) {100 - x*100},
                      conf.int = T,
                      risk.table=TRUE,
                      ylim = c(0, 15),
                      xlab = "Years",
                      ylab = "Heart failure incidence (%)",
                      break.time.by = 1,
                      ggtheme = theme_classic(),
                      font.x = c(18),font.y = c(18),font.tickslab = c(16),
                      axes.offset = FALSE, # start the axis at the origin
                      palette = c("#782e6a","#E69F00"))


p1[[1]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit: 0.8%"
p1[[2]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit: 1.8%"
p1[[3]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit: 1.0%"
p1[[4]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit: 3.1%"

plot <- grid.arrange(arrangeGrob(
  p1[[1]][["plot"]] + theme(legend.position=c(0.5, 0.85), legend.title=element_blank(), legend.text=element_text(size=16, face="bold"), plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 18), plot.subtitle=element_text(hjust=0.2, vjust=-20, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=5)) + guides(colour = guide_legend(nrow = 1)),
  p1[[2]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 18), plot.subtitle=element_text(hjust=0.2, vjust=-20, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=5)),
  p1[[3]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 18), plot.subtitle=element_text(hjust=0.2, vjust=-20, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=5)),
  p1[[4]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 18), plot.subtitle=element_text(hjust=0.2, vjust=-20, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=5)), ncol=2, nrow=2, widths=c(1,1)))


tiff("/slade/CPRD_data/Katie SGLT2/Plots/adaeasd_vs_qdhf_km_plots.tiff", width=10, height=10, units = "in", res=800)

as_ggplot(plot) +    
  draw_plot_label(label = c("a)", "b)"), size = 16,  x = c(0, 0), y = c(1, 0.5)) #+    
  #draw_plot_label(label = c("(Predicted absolute benefit <2.4%)", "(Predicted absolute benefit \u22652.4%)"), size = 16,  x = c(-0.083, 0.417), y = c(0.45, 0.45))

dev.off()




## A ADA/EASD vs QRISK2

cohort %>% group_by(ada_easd) %>% dplyr::summarise(n=n(),
                                                   pc=n()/nrow(.),
                                                   hfb.mean=mean(qrisk2_10yr_score),
                                                   hfb=median(qrisk2_10yr_score),
                                                   hfb.l=min(qrisk2_10yr_score),
                                                   hfb.u=max(qrisk2_10yr_score),
                                                   hfb.l.iqr=quantile(qrisk2_10yr_score,0.25),
                                                   hfb.u.iqr=quantile(qrisk2_10yr_score,0.75))


source("full_covariate_set.R")
ps.formula <- formula(paste0("studydrug ~ ", return_cov_set("weight")))

#ADA
cohort_adaN <- cohort %>% filter(ada_easd==0) %>% mutate(subgp="cohort_adaN")
overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort_adaN), weight="overlap")
cohort_adaN$overlap_weights <- overlap$ps.weights$overlap

cohort_adaY <- cohort %>% filter(ada_easd==1) %>% mutate(subgp="cohort_adaY")
overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort_adaY), weight="overlap")
cohort_adaY$overlap_weights <- overlap$ps.weights$overlap


#ADA QDHF 75th percentile
ada.threshold <- 34.2 #hfb.u.iqr in ada/easd treat group

cohort_ada_qriskN <- cohort %>% filter(qrisk2_10yr_score<ada.threshold) %>% mutate(subgp="cohort_ada_qriskN")
overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort_ada_qriskN), weight="overlap")
cohort_ada_qriskN$overlap_weights <- overlap$ps.weights$overlap

cohort_ada_qriskY <- cohort %>% filter(qrisk2_10yr_score>=ada.threshold) %>% mutate(subgp="cohort_ada_qriskY")
overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort_ada_qriskY), weight="overlap")
cohort_ada_qriskY$overlap_weights <- overlap$ps.weights$overlap

describe(cohort_ada_qriskY$qdhf_sglt2_benefit*100)
describe(cohort_ada_qriskY$qdhf_sglt2_benefit*100)


#Observed benefits

##cohort_adaN
ddist <- datadist(cohort_adaN);options(datadist='ddist')
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, weights=overlap_weights, data=cohort_adaN, x=TRUE, y=TRUE, surv=TRUE)
obs_SGLT2 <- survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv
obs_DPP4 <- survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv           
(obs_SGLT2-obs_DPP4)*100 #0.79%
(1-survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv) * 100 #1.84%
(1-survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv) * 100 #2.63%   

##cohort_adaY
ddist <- datadist(cohort_adaY);options(datadist='ddist')
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, weights=overlap_weights, data=cohort_adaY, x=TRUE, y=TRUE, surv=TRUE)
obs_SGLT2 <- survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv
obs_DPP4 <- survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv           
(obs_SGLT2-obs_DPP4)*100 #1.81%
(1-survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv) * 100 #4.32%
(1-survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv) * 100 #6.12%   

##cohort_ada_qriskN
ddist <- datadist(cohort_ada_qriskN);options(datadist='ddist')
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, weights=overlap_weights, data=cohort_ada_qriskN, x=TRUE, y=TRUE, surv=TRUE)
obs_SGLT2 <- survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv
obs_DPP4 <- survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv           
(obs_SGLT2-obs_DPP4)*100 #1.15%
(1-survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv) * 100 #2.42%
(1-survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv) * 100 #3.57%   

##cohort_ada_qdhfY
ddist <- datadist(cohort_ada_qriskY);options(datadist='ddist')
model <- cph(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, weights=overlap_weights, data=cohort_ada_qriskY, x=TRUE, y=TRUE, surv=TRUE)
obs_SGLT2 <- survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv
obs_DPP4 <- survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv           
(obs_SGLT2-obs_DPP4)*100  #2.40%
(1-survest(model, newdata=expand.grid(studydrug="SGLT2"), times=5)$surv) * 100 #7.49%
(1-survest(model, newdata=expand.grid(studydrug="DPP4SU"), times=5)$surv) * 100 #9.89%


#Km plots

new_cohort <- rbind(cohort_adaN, cohort_adaY, cohort_ada_qriskN, cohort_ada_qriskY)

table(new_cohort$subgp)
# cohort_ada_qriskN cohort_ada_qriskY       cohort_adaN       cohort_adaY 
# 143184             25857             84359             84682 

new_cohort <- new_cohort %>% mutate(subgp=factor(subgp, levels=c("cohort_adaN", "cohort_adaY", "cohort_ada_qriskN", "cohort_ada_qriskY")))

survfit_list <- lapply(split(new_cohort, f=new_cohort$subgp), function(x) survfit(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, data=x, weights=x$overlap_weights))

names(survfit_list[[1]][["strata"]]) <- c("DPP4i/SU   ", "SGLT2i")


p1 <- ggsurvplot_list(survfit_list,
                      data=new_cohort,
                      title=c(
                        paste0("ADA/EASD: SGLT2i not\nrecommended (n=84,359 [",round(100*(nrow(cohort_adaN)/nrow(cohort)),1),"%])"), 
                        paste0("ADA/EASD: SGLT2i\nrecommended (n=84,682 [",round(100*(nrow(cohort_adaY)/nrow(cohort)),1),"%])"), 
                        paste0("QRISK2: SGLT2i not\nrecommended (n=143,184 [",round(100*(nrow(cohort_ada_qdhfN)/nrow(cohort)),1),"%])"),
                        paste0("QRISK2: SGLT2i\nrecommended (n=25,857 [",round(100*(nrow(cohort_ada_qdhfY)/nrow(cohort)),1),"%])")),
                      size= 1.5,
                      fun = function(x) {100 - x*100},
                      conf.int = T,
                      risk.table=TRUE,
                      ylim = c(0, 15),
                      xlab = "Years",
                      ylab = "Heart failure incidence (%)",
                      break.time.by = 1,
                      ggtheme = theme_classic(),
                      font.x = c(18),font.y = c(18),font.tickslab = c(16),
                      axes.offset = FALSE, # start the axis at the origin
                      palette = c("#782e6a","#E69F00"))


p1[[1]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit: 0.8%"
p1[[2]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit: 1.8%"
p1[[3]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit: 1.2%"
p1[[4]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit: 2.4%"

plot <- grid.arrange(arrangeGrob(
  p1[[1]][["plot"]] + theme(legend.position=c(0.5, 0.85), legend.title=element_blank(), legend.text=element_text(size=16, face="bold"), plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 18), plot.subtitle=element_text(hjust=0.2, vjust=-20, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=5)) + guides(colour = guide_legend(nrow = 1)),
  p1[[2]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 18), plot.subtitle=element_text(hjust=0.2, vjust=-20, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=5)),
  p1[[3]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 18), plot.subtitle=element_text(hjust=0.2, vjust=-20, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=5)),
  p1[[4]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 18), plot.subtitle=element_text(hjust=0.2, vjust=-20, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=5)), ncol=2, nrow=2, widths=c(1,1)))


tiff("/slade/CPRD_data/Katie SGLT2/Plots/adaeasd_vs_qrisk_km_plots.tiff", width=10, height=10, units = "in", res=800)

as_ggplot(plot) +    
  draw_plot_label(label = c("a)", "b)"), size = 16,  x = c(0, 0), y = c(1, 0.5)) #+    
#draw_plot_label(label = c("(Predicted absolute benefit \u22642.4%)", "(Predicted absolute benefit >2.4%)"), size = 16,  x = c(-0.083, 0.417), y = c(0.45, 0.45))

dev.off()



############################################################################################

# 5 Missed cases

## Total number of HF cases based on mean HF risk
describe(cohort$qdiabeteshf_5yr_score) #Mean HR risk 3.78%
hf.notx <- round(nrow(cohort)*mean(cohort$qdiabeteshf_5yr_score/100)) 
hf.notx
#6382

## If treat everyone with SGLT2i, new risk = qdiabeteshf_5yr_score_sglt2
describe(cohort$qdiabeteshf_5yr_score_sglt2) #Mean HR risk 2.41%  (same as 1-qdhf_survival_sglt2)
hf.tx <- round(nrow(cohort)*mean(cohort$qdiabeteshf_5yr_score_sglt2/100)) 
hf.tx
#4075

### So no. of HF events prevented
hf.notx - hf.tx #2307

##treated=169041
##cases prevented=2307
169041/2307 #73



#Treat NICE
describe(cohort$nice_qrisk2_10)

## New HF risk: lowered in those treated
cohort <- cohort %>% mutate(qdiabeteshf_5yr_score.applied = 
                              ifelse(nice_qrisk2_10==0, qdiabeteshf_5yr_score, qdiabeteshf_5yr_score_sglt2))
hf.tx <- round(nrow(cohort)*mean(cohort$qdiabeteshf_5yr_score.applied/100)) 
hf.tx #4209 cases
hf.notx - hf.tx   #2173 cases prevented

##treated=134894
##cases=2173
134894/2173 #62


#Treat ADA
describe(cohort$ada_easd)
cohort <- cohort %>% mutate(qdiabeteshf_5yr_score.applied = 
                              ifelse(ada_easd==0, qdiabeteshf_5yr_score, qdiabeteshf_5yr_score_sglt2))
hf.tx <- round(nrow(cohort)*mean(cohort$qdiabeteshf_5yr_score.applied/100)) 
hf.tx #4712
hf.notx - hf.tx #1670

##treated=84359
##cases=1670
84359/1670 #51


#Treat model 2.44
ada.threshold <- 2.44
cohort_modelADA.Y <- cohort %>% filter((100*qdhf_sglt2_benefit)>=ada.threshold) %>% mutate(subgp="cohort_modelADA.Y")
nrow(cohort_modelADA.Y)/nrow(cohort)

cohort <- cohort %>% mutate(qdiabeteshf_5yr_score.applied = 
                              ifelse(qdhf_sglt2_benefit<(2.44/100),qdiabeteshf_5yr_score,qdiabeteshf_5yr_score_sglt2))
mean(cohort$qdiabeteshf_5yr_score.applied)
hf.tx <- round(nrow(cohort)*mean(cohort$qdiabeteshf_5yr_score.applied/100)) 
hf.tx
hf.notx - hf.tx

#NNT
23954/877




#Treat qrisk2 34.2
ada.threshold <- 34.2
cohort_modelADA.Y <- cohort %>% filter(qrisk2_10yr_score>=ada.threshold) %>% mutate(subgp="cohort_modelADA.Y")
nrow(cohort_modelADA.Y)/nrow(cohort)
#15.3%
cohort_modelADA.Y %>% count()
#25,857

cohort <- cohort %>% mutate(qdiabeteshf_5yr_score.applied = 
                              ifelse(qrisk2_10yr_score<34.2, qdiabeteshf_5yr_score, qdiabeteshf_5yr_score_sglt2))
mean(cohort$qdiabeteshf_5yr_score.applied)
hf.tx <- round(nrow(cohort)*mean(cohort$qdiabeteshf_5yr_score.applied/100)) 
hf.tx
hf.notx - hf.tx

#NNT
25857/805




