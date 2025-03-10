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
library(table1)

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


# 4. Find SABRE and QRISK2 thresholds for matched and restricted strategies

cohort %>% group_by(ada_easd) %>% summarise(n=n(),
                                            pc=n()/nrow(.),
                                            hfb.mean=mean(qdhf_sglt2_benefit*100),
                                            hfb=median(qdhf_sglt2_benefit*100),
                                            hfb.l=min(qdhf_sglt2_benefit*100),
                                            hfb.u=max(qdhf_sglt2_benefit*100),
                                            hfb.l.iqr=quantile(qdhf_sglt2_benefit*100,0.25),
                                            hfb.u.iqr=quantile(qdhf_sglt2_benefit*100,0.75),
                                            cvdb.u.iqr=quantile(qrisk2_10yr_score,0.75))
                                                   
# 50.1% treated by ADA-EASD
# SABRE matched: 1.01%
# SABRE: restricted: use upper IQR: 2.44%
# QRISK2 restricted: use upper IQR: 34.2%

cohort <- cohort %>% mutate(sabre_matched=ifelse((qdhf_sglt2_benefit*100)>1.01, 1, 0),
                            sabre_restricted=ifelse((qdhf_sglt2_benefit*100)>2.44, 1, 0),
                            qrisk2_restricted=ifelse(qrisk2_10yr_score>34.2, 1, 0))


############################################################################################

# 5 Figures for table

# ## Total number of HF cases in whole study population over 5 years based on mean HF risk
# describe(cohort$qdiabeteshf_5yr_score) #Mean HR risk 3.78%
# hf.notx <- round(nrow(cohort)*mean(cohort$qdiabeteshf_5yr_score/100)) 
# hf.notx
# #6382
# 
# cohort <- cohort %>% mutate(treat_all=1)
# 
# strategies <- c("treat_all", "nice_qrisk2_10", "ada_easd", "sabre_matched", "sabre_restricted", "qrisk2_restricted")
# 
# table <- data.frame(strategy=character(), treat=character(), hf_cases=character(), hf_prevented=character(), nnt=numeric())
# 
# 
# for (i in strategies) {
#   
#   treat_n <- unlist(cohort %>% filter(!!(as.symbol(i))==1) %>% count())
#   treat_perc <- round_pad((treat_n/169041)*100,1)
#   treat <- paste0(treat_n, " (", treat_perc, "%)")
#   
#   cohort <- cohort %>% mutate(qdiabeteshf_5yr_score.applied=ifelse(!!(as.symbol(i))==0, qdiabeteshf_5yr_score, qdiabeteshf_5yr_score_sglt2))
#   hf_cases_n <- round(169041*mean(cohort$qdiabeteshf_5yr_score.applied/100))
#   hf_cases_perc <- round_pad((hf_cases_n/169041)*100,1)
#   hf_cases <- paste0(hf_cases_n, " (", hf_cases_perc, "%)")
#   
#   hf_prevented_n <- 6382-hf_cases_n
#   hf_prevented_perc <- round_pad((hf_prevented_n/2307)*100,1)
#   hf_prevented <- paste0(hf_prevented_n, " (", hf_prevented_perc, "%)")
#   
#   nnt <- round(treat_n/hf_prevented_n, 0)
#   
#   data <- cbind(strategy=i, treat, hf_cases, hf_prevented, nnt)
#   
#   table <- table %>% rbind(data)
#   
# }


############################################################################################

# 6 Km plots - weight within subgroups

source("full_covariate_set.R")
ps.formula <- formula(paste0("studydrug ~ ", return_cov_set("weight")))


## Divide into subgroups and weight

strategies <- c("nice_qrisk2_10", "ada_easd", "sabre_matched", "sabre_restricted", "qrisk2_restricted")

for (i in strategies) {
  
  not_treated_group_name <- paste0(i, "_N")
  treated_group_name <- paste0(i, "_Y")
  
  treated_group <- cohort %>% filter(!!(as.symbol(i))==1) %>% mutate(subgp=treated_group_name)
  overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(treated_group), weight="overlap")
  treated_group$overlap_weights <- overlap$ps.weights$overlap
  assign(treated_group_name, treated_group)
  
  not_treated_group <- cohort %>% filter(!!(as.symbol(i))==0) %>% mutate(subgp=not_treated_group_name)
  overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(not_treated_group), weight="overlap")
  not_treated_group$overlap_weights <- overlap$ps.weights$overlap
  assign(not_treated_group_name, not_treated_group)
}


## ADA-EASD vs SABRE-restricted

main_strategies <- rbind(ada_easd_N, ada_easd_Y, sabre_restricted_N, sabre_restricted_Y) %>%
  mutate(subgp=factor(subgp, levels=c("ada_easd_N", "ada_easd_Y", "sabre_restricted_N", "sabre_restricted_Y")))


table(main_strategies$subgp)
# ada_easd_N         ada_easd_Y sabre_restricted_N sabre_restricted_Y 
# 84359              84682             145087              23954 

survfit_list <- lapply(split(main_strategies, f=main_strategies$subgp), function(x) survfit(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, data=x, weights=x$overlap_weights))


### Difference in survival at 5 years with confidence intervals

table <- data.frame(group=character(), surv_diff=character())

for (i in unique(main_strategies$subgp)) {
  test <- data.frame(time=summary(survfit_list[[i]])$time, surv=summary(survfit_list[[i]])$surv, lower=summary(survfit_list[[i]])$lower, upper=summary(survfit_list[[i]])$upper, strata=summary(survfit_list[[i]])$strata)
  
  closest_index_dpp4su <- which.max((test %>% filter(strata=="studydrug=DPP4SU"))$time)
  surv_dpp4su <- (1-test$surv[closest_index_dpp4su])*100
  surv_upper_dpp4su <- (1-test$lower[closest_index_dpp4su])*100
  surv_lower_dpp4su <- (1-test$upper[closest_index_dpp4su])*100
  
  closest_index_sglt2 <- which.max((test %>% filter(strata=="studydrug=SGLT2"))$time) + closest_index_dpp4su
  surv_sglt2 <- (1-test$surv[closest_index_sglt2])*100
  surv_upper_sglt2 <- (1-test$lower[closest_index_sglt2])*100
  surv_lower_sglt2 <- (1-test$upper[closest_index_sglt2])*100

  data <- data.frame(group=i, surv_dpp4su, surv_lower_dpp4su, surv_upper_dpp4su, surv_sglt2, surv_lower_sglt2, surv_upper_sglt2) %>%
    mutate(surv_diff=surv_dpp4su-surv_sglt2,
           SE1=(surv_upper_dpp4su - surv_lower_dpp4su) / (2 * 1.96),
           SE2=(surv_upper_sglt2 - surv_lower_sglt2) / (2 * 1.96),
           SE_diff=sqrt(SE1^2 + SE2^2),
           CI_lower=surv_diff - 1.96 * SE_diff,
           CI_upper=surv_diff + 1.96 * SE_diff,
           surv_diff=paste0(round_pad(surv_diff, 1), "% (95% CI: ", round_pad(CI_lower, 1), " to ",  round_pad(CI_upper, 1), "%)")) %>%
    select(group, surv_diff)
  
  table <- table %>% rbind(data)
  
}

### Km plot

names(survfit_list$ada_easd_N$strata) <- c("DPP4/SU", "SGLT2")

p1 <- ggsurvplot_list(survfit_list,
                      data=main_strategies,
                      title=c(
                        paste0("SGLT2i not recommended\n(n=84,359 [", round(100*(nrow(ada_easd_N)/nrow(cohort)),1),"%])"), 
                        paste0("SGLT2i recommended\n(n=84,682 [", round(100*(nrow(ada_easd_Y)/nrow(cohort)),1),"%])"), 
                        paste0("SGLT2i not recommended\n(n=145,087 [", round(100*(nrow(sabre_restricted_N)/nrow(cohort)),1),"%])"),
                        paste0("SGLT2i recommended\n(n=23,954 [", round(100*(nrow(sabre_restricted_Y)/nrow(cohort)),1),"%])")),
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


p1[[1]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n0.7% (95% CI: 0.1 to 1.3%)"
p1[[2]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n1.7% (95% CI: 0.9 to 2.4%)"
p1[[3]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n0.9% (95% CI: 0.5 to 1.3%)"
p1[[4]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n3.0% (95% CI: 1.0 to 5.0%)"


plot <- grid.arrange(arrangeGrob(
  p1[[1]][["plot"]] + theme(legend.position=c(0.5, 0.5), legend.title=element_blank(), legend.text=element_text(size=16, face="bold"), plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5, margin=margin(b=5, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)) + guides(colour = guide_legend(nrow = 1)),
  p1[[2]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5, margin=margin(b=5, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
  p1[[3]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)),
  p1[[4]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)), ncol=2, nrow=2, widths=c(1,1)))


tiff("/slade/CPRD_data/Katie SGLT2/Plots/main_strategies.tiff", width=10, height=12, units = "in", res=800)

as_ggplot(plot) +    
  draw_plot_label(label = c("a) International guidance strategy", "b) SABRE model restricted strategy"), size = 20,  hjust=0, x = c(0, 0), y = c(1, 0.475))

dev.off()





## NICE, SABRE-matched, QRISK2 restricted

other_strategies <- rbind(nice_qrisk2_10_N, nice_qrisk2_10_Y, sabre_matched_N, sabre_matched_Y, qrisk2_restricted_N, qrisk2_restricted_Y) %>%
  mutate(subgp=factor(subgp, levels=c("nice_qrisk2_10_N", "nice_qrisk2_10_Y", "sabre_matched_N", "sabre_matched_Y", "qrisk2_restricted_N", "qrisk2_restricted_Y")))


table(other_strategies$subgp)
# nice_qrisk2_10_N    nice_qrisk2_10_Y     sabre_matched_N     sabre_matched_Y qrisk2_restricted_N qrisk2_restricted_Y 
# 34161              134880               84524               84517              143184               25857 

survfit_list <- lapply(split(other_strategies, f=other_strategies$subgp), function(x) survfit(Surv(hf_censtime_yrs, hf_censvar) ~ studydrug, data=x, weights=x$overlap_weights))


### Difference in survival at 5 years with confidence intervals

table <- data.frame(group=character(), surv_diff=character())

for (i in unique(other_strategies$subgp)) {
  test <- data.frame(time=summary(survfit_list[[i]])$time, surv=summary(survfit_list[[i]])$surv, lower=summary(survfit_list[[i]])$lower, upper=summary(survfit_list[[i]])$upper, strata=summary(survfit_list[[i]])$strata)
  
  closest_index_dpp4su <- which.max((test %>% filter(strata=="studydrug=DPP4SU"))$time)
  surv_dpp4su <- (1-test$surv[closest_index_dpp4su])*100
  surv_upper_dpp4su <- (1-test$lower[closest_index_dpp4su])*100
  surv_lower_dpp4su <- (1-test$upper[closest_index_dpp4su])*100
  
  closest_index_sglt2 <- which.max((test %>% filter(strata=="studydrug=SGLT2"))$time) + closest_index_dpp4su
  surv_sglt2 <- (1-test$surv[closest_index_sglt2])*100
  surv_upper_sglt2 <- (1-test$lower[closest_index_sglt2])*100
  surv_lower_sglt2 <- (1-test$upper[closest_index_sglt2])*100
  
  data <- data.frame(group=i, surv_dpp4su, surv_lower_dpp4su, surv_upper_dpp4su, surv_sglt2, surv_lower_sglt2, surv_upper_sglt2) %>%
    mutate(surv_diff=surv_dpp4su-surv_sglt2,
           SE1=(surv_upper_dpp4su - surv_lower_dpp4su) / (2 * 1.96),
           SE2=(surv_upper_sglt2 - surv_lower_sglt2) / (2 * 1.96),
           SE_diff=sqrt(SE1^2 + SE2^2),
           CI_lower=surv_diff - 1.96 * SE_diff,
           CI_upper=surv_diff + 1.96 * SE_diff,
           surv_diff=paste0(round_pad(surv_diff, 1), "% (95% CI: ", round_pad(CI_lower, 1), " to ",  round_pad(CI_upper, 1), "%)")) %>%
    select(group, surv_diff)
  
  table <- table %>% rbind(data)
  
}

### Km plot

names(survfit_list$nice_qrisk2_10_N$strata) <- c("DPP4/SU", "SGLT2")

p1 <- ggsurvplot_list(survfit_list,
                      data=main_strategies,
                      title=c(
                        paste0("SGLT2i not recommended\n(n=34,161 [", round_pad(100*(nrow(nice_qrisk2_10_N)/nrow(cohort)),1),"%])"), 
                        paste0("SGLT2i recommended\n(n=134,880 [", round_pad(100*(nrow(nice_qrisk2_10_Y)/nrow(cohort)),1),"%])"), 
                        paste0("SGLT2i not recommended\n(n=84,524 [", round_pad(100*(nrow(sabre_matched_N)/nrow(cohort)),1),"%])"),
                        paste0("SGLT2i recommended\n(n=84,517 [", round_pad(100*(nrow(sabre_matched_Y)/nrow(cohort)),1),"%])"),
                        paste0("SGLT2i not recommended\n(n=143,184 [", round_pad(100*(nrow(qrisk2_restricted_N)/nrow(cohort)),1),"%])"),
                        paste0("SGLT2i recommended\n(n=25,857 [", round_pad(100*(nrow(qrisk2_restricted_Y)/nrow(cohort)),1),"%])")),
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


p1[[1]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n0.2% (95% CI: -0.5 to 0.9%)"
p1[[2]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n1.5% (95% CI: 1.0 to 2.1%)"
p1[[3]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n0.5% (95% CI: 0.0 to 1.0%)"
p1[[4]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n1.8% (95% CI: 1.0 to 2.6%)"
p1[[5]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n1.0% (95% CI: 0.5 to 1.5%)"
p1[[6]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n2.6% (95% CI: 0.7 to 4.4%)"

plot <- grid.arrange(arrangeGrob(
  p1[[1]][["plot"]] + theme(legend.position=c(0.5, 0.5), legend.title=element_blank(), legend.text=element_text(size=16, face="bold"), plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5, margin=margin(b=5, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)) + guides(colour = guide_legend(nrow = 1)),
  p1[[2]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5, margin=margin(b=5, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
  p1[[3]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)),
  p1[[4]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)),
  p1[[5]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)),
  p1[[6]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)), ncol=2, nrow=3, widths=c(1,1)))


tiff("/slade/CPRD_data/Katie SGLT2/Plots/other_strategies.tiff", width=10, height=14, units = "in", res=800)

as_ggplot(plot) +    
  draw_plot_label(label = c("a) UK NICE guidance strategy", "b) SABRE model matched strategy", "c) QRISK2 restricted strategy"), size = 20,  hjust=0, x = c(0, 0, 0), y = c(1, 0.647, 0.315))

dev.off()

