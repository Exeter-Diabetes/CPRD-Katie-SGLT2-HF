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
library(forestplot)
library(boot)

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

## Matched to ADA/EASD

table(cohort$ada_easd)
#84682 treated

table((cohort %>% mutate(sabre=ifelse((qdhf_sglt2_benefit*100)>1.00787, 1, 0)))$sabre)
#84682 treated

table((cohort %>% mutate(qrisk2=ifelse(qrisk2_10yr_score>19.0686, 1, 0)))$qrisk2)
#84682 treated

test <- cohort %>% group_by(ada_easd) %>% summarise(hfb.u.iqr=quantile(qdhf_sglt2_benefit*100,0.75),
                                                    cvdb.u.iqr=quantile(qrisk2_10yr_score,0.75))

nrow(cohort %>% filter(ada_easd==1))/4
#21170.5

table((cohort %>% mutate(sabre=ifelse((qdhf_sglt2_benefit*100)>2.582255, 1, 0)))$sabre)
#21171 treated

table((cohort %>% mutate(sabre=ifelse(qrisk2_10yr_score>36.262, 1, 0)))$sabre)
#22171 treated

cohort <- cohort %>% mutate(sabre_matched_adaeasd=ifelse((qdhf_sglt2_benefit*100)>1.00787, 1, 0),
                            qrisk2_matched_adaeasd=ifelse(qrisk2_10yr_score>19.0686, 1, 0),
                            sabre_restricted=ifelse((qdhf_sglt2_benefit*100)>2.582255, 1, 0),
                            qrisk2_restricted=ifelse(qrisk2_10yr_score>36.262, 1, 0))

## Matched to NICE

table(cohort$nice_qrisk2_10)
#134880 treated

table((cohort %>% mutate(sabre=ifelse((qdhf_sglt2_benefit*100)>0.483705, 1, 0)))$sabre)
#134880 treated

cohort <- cohort %>% mutate(sabre_matched_nice=ifelse((qdhf_sglt2_benefit*100)>0.483705, 1, 0))


############################################################################################

# 5 Figures for plot

## Absolute risk reduction (benefit) per 100 patient-years and NNT as per: https://dom-pubs.pericles-prod.literatumonline.com/doi/full/10.1111/dom.14893)

## Bootstrap for 95% CIs

cohort <- cohort %>% mutate(treat_all=1)

strategies <- c("treat_all", "nice_qrisk2_10", "ada_easd", "sabre_matched_nice", "sabre_matched_adaeasd", "sabre_restricted", "qrisk2_matched_adaeasd", "qrisk2_restricted")

table <- data.frame(strategy=character(), treat=character(), arr=numeric(), arr_lower_ci=numeric(), arr_upper_ci=numeric(), nnt=character())

boot_fn <- function(data, indices) {

  d <- data[indices, ]
  
  n <- nrow(d)
  
  # Compute events under untreated and treated scenarios
  events_if_untreated <- round(n * mean(d$qdiabeteshf_5yr_score / 100))
  events_if_treated <- round(n * mean(d$qdiabeteshf_5yr_score_sglt2 / 100))
  
  # Absolute risk reduction (ARR) per 100 person-years
  arr <- ((events_if_untreated - events_if_treated) / (5 * n)) * 100
  
  # NNT per year to prevent 1 event
  nnt <- (1/(arr))*100
  
  return(c(arr = arr, nnt = nnt))
}


for (i in strategies) {
  
  cohort_to_treat <- cohort %>% filter(!!(as.symbol(i))==1)
  
  # % treated
  treat_n <- unlist(cohort_to_treat %>% count())
  treat <- paste0(round_pad((treat_n/169041)*100,1), "%")
  
  # Bootstrap for CIs
  set.seed(123)  # For reproducibility
  boot_results <- boot(data = cohort_to_treat, statistic = boot_fn, R = 1000)
  
  arr <- boot.ci(boot_results, type = "perc", index = 1)$t0
  arr_lower_ci <-  boot.ci(boot_results, type = "perc", index = 1)$percent[4]
  arr_upper_ci <-  boot.ci(boot_results, type = "perc", index = 1)$percent[5]
  
  nnt <- boot.ci(boot_results, type = "perc", index = 2)$t0
  #nnt_lower_ci <-  boot.ci(boot_results, type = "perc", index = 2)$percent[4]
  #nnt_upper_ci <-  boot.ci(boot_results, type = "perc", index = 2)$percent[5]
  
  nnt <- round_pad(nnt, 0) #paste0(round_pad(nnt, 0)," (", round_pad(nnt_lower_ci, 0), "-", round_pad(nnt_upper_ci, 0), ")")
  
  data <- data.frame(strategy=i, treat, arr, arr_lower_ci, arr_upper_ci, nnt)
  
  table <- table %>% rbind(data)
  
}


# Plot

table <- table %>% mutate(arr_text=paste0(round_pad(arr, 2), " (", round_pad(arr_lower_ci, 2), "-", round_pad(arr_upper_ci, 2), ")")) %>% rename(mean=arr, lower=arr_lower_ci, upper=arr_upper_ci)

tiff("/slade/CPRD_data/Katie SGLT2/Plots/strategy_benefit.tiff", width=16, height=6, units = "in", res=400) 

table %>% forestplot(labeltext = list(analysis_text=list("Treat all", expression("UK NICE guidance "^1), expression("ADA/EASD guidance "^2), expression("SABRE model matched to NICE "^3), expression("SABRE model matched to ADA/EASD "^4), expression("SABRE model restricted strategy "^5), expression("QRISK2 matched to ADA/EASD "^6), expression("QRISK2 restricted strategy "^7)), treat_text=table$treat, arr_text=table$arr_text, nnt=table$nnt),
                     ci.vertices = TRUE,
                     #xlab = "Events prevented per 100 patient-years",
                     fn.ci_norm = fpDrawCircleCI,
                     xticks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     boxsize = 0.25,
                     graph.pos=3,
                     lwd.zero = 1.2,
                     lwd.ci=4,
                     ci.vertices.height = 0.2,
                     colgap = unit(5, "mm"),
                     mar = unit(c(15,2,2,2), "mm"),
                     txt_gp = fpTxtGp(xlab=gpar(cex=1.5, fontface="bold", lineheight=10), ticks=gpar(cex=1.4), label=gpar(cex=1.4))) %>%
  fp_add_header(analysis_text = "Strategy", treat_text = "Proportion of\npopulation treated", arr_text = "", nnt = "NNT")


grid.text("Absolute risk for HF in SGLT2i treated group", 
          x = 0.7, 
          y = unit(1, "npc") - unit(1.1, "lines"), 
          gp = gpar(fontsize = 19, fontface = "bold"))

grid.text("Events prevented per 100 patient-year (95% CI)", 
          x = 0.7, 
          y = 0.05, 
          gp = gpar(fontsize = 19, fontface = "bold"))

dev.off()




############################################################################################

# 6 Km plots - weight within subgroups

source("full_covariate_set.R")
ps.formula <- formula(paste0("studydrug ~ ", return_cov_set("weight")))


## Divide into subgroups and weight

strategies <- c("nice_qrisk2_10", "ada_easd", "sabre_matched_nice", "sabre_matched_adaeasd", "sabre_restricted", "qrisk2_matched_adaeasd", "qrisk2_restricted")

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
# 84359              84682             147870              21170 

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
                        paste0("SGLT2i not recommended\n(n=147,870 [", round(100*(nrow(sabre_restricted_N)/nrow(cohort)),1),"%])"),
                        paste0("SGLT2i recommended\n(n=21,170 [", round(100*(nrow(sabre_restricted_Y)/nrow(cohort)),1),"%])")),
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
p1[[3]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n0.9% (95% CI: 0.5 to 1.4%)"
p1[[4]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n3.4% (95% CI: 1.3 to 5.6%)"


plot <- grid.arrange(arrangeGrob(
  p1[[1]][["plot"]] + theme(legend.position=c(0.5, 0.5), legend.title=element_blank(), legend.text=element_text(size=16, face="bold"), plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5, margin=margin(b=5, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)) + guides(colour = guide_legend(nrow = 1)),
  p1[[2]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5, margin=margin(b=5, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
  p1[[3]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)),
  p1[[4]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)), ncol=2, nrow=2, widths=c(1,1)))


tiff("/slade/CPRD_data/Katie SGLT2/Plots/main_strategies.tiff", width=10, height=12, units = "in", res=800)

as_ggplot(plot) +    
  draw_plot_label(label = c("a) ADA/EASD guidance strategy", "b) SABRE model restricted strategy"), size = 20,  hjust=0, x = c(0, 0), y = c(1, 0.475))

dev.off()





## NICE, SABRE-matched NICE, SABRE matched ADAEASD, QRISK2 matched ADAEASD, QRISK2 restricted

other_strategies <- rbind(nice_qrisk2_10_N, nice_qrisk2_10_Y, sabre_matched_nice_N, sabre_matched_nice_Y, sabre_matched_adaeasd_N, sabre_matched_adaeasd_Y, qrisk2_matched_adaeasd_N, qrisk2_matched_adaeasd_Y, qrisk2_restricted_N, qrisk2_restricted_Y) %>%
  mutate(subgp=factor(subgp, levels=c("nice_qrisk2_10_N", "nice_qrisk2_10_Y", "sabre_matched_nice_N", "sabre_matched_nice_Y", "sabre_matched_adaeasd_N", "sabre_matched_adaeasd_Y", "qrisk2_matched_adaeasd_N", "qrisk2_matched_adaeasd_Y", "qrisk2_restricted_N", "qrisk2_restricted_Y")))


table(other_strategies$subgp)
# nice_qrisk2_10_N         nice_qrisk2_10_Y     sabre_matched_nice_N     sabre_matched_nice_Y  sabre_matched_adaeasd_N  sabre_matched_adaeasd_Y 
# 34161                   134880                    34161                   134880                    84359                    84682 
# qrisk2_matched_adaeasd_N qrisk2_matched_adaeasd_Y      qrisk2_restricted_N      qrisk2_restricted_Y 
# 84359                    84682                   147870                    21171 

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
names(survfit_list$sabre_matched_adaeasd_N$strata) <- c("DPP4/SU", "SGLT2")

p1 <- ggsurvplot_list(survfit_list,
                      data=main_strategies,
                      title=c(
                        paste0("SGLT2i not recommended\n(n=34,161 [", round_pad(100*(nrow(nice_qrisk2_10_N)/nrow(cohort)),1),"%])"), 
                        paste0("SGLT2i recommended\n(n=134,880 [", round_pad(100*(nrow(nice_qrisk2_10_Y)/nrow(cohort)),1),"%])"), 
                        paste0("SGLT2i not recommended\n(n=34,161 [", round_pad(100*(nrow(sabre_matched_nice_N)/nrow(cohort)),1),"%])"),
                        paste0("SGLT2i recommended\n(n=134,880 [", round_pad(100*(nrow(sabre_matched_nice_Y)/nrow(cohort)),1),"%])"),
                        paste0("SGLT2i not recommended\n(n=84,359 [", round_pad(100*(nrow(sabre_matched_adaeasd_N)/nrow(cohort)),1),"%])"),
                        paste0("SGLT2i recommended\n(n=84,682 [", round_pad(100*(nrow(sabre_matched_adaeasd_Y)/nrow(cohort)),1),"%])"),
                        paste0("SGLT2i not recommended\n(n=84,359 [", round_pad(100*(nrow(qrisk2_matched_adaeasd_N)/nrow(cohort)),1),"%])"),
                        paste0("SGLT2i recommended\n(n=84,682 [", round_pad(100*(nrow(qrisk2_matched_adaeasd_Y)/nrow(cohort)),1),"%])"),
                        paste0("SGLT2i not recommended\n(n=147,870 [", round_pad(100*(nrow(qrisk2_restricted_N)/nrow(cohort)),1),"%])"),
                        paste0("SGLT2i recommended\n(n=21,171 [", round_pad(100*(nrow(qrisk2_restricted_Y)/nrow(cohort)),1),"%])")),
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
p1[[3]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n-0.2% (95% CI: -0.8 to 0.4%)"
p1[[4]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n1.6% (95% CI: 1.0 to 2.2%)"
p1[[5]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n0.5% (95% CI: 0.0 to 1.0%)"
p1[[6]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n1.8% (95% CI: 1.0 to 2.6%)"
p1[[7]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n0.6% (95% CI: 0.1 to 1.0%)"
p1[[8]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n1.8% (95% CI: 1.0 to 2.6%)"
p1[[9]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n1.1% (95% CI: 0.6 to 1.5%)"
p1[[10]][["plot"]][["labels"]][["subtitle"]] <- "Observed 5-year benefit with SGLT2i:\n2.6% (95% CI: 0.5 to 4.7%)"

plot <- grid.arrange(arrangeGrob(
  p1[[1]][["plot"]] + theme(legend.position=c(0.5, 0.5), legend.title=element_blank(), legend.text=element_text(size=16, face="bold"), plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5, margin=margin(b=5, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)) + guides(colour = guide_legend(nrow = 1)),
  p1[[2]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5, margin=margin(b=5, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
  p1[[3]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)),
  p1[[4]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)),
  ncol=2, nrow=2, widths=c(1,1)))
 

tiff("/slade/CPRD_data/Katie SGLT2/Plots/other_strategies_A.tiff", width=10, height=10, units = "in", res=800)

as_ggplot(plot) +    
  draw_plot_label(label = c("a) UK NICE guidance strategy", "b) SABRE model matched to NICE"), size = 20,  hjust=0, x = c(0, 0), y = c(1, 0.47))

dev.off()



plot <- grid.arrange(arrangeGrob(
  p1[[5]][["plot"]] + theme(legend.position=c(0.5, 0.5), legend.title=element_blank(), legend.text=element_text(size=16, face="bold"), plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5, margin=margin(b=5, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)) + guides(colour = guide_legend(nrow = 1)),
  p1[[6]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)),
  p1[[7]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)),
  p1[[8]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)),
  p1[[9]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)),
  p1[[10]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, vjust=-0.5,  margin=margin(b=-6, t=40), size = 18), plot.subtitle=element_text(hjust=0, vjust=-10, margin=margin(b=-25, l=10), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5, t=24)),
  ncol=2, nrow=3, widths=c(1,1)))



tiff("/slade/CPRD_data/Katie SGLT2/Plots/other_strategies_B.tiff", width=10, height=15, units = "in", res=800)

as_ggplot(plot) +    
  draw_plot_label(label = c("c) SABRE model matched to ADA/EASD", "d) QRISK2 model matched to ADA/EASD", "e) QRISK2 restricted strategy"), size = 20,  hjust=0, x = c(0, 0, 0), y = c(1, 0.647, 0.315))

dev.off()
