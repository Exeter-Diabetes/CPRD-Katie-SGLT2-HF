
# Baseline characteristics: unweighted and weighted


############################################################################################

# Setup
library(tidyverse)
library(table1)
library(PSweight)
library(cobalt)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# 1 Treatment outcome cohort

setwd("/slade/CPRD_data/Katie SGLT2/Processed data/")
load("treatment_outcome_cohort_jun24.rda")
#169,041


## Overall median age
summary(cohort$dstartdate_age)
#58.0 (50.4-66.1)

## Overall sex split
prop.table(table(cohort$malesex))
#57.7% male


## Overall ethnicity split
prop.table(table(cohort$ethnicity_qrisk2_decoded))
#74.6% white



## DPP4SU split
table(cohort$drugclass)
# DPP4 SGLT2    SU 
# 68708 57368 42965 


cat <- function(x, ...) {
  c("", sapply(stats.apply.rounding(stats.default(x)), function(y) with(y, sprintf("%s (%s%%)", prettyNum(FREQ, big.mark=","), PCT))))
}

cont <- function(x) {
  with(stats.apply.rounding(stats.default(x)), c("Median (IQR)"=sprintf("%s (%s-%s)", round_pad(as.numeric(MEDIAN),1), round_pad(as.numeric(Q1),1), round_pad(as.numeric(Q3),1))))
}

missing <- function(x, ...) {
  with(stats.apply.rounding(stats.default(x)), c("Missing"=sprintf("%s", prettyNum(NMISS, big.mark=","))))
}

strat <- function (label, n, ...) {
  sprintf("<span class='stratlabel'>%s<br><span class='stratn'>(N=%s)</span></span>", 
          label, prettyNum(n, big.mark=","))
}

rndr <- function(x, name, ...) {
  y <- render.default(x, name, ...)
  if (is.logical(x)) {
    y[2]
  } else {
    y
  }
}

# Reduce ethnicity to 5 category, smoking to 3 category, and IMD to quintiles for this

cohort <- cohort %>%
  mutate(ethnicity_decoded=factor(case_when(ethnicity_qrisk2_decoded=="missing" ~"Missing",
                                     ethnicity_qrisk2_decoded=="White" ~"White",
                                     ethnicity_qrisk2_decoded=="Indian" | ethnicity_qrisk2_decoded=="Pakistani" |  ethnicity_qrisk2_decoded=="Bangladeshi" | ethnicity_qrisk2_decoded=="Other Asian" ~"Asian",
                                     ethnicity_qrisk2_decoded=="Black Caribbean" | ethnicity_qrisk2_decoded=="Black African" ~"Black",
                                     ethnicity_qrisk2_decoded=="Chinese" | ethnicity_qrisk2_decoded=="Other" ~"Other"), levels=c("White", "Asian", "Black", "Other", "Missing")),
         smoker_decoded=factor(case_when(is.na(qrisk2_smoking_cat) ~as.character(NA),
                                  qrisk2_smoking_cat==0 ~"Non-smoker",
                                  qrisk2_smoking_cat==1 ~"Ex-smoker",
                                  qrisk2_smoking_cat==2 | qrisk2_smoking_cat==3 | qrisk2_smoking_cat==4 ~"Active smoker"), levels=c("Non-smoker", "Active smoker", "Ex-smoker")),
         imd_quintiles=as.factor(case_when(is.na(imd2015_10) ~as.character(NA),
                                 imd2015_10==1 | imd2015_10==2 ~"1 (least deprived)",
                                 imd2015_10==3 | imd2015_10==4 ~"2",
                                 imd2015_10==5 | imd2015_10==6 ~"3",
                                 imd2015_10==7 | imd2015_10==8 ~"4",
                                 imd2015_10==9 | imd2015_10==10 ~"5 (most deprived)")))
  

# Add labels
label(cohort$malesex)               <- "Sex (% male)"
label(cohort$dstartdate_age)        <- "Age at drug initiation (years)"
label(cohort$dstartdate_dm_dur_all) <- "Diabetes duration (years)"
label(cohort$ethnicity_decoded)     <- "Ethnicity"
label(cohort$imd_quintiles)         <- "Index of Multiple Deprivation quintile"
label(cohort$smoker_decoded)        <- "Smoking status"
label(cohort$hypertension)          <- "Hypertension"
label(cohort$predrug_af)            <- "Atrial fibrillation"
label(cohort$hosp_admission_prev_year_count) <- "Number of hospital admissions in previous year"
label(cohort$prebmi)                <- "BMI (kg/m2)"
label(cohort$prehba1c2yrs)          <- "HbA1c (mmol/mol)"
label(cohort$presbp)                <- "SBP (mmHg)"
label(cohort$precholhdl)            <- "Cholesterol:HDL"
label(cohort$drugline_all)          <- "Drug line"
label(cohort$ncurrtx_cat)           <- "Number of other current non-insulin glucose-lowering medications"
label(cohort$INS)                   <- "Current insulin use"
label(cohort$initiation_year)       <- "Year of drug initiation"
label(cohort$qdiabeteshf_5yr_score) <- "QDiabetes-Heart Failure 5-year score (%)"


setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions/")
source("full_covariate_set.R")

return_cov_set("weight")
"malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_qrisk2_decoded + imd2015_10 + qrisk2_smoking_cat + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score"

# Include all variables used for weighting/adjustment, plus chol:HDL (not used as has missingness)
# Use simplified versions of ethnicity, IMD and smoking

table1(~ malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_decoded + imd_quintiles + smoker_decoded + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + precholhdl + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score | studydrug, data=cohort, overall=F, render=rndr, render.categorical=cat, render.continuous=cont, render.strat=strat)
#superscript kg/m2
#add comma to missing cholesterol:HDL count
#Add DPP4 vs SU count


# Table with overlap weighting
# Did initially manage to code up for proportions, but easier to do this way for medians

ps.formula <- formula(paste0("studydrug ~ ", return_cov_set("weight")))
                      
overlap <- SumStat(ps.formula=ps.formula, data=as.data.frame(cohort), weight="overlap")
cohort$overlap_weight <- overlap$ps.weights$overlap
rm(overlap)
gc()

weighted_cohort <- cohort %>%
  select(studydrug, malesex, dstartdate_age, dstartdate_dm_dur_all, ethnicity_decoded, imd_quintiles, smoker_decoded, hypertension, predrug_af, hosp_admission_prev_year_count, prebmi, prehba1c2yrs, presbp, precholhdl, drugline_all, ncurrtx_cat, INS, initiation_year, qdiabeteshf_5yr_score, overlap_weight) %>%
  mutate(count=round(overlap_weight*10000000),0) %>% 
  uncount(count)
  

strat <- function (label, n, ...) {
  sprintf("<span class='stratlabel'>%s</span>", 
          label, prettyNum(n, big.mark=","))
}

cat <- function(x, ...) {
  c("", sapply(stats.apply.rounding(stats.default(x)), function(y) with(y, sprintf("%s%%", PCT))))
}


# Same variables as above
table1(~ malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_decoded + imd_quintiles + smoker_decoded + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + precholhdl + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score | studydrug, data=weighted_cohort, overall=F, render=rndr, render.categorical=cat, render.continuous=cont, render.strat=strat)
#need to remove missing counts from Chol:HDL afterwards
#superscript kg/m2


# Love plot

cohort <- cohort %>%
  select(studydrug, malesex, dstartdate_age, dstartdate_dm_dur_all, ethnicity_decoded, imd_quintiles, smoker_decoded, hypertension, predrug_af, hosp_admission_prev_year_count, prebmi, prehba1c2yrs, presbp, precholhdl, drugline_all, ncurrtx_cat, INS, initiation_year, qdiabeteshf_5yr_score, overlap_weight)

label(cohort$studydrug) <- "studydrug"
label(cohort$overlap_weight) <- "overlap_weight"

new_names <- data.frame(names=names(cohort), labels=(cohort %>% map_chr(attr_getter("label"))))
    

tiff("/slade/CPRD_data/Katie SGLT2/Plots/love_plot.tiff", width=8, height=10, units = "in", res=600)          
       
love.plot(studydrug ~ malesex + dstartdate_age + dstartdate_dm_dur_all + ethnicity_decoded + imd_quintiles + smoker_decoded + hypertension + predrug_af + hosp_admission_prev_year_count + prebmi + prehba1c2yrs + presbp + precholhdl + drugline_all + ncurrtx_cat + INS + initiation_year + qdiabeteshf_5yr_score, data = cohort, weights = cohort$overlap_weight,
          method = "weighting", estimand = "ATE",  var.names = new_names, sample.names = c("Unweighted", "Overlap weighted"), binary="std") +
  theme(legend.position="bottom",
        axis.text=element_text(size=11),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        plot.title=element_text(size=16))

dev.off()


############################################################################################

# Look at how close biomarkers are to drug start date on average

data <- cohort %>% 
  mutate(preacrdrugdiff_new=ifelse(!is.na(preacr), preacrdrugdiff, preacr_from_separatedrugdiff)) %>%
  select(patid, dstartdate, studydrug, prehba1cdrugdiff, prebmidrugdiff, pretotalcholesteroldrugdiff, prehdldrugdiff, presbpdrugdiff, predbpdrugdiff, preacrdrugdiff_new)
  
summary(data$prehba1cdrugdiff)
#median=15 days

summary(data$prebmidrugdiff)
#median=27 days

summary(data$pretotalcholesteroldrugdiff)
#median=32 days

summary(data$prehdldrugdiff)
#median=35 days

summary(data$presbpdrugdiff)
#median=13 days

summary(data$predbpdrugdiff)
#median=13 days

summary(data$preacrdrugdiff_new)
#median=123 days

data <- data %>%
  pivot_longer(cols=c(prehba1cdrugdiff, prebmidrugdiff, pretotalcholesteroldrugdiff, prehdldrugdiff, presbpdrugdiff, predbpdrugdiff, preacrdrugdiff_new))

summary(data$value)
#median=24 days

data %>% filter(!is.na(value)) %>% count()
#1129305

data %>% filter(!is.na(value) & value>=-183) %>% count()
#944122

944122/1129305
#83.6%



