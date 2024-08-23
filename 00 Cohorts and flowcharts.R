
# Flowchart for study inclusion for both treatment response and contemporary cohort


############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(DiagrammeR)
library(DiagrammeRsvg)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# Treatment outcome cohort

# Have to do initial steps on MySQL

cprd = CPRDData$new(cprdEnv = "diabetes-2020",cprdConf = "~/.aurum.yaml")

## T2Ds (with known diabetes duration)
## With HES linkage
## 1st instance
## Aged 18+
## SGLT2/DPP4/SU
## Initiated between 01/01/2013 and end of data (31/10/2020)

analysis = cprd$analysis("all")
diabetes_cohort <- diabetes_cohort %>% analysis$cached("diabetes_cohort")

analysis = cprd$analysis("mm")
drug_start_stop <- drug_start_stop %>% analysis$cached("drug_start_stop")
combo_start_stop <- combo_start_stop %>% analysis$cached("combo_start_stop")
baseline_biomarkers <- baseline_biomarkers %>% analysis$cached("baseline_biomarkers")
comorbidities <- comorbidities %>% analysis$cached("comorbidities")
ckd_stages <- ckd_stages %>% analysis$cached("ckd_stages")

original_cohort <- diabetes_cohort %>%
  select(patid, diabetes_type, dm_diag_date_all, with_hes, regstartdate, dob, imd2015_10) %>%
  inner_join((drug_start_stop %>% select(patid, dstartdate, drugclass, druginstance, drugline_all)), by="patid") %>%
  inner_join((combo_start_stop %>% select(patid, dcstartdate, INS, GLP1, SGLT2, TZD)), by=c("patid", "dstartdate"="dcstartdate")) %>%
  inner_join((baseline_biomarkers %>% select(patid, dstartdate, drugclass, preacr, preacr_from_separate)), by=c("patid", "dstartdate", "drugclass")) %>%
  inner_join((comorbidities %>% select(patid, dstartdate, drugclass, predrug_angina, predrug_ihd, predrug_myocardialinfarction, predrug_pad, predrug_revasc, predrug_stroke, predrug_tia, predrug_heartfailure)), by=c("patid", "dstartdate", "drugclass")) %>%
  inner_join((ckd_stages %>% select(patid, dstartdate, drugclass, preckdstage)), by=c("patid", "dstartdate", "drugclass")) %>%
  mutate(dstartdate_age=datediff(dstartdate, dob)/365.25) %>%
  filter(!is.na(dm_diag_date_all) &
           diabetes_type=="type 2" &
           with_hes==1 &
           druginstance==1 &
           dstartdate_age>=25 & dstartdate_age<=84 &
           (drugclass=="SGLT2" | drugclass=="DPP4" | drugclass=="SU") &
           dstartdate>=as.Date("2013-01-01")) %>%
  collect()
# 355,667 rows

original_cohort %>% distinct(patid) %>% count()
# 256,264 individuals

table(original_cohort$drugclass)
# DPP4  SGLT2     SU 
# 150207  97446 108014



## No CVD (I'm using a broad definition: angina, IHD, MI, PAD, revasc, stroke, TIA [as per NICE but with TIA])
## No HF before index date
## No CKD (stage 3a-5 based on eGFR; or ACR>=3) before index date (assume CKD stage coded on all, so if missing assume negative)

original_cohort %>% filter(predrug_angina==1 | predrug_ihd==1 | predrug_myocardialinfarction==1 | predrug_pad==1 | predrug_revasc==1 | predrug_stroke==1 | predrug_tia==1) %>% count()
#94258

original_cohort %>% filter(predrug_heartfailure==1) %>% count()
#24978

original_cohort %>% mutate(uacr=ifelse(!is.na(preacr), preacr, preacr_from_separate)) %>% filter((!is.na(preckdstage) &  (preckdstage=="stage_3a" | preckdstage=="stage_3b" | preckdstage=="stage_4" | preckdstage=="stage_5")) | (!is.na(uacr) & uacr>=30)) %>% count()
#50571

original_cohort %>% mutate(uacr=ifelse(!is.na(preacr), preacr, preacr_from_separate)) %>% filter(predrug_angina==1 | predrug_ihd==1 | predrug_myocardialinfarction==1 | predrug_pad==1 | predrug_revasc==1 | predrug_stroke==1 | predrug_tia==1 | predrug_heartfailure==1 | (!is.na(preckdstage) &  (preckdstage=="stage_3a" | preckdstage=="stage_3b" | preckdstage=="stage_4" | preckdstage=="stage_5")) | (!is.na(uacr) & uacr>=30)) %>% count()
#124865


cohort <- original_cohort %>%
  mutate(uacr=ifelse(!is.na(preacr), preacr, preacr_from_separate)) %>%
  filter(predrug_angina==0 & predrug_ihd==0 & predrug_myocardialinfarction==0 & predrug_pad==0 & predrug_revasc==0 & predrug_stroke==0 & predrug_tia==0 & predrug_heartfailure==0 & (is.na(preckdstage) | (preckdstage!="stage_3a" & preckdstage!="stage_3b" & preckdstage!="stage_4" & preckdstage!="stage_5")) & (is.na(uacr) | uacr<30))
#230802

cohort %>% distinct(patid) %>% count()
#163,476



## Exclude if start drug within 90 days of registration
cohort <- cohort %>% filter(difftime(dstartdate, regstartdate, units="days")>90)
#206,897
230802-206897
#23,905

cohort %>% distinct(patid) %>% count()
#150,954



## Exclude if first line
cohort <- cohort %>%
  filter(drugline_all!=1)
#195809
206897-195809
#11088

cohort %>% distinct(patid) %>% count()
#143365



# Remove if on GLP1/SGLT2 (except SGLT2 arm)/TZD (i above)
cohort <- cohort %>%
  filter(GLP1==0 & TZD==0 & (drugclass=="SGLT2" | SGLT2==0))
#181876
195809-181876
#13933

cohort %>% distinct(patid) %>% count()
#136775



## Use R data from here
setwd("/slade/CPRD_data/mastermind_2022/")
load("20240219_t2d_1stinstance.Rda")


# Keep those aged >=18 and within study period and second line or later (e-h above)
cohort <- t2d_1stinstance %>%
  mutate(uacr=ifelse(!is.na(preacr), preacr, preacr_from_separate)) %>%
  filter(dstartdate_age>=25 & dstartdate_age<=84 &
           (drugclass=="SGLT2" | drugclass=="DPP4" | drugclass=="SU") &
           dstartdate>=as.Date("2013-01-01") &
           drugline_all!=1 &
           GLP1==0 & TZD==0 & (drugclass=="SGLT2" | SGLT2==0) &
           predrug_angina==0 & predrug_ihd==0 & predrug_myocardialinfarction==0 & predrug_pad==0 & predrug_revasc==0 & predrug_stroke==0 & predrug_tia==0 & predrug_heartfailure==0 & (is.na(preckdstage) | (preckdstage!="stage_3a" & preckdstage!="stage_3b" & preckdstage!="stage_4" & preckdstage!="stage_5")) & (is.na(uacr) | uacr<30) & 
           !is.na(dm_diag_date_all))
#181876 as above



# Remove if don't have QDHF variables (also removes those without QRISK2)
## Need smoking and HbA1c
## And not to have HbA1c / age / BMI / SBP / cholHDL outside of normal range / range for model

# or variables for adjustment: 
## imd2015_10 (in QDHF but can be missing)
## prebmi (in QDHF but imputed if missing)
## prehba1c2yrs (in QDHF)
## presbp (in QDHF but imputed if missing)
## qrisk2_smoking_cat (in QDHF)


cohort %>% filter(is.na(prebmi)) %>% count()
#8876
cohort %>% filter(is.na(qrisk2_smoking_cat)) %>% count()
#1576
cohort %>% filter(is.na(prehba1c2yrs)) %>% count()
#1020
cohort %>% filter(is.na(presbp)) %>% count()
#1011
cohort %>% filter(is.na(imd2015_10)) %>% count()
#89

#missing anything
cohort %>% filter(is.na(prebmi) | is.na(qrisk2_smoking_cat) | is.na(prehba1c2yrs) | is.na(presbp) | is.na(imd2015_10)) %>% count()
#10,726

#anything out of range
cohort %>% mutate(precholhdl=pretotalcholesterol/prehdl) %>% filter((!is.na(presbp) & (presbp<70 | presbp>210)) | (!is.na(prehba1c2yrs) & (prehba1c2yrs<40 | prehba1c2yrs>150)) | (!is.na(precholhdl) & (precholhdl<1 | precholhdl>11.00000001)) | (!is.na(prebmi) & prebmi<20)) %>% count()
#2,262

#adds up to more than 12835 as some people have 1 thing missing and another thing out of range


cohort %>% filter(is.na(qdiabeteshf_5yr_score)) %>% count()
#4664
#missing QDHF

cohort %>% filter(!is.na(qdiabeteshf_5yr_score) & (is.na(imd2015_10) | is.na(prebmi) | is.na(presbp))) %>% count()
#8171
#not missing QDHF but missing BMI or SBP or deprivation

4664+8171
#12835


cohort %>% mutate(precholhdl=pretotalcholesterol/prehdl) %>% filter(is.na(qrisk2_smoking_cat) | is.na(prehba1c2yrs) | dstartdate_age<25 | dstartdate_age>84 | prebmi<20 | presbp<70 | presbp>210 | prehba1c2yrs<40 | prehba1c2yrs>150 | precholhdl<1 | precholhdl>11.00001 | is.na(imd2015_10) | is.na(prebmi) | is.na(presbp)) %>% count()
#12835 as above



cohort <- cohort %>% filter(!is.na(qdiabeteshf_5yr_score) & !is.na(imd2015_10) & !is.na(prebmi) & !is.na(presbp))
#169041
181876-12835 #169041

cohort %>% distinct(patid) %>% count()
#128270

table(cohort$drugclass)
# DPP4 SGLT2    SU 
# 68708 57368 42965 

cohort <- cohort %>% mutate(studydrug=relevel(factor(ifelse(drugclass=="SGLT2", "SGLT2", "DPP4SU")), ref="DPP4SU"))

table(cohort$studydrug)
# DPP4SU 111673      
# SGLT2 57368 


# Check matches cohort definition function
# Also formats variables correctly
setwd("/slade/CPRD_data/mastermind_2022/")
load("20240219_t2d_all_drug_periods.Rda")

setwd("/slade/CPRD_data/Katie SGLT2/Scripts/Functions")
source("cohort_definition.R")

cohort <- define_cohort(t2d_1stinstance, t2d_all_drug_periods)

cohort <- cohort %>% mutate(studydrug=relevel(factor(ifelse(drugclass=="SGLT2", "SGLT2", "DPP4SU")), ref="DPP4SU"))

table(cohort$studydrug)
# DPP4SU 111673   
# SGLT2 57368


# If patient has DPP4i and/or SU period, need to add dstartdate to all drug rows so that can censor on this

cohort <- cohort %>%
  left_join(cohort %>% filter(drugclass=="DPP4") %>% select(patid, dpp4_start_later = dstartdate), by="patid") %>%
  left_join(cohort %>% filter(drugclass=="SU") %>% select(patid, su_start_later = dstartdate), by="patid") %>%
  mutate(dpp4_start_later=if_else(dpp4_start_later<dstartdate | drugclass=="DPP4", as.Date(NA), dpp4_start_later),
         su_start_later=if_else(su_start_later<dstartdate | drugclass=="SU", as.Date(NA), su_start_later))


# Add in current statin treatment and hypertension definition

cohort <- cohort %>%
  mutate(statins=ifelse((!is.na(predrug_latest_statins) & predrug_latest_statins==dstartdate) | (!is.na(predrug_latest_statins) & !is.na(postdrug_first_statins) & as.numeric(difftime(dstartdate, predrug_latest_statins, units="days"))<=183 & as.numeric(difftime(postdrug_first_statins, predrug_latest_statins, units="days"))<=183), TRUE, FALSE),
         hypertension=ifelse(predrug_hypertension==1 | presbp>190 | (!is.na(predbp) & predbp>90), TRUE, FALSE))


setwd("/slade/CPRD_data/Katie SGLT2/Processed data/")
save(cohort, file="treatment_outcome_cohort_jun24.rda")



# Flowchart

plot <- "
digraph cohort_flow_chart
{
node [fontname = Helvetica, fontsize = 12, shape = box, width = 4]
a[label = 'T2D aged >=25 and <=84 years initiating drug\nSGLT2i, DPP4i or SU from 01/01/2013-31/10/2020\n(n=355,667 drug periods for 256,264 individuals;\nSGLT2i: 97,446; DPP4i: 150,207; SU:108,014)',  fontsize=14]
d[label = '\u2008\u2008\u2008\u2008\u2008\u2008CVD, HF and/or CKD at drug initiation:\u2008\u2008\u2008\u2008\u2008\u2008\n\u2022 CVD (n=94,258 drug periods)\n\u2022 Heart failure (n=24,978 drug periods)\n\u2022 CKD (n=50,571 drug periods)\n\u2022 Any (n=124,865 drug periods)',  fontsize=14]
e[label = 'No CVD/heart failure/CKD at drug initiation\n(n=230,802 drug periods for 163,476 individuals)',  fontsize=14]
f[label = 'Drug initiation date unreliable as within 90 days of registration\n(n=23,905 drug periods)',  fontsize=14]
g[label = 'Reliable date of drug initiation date\n(n=206,897 drug periods for 150,954 individuals)',  fontsize=14]
h[label = '\u2008\u2008\u2008Drug of interest taken as first-line diabetes therapy\u2008\u2008\u2008\n(n=11,088 drug periods)',  fontsize=14]
i[label = 'Drug of interest (SGLT2i/DPP4i/SU)\nnot taken as first-line diabetes therapy\n(n=195,809 drug periods for 143,365 individuals)',  fontsize=14]
j[label = 'Patient on TZD/GLP1/SGLT2i (excluding SGLT2i arm) at drug\ninitiation\n(n=13,933 drug periods)',  fontsize=14]
k[label = 'Patient not on another diabetes drug\nwith cardiovascular effect at drug initiation\n(TZD/GLP1/SGLT2i except if in SGLT2i arm;\nn=181,876 drug periods for 136,775 individuals)',  fontsize=14]
l[label = 'Baseline variables missing or out of range for cardiovascular\nrisk scores: \n\u2022 Missing BMI (n=8,876 drug periods)\n\u2022 Missing HbA1c (n=1,020 drug periods)\n\u2022 Missing SBP (n=1,011 drug periods)\n\u2022 Missing smoking status (n=1,576 drug periods)\n\u2022 Missing IMD (n=89 drug periods)\n\u2022 Any missing (n=10,726 drug periods)\n\u2022 Variable out of range (n=2,262 drug periods)\n\u2022 Any missing or out of range (n=12,835 drug periods)',  fontsize=14]
m[label = 'Final treatment outcome cohort\n(n=169,041 drug periods for 128,270 individuals;\nSGLT2i: 57,368, DPP4i: 68,708, SU: 42,965)',  fontsize=14]

{ rank = same; a, d}
{ rank = same; e, f}
{ rank = same; g, h}
{ rank = same; i, j}
{ rank = same; k, l}

a -> d [ minlen = 2 ];
a -> e;
e -> f [ minlen = 2 ];
e -> g;
g -> h [ minlen = 2 ];
g -> i;
i -> j [ minlen = 2 ];
i -> k;
k -> l [ minlen = 2 ];
k -> m;
}

[1]: excluded
"

setwd("/slade/CPRD_data/Katie SGLT2/Plots/")

DPI=800
grViz(plot) %>%
  export_svg() %>%
  charToRaw %>%
  rsvg::rsvg(width = 6*DPI, height = 5*DPI) %>% tiff::writeTIFF("flowdiagram.tiff", bits.per.sample = 8L)


############################################################################################

# Contemporary cohort

cprd = CPRDData$new(cprdEnv = "diabetes-2020",cprdConf = "~/.aurum.yaml")
analysis = cprd$analysis("prev")
prevalent_cohort <- prevalent_cohort %>% analysis$cached("final_merge")


prevalent_cohort <- prevalent_cohort %>%
  filter(!is.na(dm_diag_date_all) &
           diabetes_type=="type 2" &
           with_hes==1 &
           index_date_age>=25 & index_date_age<=84) %>%
  collect()
#503,913 
  

prevalent_cohort %>% filter(pre_index_date_angina==1 | pre_index_date_ihd==1 | pre_index_date_myocardialinfarction==1 | pre_index_date_pad==1 | pre_index_date_revasc==1 | pre_index_date_stroke==1 | pre_index_date_tia==1) %>% count()
#153882

prevalent_cohort %>% filter(pre_index_date_heartfailure==1) %>% count()
#45195

prevalent_cohort %>% mutate(uacr=ifelse(!is.na(preacr), preacr, preacr_from_separate)) %>% filter((!is.na(preckdstage) &  (preckdstage=="stage_3a" | preckdstage=="stage_3b" | preckdstage=="stage_4" | preckdstage=="stage_5")) | (!is.na(uacr) & uacr>=30)) %>% count()
#90849

prevalent_cohort %>% mutate(uacr=ifelse(!is.na(preacr), preacr, preacr_from_separate)) %>% filter(pre_index_date_angina==1 | pre_index_date_ihd==1 | pre_index_date_myocardialinfarction==1 | pre_index_date_pad==1 | pre_index_date_revasc==1 | pre_index_date_stroke==1 | pre_index_date_tia==1 | pre_index_date_heartfailure==1 | (!is.na(preckdstage) &  (preckdstage=="stage_3a" | preckdstage=="stage_3b" | preckdstage=="stage_4" | preckdstage=="stage_5")) | (!is.na(uacr) & uacr>=30)) %>% count()
#206151


cohort <- prevalent_cohort %>%
  mutate(uacr=ifelse(!is.na(preacr), preacr, preacr_from_separate)) %>%
  filter(pre_index_date_angina==0 & pre_index_date_ihd==0 & pre_index_date_myocardialinfarction==0 & pre_index_date_pad==0 & pre_index_date_revasc==0 & pre_index_date_stroke==0 & pre_index_date_tia==0 & pre_index_date_heartfailure==0 & (is.na(preckdstage) | (preckdstage!="stage_3a" & preckdstage!="stage_3b" & preckdstage!="stage_4" & preckdstage!="stage_5")) & (is.na(uacr) | uacr<30))
#297762


cohort %>% filter(is.na(prebmi)) %>% count()
#23408
cohort %>% filter(is.na(qrisk2_smoking_cat)) %>% count()
#3838
cohort %>% filter(is.na(prehba1c2yrs)) %>% count()
#8599
cohort %>% filter(is.na(presbp)) %>% count()
#7357
cohort %>% filter(is.na(imd2015_10)) %>% count()
#147

#missing anything
cohort %>% filter(is.na(prebmi) | is.na(qrisk2_smoking_cat) | is.na(prehba1c2yrs) | is.na(presbp) | is.na(imd2015_10)) %>% count()
#27924

#anything out of range
cohort %>% mutate(precholhdl=pretotalcholesterol/prehdl) %>% filter((!is.na(presbp) & (presbp<70 | presbp>210)) | (!is.na(prehba1c2yrs) & (prehba1c2yrs<40 | prehba1c2yrs>150)) | (!is.na(precholhdl) & (precholhdl<1 | precholhdl>11.00000001)) | (!is.na(prebmi) & prebmi<20)) %>% count()
#16830

#adds up to more than 12835 as some people have 1 thing missing and another thing out of range


cohort %>% filter(is.na(qdiabeteshf_5yr_score)) %>% count()
#27907
#missing QDHF

cohort %>% filter(!is.na(qdiabeteshf_5yr_score) & (is.na(imd2015_10) | is.na(prebmi) | is.na(presbp))) %>% count()
#15356
#not missing QDHF but missing BMI or SBP or deprivation

27907+15356
#43263



cohort %>% mutate(precholhdl=pretotalcholesterol/prehdl) %>% filter(is.na(qrisk2_smoking_cat) | is.na(prehba1c2yrs) | index_date_age<25 | index_date_age>84 | prebmi<20 | presbp<70 | presbp>210 | prehba1c2yrs<40 | prehba1c2yrs>150 | precholhdl<1 | precholhdl>11.00001 | is.na(imd2015_10) | is.na(prebmi) | is.na(presbp)) %>% count()
#43263 as above





cohort <- cohort %>% filter(!is.na(qdiabeteshf_5yr_score) & !is.na(imd2015_10) & !is.na(prebmi) & !is.na(presbp))
#254499
297762-43263 #254499


# Add in variables for current diabetes treatment - easier to do this way than by using latest pre_index_date and earliest post_index_date values (have checked using same method as for statins below, and gives same results)
setwd("/slade/CPRD_data/mastermind_2022/")
load("20240219_t2d_all_drug_periods.Rda")


treatment <- cohort %>%
  select(patid) %>%
  mutate(patid=as.character(patid)) %>%
  inner_join(t2d_all_drug_periods, by="patid") %>%
  filter(dstartdate<=as.Date("2020-02-01") & dstopdate>=as.Date("2020-02-01")) %>%
  select(patid, drugclass) %>%
  pivot_wider(patid, names_from=drugclass, values_from=drugclass, values_fn=~1, values_fill=0) %>%
  select(-c(Acarbose, Glinide)) 

cohort <- cohort %>%
  mutate(patid=as.character(patid)) %>%
  left_join(treatment, by="patid") %>%
  replace_na(list(DPP4=0, GLP1=0, MFN=0, SU=0, SGLT2=0, TZD=0, INS=0)) %>%
  mutate(ncurrtx=DPP4+GLP1+MFN+SU+SGLT2+TZD,
         ncurrtx_cat=as.factor(ifelse(ncurrtx<=1, ncurrtx, ifelse(ncurrtx>=2, "2+", NA))),
         DPP4=as.logical(DPP4),
         GLP1=as.logical(GLP1),
         MFN=as.logical(MFN),
         SU=as.logical(SU),
         SGLT2=as.logical(SGLT2),
         TZD=as.logical(TZD),
         INS=as.logical(INS))


# Add in current statin treatment and hypertension

cohort <- cohort %>%
  mutate(statins=ifelse((!is.na(pre_index_date_latest_statins) & pre_index_date_latest_statins==as.Date("2020-02-01")) | (!is.na(pre_index_date_latest_statins) & !is.na(post_index_date_first_statins) & as.numeric(difftime(as.Date("2020-02-01"), pre_index_date_latest_statins, units="days"))<=183 & as.numeric(difftime(post_index_date_first_statins, pre_index_date_latest_statins, units="days"))<=183), TRUE, FALSE),
         hypertension=ifelse(pre_index_date_hypertension==1 | presbp>190 | (!is.na(predbp) & predbp>90), TRUE, FALSE))


setwd("/slade/CPRD_data/Katie SGLT2/Processed data/")
save(cohort, file="contemporary_cohort_jun24.rda")



# Flowchart

plot <- "
digraph cohort_flow_chart
{
node [fontname = Helvetica, fontsize = 12, shape = box, width = 4]
a[label = 'T2D aged >=25 and <=84 years registered in\nprimary care on 1st February 2020\n(n=503,913)',  fontsize=14]
d[label = '\u2008\u2008\u2008\u2008\u2008CVD, HF and/or CKD at 1st February 2020:\u2008\u2008\u2008\u2008\u2008\n\u2022 CVD (n=153,882\n\u2022 Heart failure (n=45,195)\n\u2022 CKD (n=90,849)\n\u2022 Any (n=206,151)',  fontsize=14]
e[label = 'No CVD/heart failure/CKD at 1st February 2020\n(n=297,762)',  fontsize=14]
l[label = 'Baseline variables missing or out of range for cardiovascular\nrisk scores: \n\u2022 Missing BMI (n=23,408)\n\u2022 Missing HbA1c (n=8,599)\n\u2022 Missing SBP (n=7,357)\n\u2022 Missing smoking status (n=3,838)\n\u2022 Missing IMD (n=147)\n\u2022 Any missing (n=27,924)\n\u2022 Variable out of range (n=16,830)\n\u2022 Any missing or out of range (n=43,263)',  fontsize=14]
m[label = 'Final contemporary cohort\n(n=254,499)',  fontsize=14]

{ rank = same; a, d}
{ rank = same; e, l}

a -> d [ minlen = 2 ];
a -> e;
e -> l [ minlen = 2 ];
e -> m;
}

[1]: excluded
"

setwd("/slade/CPRD_data/Katie SGLT2/Plots/")

DPI=800
grViz(plot) %>%
  export_svg() %>%
  charToRaw %>%
  rsvg::rsvg(width = 6*DPI, height = 3*DPI) %>% tiff::writeTIFF("flowdiagram_contemporary.tiff", bits.per.sample = 8L)


  



