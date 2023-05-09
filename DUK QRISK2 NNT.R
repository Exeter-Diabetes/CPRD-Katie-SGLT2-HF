
# For poster, used uncalibrated values
# For presentation, used calibrated

############################################################################################

# Setup
library(tidyverse)
library(aurum)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "test-remote",cprdConf = "~/.aurum.yaml")

analysis = cprd$analysis("prev")

index_date <- as.Date("2020-02-01")


############################################################################################

final_merge <- final_merge %>% analysis$cached("final_merge")


analysis = cprd$analysis("katie_duk23")


duk_cohort <- final_merge %>%
  filter(diabetes_type=="type 2" &
           index_date_age>=18 &
           (
    datediff(index_date, pre_index_date_latest_dpp4)<=183 |
      datediff(index_date, pre_index_date_latest_glp1)<=183 | 
      datediff(index_date, pre_index_date_latest_mfn)<=183 | 
      datediff(index_date, pre_index_date_latest_sglt2)<=183 | 
      datediff(index_date, pre_index_date_latest_su)<=183 | 
      datediff(index_date, pre_index_date_latest_tzd)<=183 | 
      datediff(index_date, pre_index_date_latest_insulin)<=183
  ) &
    !is.na(qdiabeteshf_5yr_score) &
    !is.na(qrisk2_10yr_score)) %>%
  analysis$cached("cohort_ids", unique_indexes="patid")


duk_cohort <- duk_cohort %>%
  mutate(cvd=ifelse(pre_index_date_angina==1 |
                      pre_index_date_heartfailure==1 |
                      pre_index_date_ihd==1 |
                      pre_index_date_myocardialinfarction==1 |
                      pre_index_date_pad==1 |
                      pre_index_date_revasc==1 |
                      pre_index_date_stroke==1 |
                      pre_index_date_tia==1, 1L, 0L),
         
         nice_status=ifelse(cvd==1, "cvd",
                            ifelse(qrisk2_10yr_score>10, "high_risk", "low_risk_no_cvd")))
           
duk_cohort %>% count()
#269606

duk_cohort %>% group_by(nice_status) %>% summarise(count=n())
# 1 high_risk        154785
# 2 cvd               91719
# 3 low_risk_no_cvd   23102



duk_cohort <- duk_cohort %>%
  mutate(centred_qdiabeteshf_lin_predictor=ifelse(gender==1, qdiabeteshf_lin_predictor-0.464641698462435, qdiabeteshf_lin_predictor-0.468985206202883)) %>%
  mutate(qdiabeteshf_5yr_score_cal=ifelse(gender==1, (1-(0.965209400040619^exp(centred_qdiabeteshf_lin_predictor)))*100, (1-(0.965202142600253^exp(centred_qdiabeteshf_lin_predictor)))*100))


# NNT for qrisk2>10%:

duk_cohort %>%
  filter(nice_status=="high_risk") %>%
  summarise(mean=mean(qdiabeteshf_5yr_score_cal, na.rm=TRUE))
#mean qdhf=6.95
# mean benefit=2.57
# NNT=39



#[1] "Mean linear predictor male: 0.464641698462435"
#[1] "Mean linear predictor female: 0.468985206202883"
#[1] "New male surv: 0.965209400040619"
#[1] "New female surv: 0.965202142600253"




duk_cohort <- duk_cohort %>%
  mutate(status=ifelse(nice_status=="cvd", "cvd", "other"))


# For NNT of 26, mean benefit in group: 100/26=3.846%
# Mean QDHF: 3.846/0.37=10.395%

local_cohort <- collect(duk_cohort %>% filter(status=="other") %>% select(patid, qdiabeteshf_5yr_score_cal))

test <- local_cohort %>% filter(4.96<qdiabeteshf_5yr_score_cal)
summary(test)
#mean=10.395

#test:n=82,724


duk_cohort <- duk_cohort %>%
  mutate(status=ifelse(status!="other", status, ifelse(qdiabeteshf_5yr_score_cal>4.96, "nnt_26", "not_nnt_26"))) %>%
  analysis$cached("duk_cohort_cal2", unique_indexes="patid")


duk_cohort %>%
  group_by(status) %>%
  summarise(count=n())

# 
# 1 not_nnt_26  95163
# 2 nnt_26       82724
# 3 cvd          91719


# No. of HF cases for QRISK2:
## 154785/39 = 3969

## NNT 26:82724/26= 3181



## Not NNT group:

duk_cohort %>%
  group_by(status) %>%
  summarise(mean_qdhf=mean(qdiabeteshf_5yr_score_cal))

#2.61
# mean absolute benefit = 0.37*2.61 = 0.9657
#NNT = 100/0.9657= 104



