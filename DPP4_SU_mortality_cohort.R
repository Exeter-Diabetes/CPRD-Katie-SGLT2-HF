
# Defining cohort for all-cause mortality work

library(tidyverse)
load("/slade/CPRD_data/mastermind_2022/20231205_t2d_1stinstance.Rda")

# Inclusion/exclusion criteria:
## a) T2Ds
## b) With HES linkage
## c) 1st instance
## d) Exclude if start drug within 90 days of registration

## e) Aged 18+
## f) DPP4/SU
## g) Initiated between 01/01/2013 and end of data (31/10/2020)
## h) Second line only (and not including where not sure of drugline)
## i) On metformin only at initiation
## j) No CVD (I'm using a broad definition: angina, IHD, MI, PAD, revasc, stroke, TIA [as per NICE but with TIA])
## k) No HF before index date
## l) No CKD (stage 3a-5) before index date (assume coded on all, so if missing assume negative)

# Use "t2d_1stinstance" cohort dataset which already has a)-d) applied

# e)-i) above
cohort <- t2d_1stinstance %>%
  filter(dstartdate_age>=18 &
           (drugclass=="DPP4" | drugclass=="SU") &
           dstartdate>=as.Date("2013-01-01") &
           drugline==2 &
           MFN==1 & Glinide==0 & GLP1==0 & INS==0 & SGLT2==0 & TZD==0 & ((DPP4==1 & SU==0) | (DPP4==0 & SU==1)))


# Remove if CVD before index date (j above)
cohort <- cohort %>%
  mutate(predrug_cvd=ifelse(predrug_angina==1 | predrug_ihd==1 | predrug_myocardialinfarction==1 | predrug_pad==1 | predrug_revasc==1 | predrug_stroke==1 | predrug_tia==1, 1, 0)) %>%
  filter(predrug_cvd==0)


# Remove if HF before index date (k above)
cohort <- cohort %>%
  filter(predrug_heartfailure==0)


# Remove if CKD before index date (l above)
cohort <- cohort %>%
  filter(is.na(preckdstage) | (preckdstage!="stage_3a" & preckdstage!="stage_3b" & preckdstage!="stage_4" & preckdstage!="stage_5"))

