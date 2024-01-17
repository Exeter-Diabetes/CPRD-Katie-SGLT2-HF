# check overlap




test <- cohort %>%
  select(patid, studydrug, drugclass, dstartdate, dstopdate, hf_censdate) %>%
  group_by(patid) %>%
  arrange(patid, dstartdate) %>%
  mutate(overlap=ifelse(!is.na(lag(hf_censdate)) & dstartdate<lag(hf_censdate), 1, 0))
#141,123

table(test$studydrug, test$overlap)
#0     1
#DPP4SU  84531 10811
#SGLT2 45781     0


test <- cohort %>%
  select(patid, studydrug, drugclass, dstartdate, dstopdate, hf_censdate, DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, ncurrtx) %>%
  filter((drugclass=="DPP4" & GLP1==0 & SGLT2==0 & SU==0 & TZD==0 & INS==0 & MFN==1) | (drugclass=="SU" & DPP4==0 & GLP1==0 & SGLT2==0 & TZD==0 & INS==0 & MFN==1) | (drugclass=="SGLT2" & DPP4==0 & GLP1==0 & SU==0 & TZD==0 & INS==0 & MFN==1)) %>%
  group_by(patid) %>%
  arrange(patid, dstartdate) %>%
  mutate(overlap=ifelse(!is.na(lag(hf_censdate)) & dstartdate<lag(hf_censdate), 1, 0))
#79,096

table(test$studydrug)

table(test$studydrug, test$overlap)
#0     1
#DPP4SU  57703  2619
#SGLT2 18774     0

table(test$ncurrtx)
# all are 2
