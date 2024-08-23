
############################################################################################

# Setup
library(tidyverse)
library(ggthemes)
library(gtools)
library(gridExtra)
library(table1)
library(ComplexUpset)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())


############################################################################################

# 1 Cohort selection - see script 00

setwd("/slade/CPRD_data/Katie SGLT2/Processed data/")
load(file="contemporary_cohort_jun24.rda")
contemporary_cohort <- cohort
#254499

load("treatment_outcome_cohort_jun24.rda")
#169,041


############################################################################################

# 2. Histogram of predicted benefits

both_cohorts <- contemporary_cohort %>%
  select(qdiabeteshf_5yr_score) %>%
  mutate(Cohort="Contemporary cohort") %>%
  rbind((cohort %>%  select(qdiabeteshf_5yr_score) %>% mutate(Cohort="Main study cohort (SGLT2i/DPPi/SU treated)"))) %>%
  mutate(qdhf_survival=(100-qdiabeteshf_5yr_score)/100,
         qdhf_survival_sglt2=qdhf_survival^0.63,
         qdiabeteshf_5yr_score_sglt2=100-(qdhf_survival_sglt2*100),
         qdhf_sglt2_benefit=qdhf_survival_sglt2-qdhf_survival) %>%
  group_by(Cohort) %>%
  mutate(weights = 1/n()) %>%
  ungroup()


tiff("/slade/CPRD_data/Katie SGLT2/Plots/benefit_histogram_cohorts.tiff", width=8, height=5, units = "in", res=400) 

ggplot(both_cohorts, aes(x=qdhf_sglt2_benefit*100, color=Cohort, fill=Cohort)) + 
  geom_histogram(aes(weight=weights*100), alpha = 0.6, position="identity", color=NA, binwidth=0.05) +
  ylab("Proportion of cohort (%)") + xlab("Predicted 5-year SGLT2i benefit (%)") +
  scale_x_continuous(limits=c(0, 6),breaks=c(seq(0, 6,by=1))) +
  scale_y_continuous(limits=c(0,3.5),breaks=c(seq(0, 3.5, by=0.5))) +
  theme_base() + 
  theme(plot.background = element_blank(),
        legend.position=c(0.65, 0.85)) +
  scale_fill_manual(values=c(rgb(0,60/255,140/255), rgb(240/255,180/255,0)))

dev.off()

############################################################################################

# 3. Prep for NNTs in contemporary cohort

rm(list=ls())

setwd("/slade/CPRD_data/Katie SGLT2/Processed data/")
load(file="contemporary_cohort_jun24.rda")

cohort <- cohort %>%
  mutate(qdhf_survival=(100-qdiabeteshf_5yr_score)/100,
         qdhf_survival_sglt2=qdhf_survival^0.63,
         qdiabeteshf_5yr_score_sglt2=100-(qdhf_survival_sglt2*100),
         qdhf_sglt2_benefit=qdhf_survival_sglt2-qdhf_survival)


## Add strata used previously
cohort <- cohort %>%
  mutate(qdhf_benefit_strata=cut(qdhf_sglt2_benefit*100, breaks = c(-Inf, 0.5, 1, 2, 3, Inf), labels = c(1, 2, 3, 4, 5)))

prop.table(table(cohort$qdhf_benefit_strata))
#18 26  31  15  11


## Add categories based on strata
cohort <- cohort %>%
  mutate(all=1,
         abs_hf_benefit_over_3=ifelse(qdhf_benefit_strata==5, 1, 0),
         abs_hf_benefit_over_2=ifelse(qdhf_benefit_strata==4 | qdhf_benefit_strata==5, 1, 0),
         abs_hf_benefit_over_1=ifelse(qdhf_benefit_strata==3 | qdhf_benefit_strata==4 | qdhf_benefit_strata==5, 1, 0),
         abs_hf_benefit_over_0.5=ifelse(qdhf_benefit_strata!=1, 1, 0))


# Add in ADA/EASD and NICE categories

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
         hypertension=ifelse(presbp>=140 | pre_index_date_hypertension==1 | (!is.na(predbp) & predbp>90), 1, 0),
         smoking=ifelse(qrisk2_smoking_cat>1, 1, 0),
         dyslipidemia=ifelse((!is.na(pretotalcholesterol) & pretotalcholesterol>5) | statins==1, 1, 0),
         
         uacr=ifelse(!is.na(preacr), preacr, preacr_from_separate),
         albuminuria=ifelse(!is.na(uacr) & uacr>3, 1, 0),
         
         ada_easd=ifelse(index_date_age>=55 & (obesity+hypertension+smoking+dyslipidemia+albuminuria)>=2, 1, 0))


# NICE
cohort <- cohort %>%
  mutate(nice_qrisk2_10=ifelse(qrisk2_10yr_score>10, 1, 0),
         nice_qrisk2_20=ifelse(qrisk2_10yr_score>20, 1, 0),
         nice_qrisk2_30=ifelse(qrisk2_10yr_score>30, 1, 0))


############################################################################################

# 4. Calculate NNTs


categories <- c("all", "abs_hf_benefit_over_0.5", "abs_hf_benefit_over_1", "abs_hf_benefit_over_2", "abs_hf_benefit_over_3", "ada_easd", "nice_qrisk2_10", "nice_qrisk2_20", "nice_qrisk2_30")


data <- data.frame(proportion_treated=NA, median_with_iqr=NA, nnt_with_iqr=NA, strategy=NA)

for (i in categories) {

strategy <- cohort %>%
  filter(cohort[[i]]==1) %>%
  summarise(proportion_treated=paste0(round_pad((n()*100)/254499, 1), "%"),
            
            median=round(median(qdhf_sglt2_benefit, na.rm=T)*100,1),
            lq_median=round(quantile(qdhf_sglt2_benefit, prob=0.25)*100,1),
            uq_median=round(quantile(qdhf_sglt2_benefit, prob=0.75)*100,1),
              
            nnt_hf=round(1/median(qdhf_sglt2_benefit, na.rm=T)),
            lq_nnt=round(1/quantile(qdhf_sglt2_benefit, prob=0.75)),
            uq_nnt=round(1/quantile(qdhf_sglt2_benefit, prob=0.25)),
          
            median_with_iqr=paste0(median, " (", lq_median, ", ", uq_median, ")"),
            nnt_with_iqr=paste0(nnt_hf, " (", lq_nnt, ", ", uq_nnt, ")")) %>%
  select(proportion_treated, median_with_iqr, nnt_with_iqr) %>%
  mutate(strategy=i)

data <- rbind(data, strategy)

}


#####################################################################################################################

# Upset plots
## Only difference between three plots is label vjust on QRISK2 >20% and 'descending' sort intersections on QRISK2 20%


# Redo proportions so out of those with HF benefit>1

cohort %>% filter(abs_hf_benefit_over_1==1) %>% count()
#143187


data <- cohort %>%
  select(patid, ada_easd, abs_hf_benefit_over_1, nice_qrisk2_10, nice_qrisk2_20) %>%
  rename(`ADA-EASD`=ada_easd,
         `Absolute 5-year\nHF benefit >1%`=abs_hf_benefit_over_1,
         `QRISK2 >10%`=nice_qrisk2_10,
         `QRISK2 >20%`=nice_qrisk2_20)



intersect_columns <- c("QRISK2 >10%", "Absolute 5-year\nHF benefit >1%")

setwd("/slade/CPRD_data/Katie SGLT2/Processed data/Plots/")
tiff("upset_qrisk2_10.tiff", width=11, height=7, units = "in", res=800)

print(upset(data, intersect_columns, name="SGLT2i targeting strategy",
            min_degree=1,
            sort_intersections=F,
            sort_sets=F,
            intersections=list(
              c("QRISK2 >10%", "Absolute 5-year\nHF benefit >1%"),
              "Absolute 5-year\nHF benefit >1%",
              "QRISK2 >10%"
            ),
            base_annotations = list(
              'Intersection size'=(
                intersection_size(
                  text_mapping=aes(
                    label=paste0(round(!!get_size_mode('exclusive_intersection')/143187 * 100, 1), '%')
                  ), text=list(size=6, fontface="bold", vjust=c(0,0,-1,-1))
                )
                + ylab('Intersection as % of those with\nabsolute 5-year HF benefit >1%')
                + scale_y_continuous(
                  labels=scales::percent_format(scale=100 / 143187),
                  limits=c(0,100)/ 100 * 143187,
                  breaks=c(0, 25, 50, 75, 100) / 100 * 143187
                )
              )),
            set_sizes=upset_set_size(
              geom=geom_bar(stat='count'),
              mapping=aes(y=..count../nrow(data))
            )
            + ylab('Set size as % of whole cohort')
            + scale_y_reverse(labels=scales::percent,
                              limits=c(1,0))
            + theme(axis.text.x=element_text(size=15),
                    axis.title.x=element_text(size=18, vjust=-1.5)),
            themes=upset_default_themes(text=element_text(size=20),
                                        plot.margin = unit(c(0,0,0.5,1), "cm")),
            width_ratio=0.3,
            height_ratio=0.25
)) #,sort_intersections='ascending',

dev.off()



intersect_columns <- c("QRISK2 >20%", "Absolute 5-year\nHF benefit >1%")

setwd("/slade/CPRD_data/Katie SGLT2/Processed data/Plots/")
tiff("upset_qrisk2_20.tiff", width=11, height=7, units = "in", res=800)

print(upset(data, intersect_columns, name="SGLT2i targeting strategy",
            min_degree=1,
            sort_intersections=F,
            sort_sets=F,
            intersections=list(
              c("QRISK2 >20%", "Absolute 5-year\nHF benefit >1%"),
              "Absolute 5-year\nHF benefit >1%",
              "QRISK2 >20%"
            ),
            base_annotations = list(
              'Intersection size'=(
                intersection_size(
                  text_mapping=aes(
                    label=paste0(round(!!get_size_mode('exclusive_intersection')/143187 * 100, 1), '%')
                  ), text=list(size=6, fontface="bold", vjust=c(0,-1,0,-1))
                )
                + ylab('Intersection as % of those with\nabsolute 5-year HF benefit >1%')
                + scale_y_continuous(
                  labels=scales::percent_format(scale=100 / 143187),
                  limits=c(0,100)/ 100 * 143187,
                  breaks=c(0, 25, 50, 75, 100) / 100 * 143187
                )
              )),
            set_sizes=upset_set_size(
              geom=geom_bar(stat='count'),
              mapping=aes(y=..count../nrow(data))
            )
            + ylab('Set size as % of whole cohort')
            + scale_y_reverse(labels=scales::percent,
                              limits=c(1,0))
            + theme(axis.text.x=element_text(size=15),
                    axis.title.x=element_text(size=18, vjust=-1.5)),
            themes=upset_default_themes(text=element_text(size=20),
                                        plot.margin = unit(c(0,0,0.5,1), "cm")),
            width_ratio=0.3,
            height_ratio=0.25
))

dev.off()  



intersect_columns <- c("ADA-EASD", "Absolute 5-year\nHF benefit >1%")

setwd("/slade/CPRD_data/Katie SGLT2/Processed data/Plots/")
tiff("upset_ada_easd.tiff", width=11, height=7, units = "in", res=800)

print(upset(data, intersect_columns, name="SGLT2i targeting strategy",
            min_degree=1,
            sort_intersections=F,
            sort_sets=F,
            intersections=list(
              c("ADA-EASD", "Absolute 5-year\nHF benefit >1%"),
              "Absolute 5-year\nHF benefit >1%",
              "ADA-EASD"
            ),
            base_annotations = list(
              'Intersection size'=(
                intersection_size(
                  text_mapping=aes(
                    label=paste0(round(!!get_size_mode('exclusive_intersection')/143187 * 100, 1), '%')
                  ), text=list(size=6, fontface="bold", vjust=c(0,0,-1,-1))
                )
                + ylab('Intersection as % of those with\nabsolute 5-year HF benefit >1%')
                + scale_y_continuous(
                  labels=scales::percent_format(scale=100 / 143187),
                  limits=c(0,100)/ 100 * 143187,
                  breaks=c(0, 25, 50, 75, 100) / 100 * 143187
                )
              )),
            set_sizes=upset_set_size(
              geom=geom_bar(stat='count'),
              mapping=aes(y=..count../nrow(data))
            )
            + ylab('Set size as % of whole cohort')
            + scale_y_reverse(labels=scales::percent,
                              limits=c(1,0))
            + theme(axis.text.x=element_text(size=15),
                    axis.title.x=element_text(size=18, vjust=-1.5)),
            themes=upset_default_themes(text=element_text(size=20),
                                        plot.margin = unit(c(0,0,0.5,1), "cm")),
            width_ratio=0.3,
            height_ratio=0.25
)) #,sort_intersections='ascending',

dev.off()  

