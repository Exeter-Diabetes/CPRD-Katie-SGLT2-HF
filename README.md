# CPRD-Katie-CVD-HF

## Introduction

&nbsp;

## Risk score (QDiabetes-Heart Failure [QDHF] 2015 and QRISK2 2017) information

QDiabetes-Heart Failure paper: https://bmjopen.bmj.com/content/5/9/e008503
QRISK2 paper: https://www.bmj.com/content/336/7659/1475 
QRISK2 2013 codelists: https://qrisk.org/QRISK2-2013-Code-Lists.pdf (link no longer works)

&nbsp;

### Variable coding

The following table provides a brief overview of how we have defined the variables for these risk scores in CPRD, and how this might differ to the variable definition used to validate these risk scores / the way in which QRISK2 has previously been automatically run on GP systems in the UK (where known; information from papers and QRISK 2013 codelists as linked above). All our codelists are available here: https://github.com/Exeter-Diabetes/CPRD-Codelists/.

| Variable | Which risk scores | How this was defined during validation / for QRISK2 in automatic GP systems | How we have defined it |
| --- | --- | --- | --- |
| Age | Both |  | Age from DOB; [DOB definition](https://github.com/Exeter-Diabetes/CPRD-Codelists/blob/main/readme.md#general-notes-on-implementation) |
| Sex | Both |  | Provided by CPRD |
| Ethnicity | Both | UK census Read codes, no details on how conflicting codes dealt with | Medcode list (includes all of QRISK2 Read codes + extras); hospital (HES) ethnicity used if GP missing; and our algorithm if conflicting codes are present as per [our algorithm](https://github.com/Exeter-Diabetes/CPRD-Codelists/blob/main/readme.md#ethnicity). Some missingness |
| Deprivation | Both | Uses patient postcode + lookup of Townsend Deprivation Scores | We only have IMD deciles, not Townsend Scores, so are using the median Townsend Score for LSOAs with same IMD decile as patient |
| Smoking status | Both | Read codes, but quantity (value) associated with code used to determine whether non/current light/current moderate/current heavy smoker (ignore values for ex-smoker codes). If value missing, use code type to determine category (light if unspecified) | Medcode list (includes all of QRISK2 Read codes and extras) and have tried to replicate their algorithm (see [details](https://github.com/Exeter-Diabetes/CPRD-Codelists/blob/main/readme.md#smoking)). We only include values with valid/missing numunitid. Some missingness |
| Diabetes status (none/T1/T2; or just T1/T2 for QDHF) | Both | Type 1 and Type 2 (includes unspecified type) Read codes. Uses most recent code; if has both types on that day then Type 2 | [Our algorithms](https://github.com/Exeter-Diabetes/CPRD-Codelists/blob/main/readme.md#diabetes-algorithms) |
| Diabetes duration | QDHF |  | [Our algorithm](https://github.com/Exeter-Diabetes/CPRD-Codelists/blob/main/readme.md#defining-the-diagnosis-date-of-diabetes) for diagnosis date |
| Angina or heart attack in 1st degree relative <60 | QRISK2 | Read codes | Medcode list (includes all of QRISK2 Read codes – including ones for angina/MI in relative <age 55 which differs from <60 stated on tool; and extras including CVD in 1st degree relative <60 [not specifically MI or stroke]) |
| Heart attack, angina or stroke | QDHF |  | Medcode and ICD10 (HES) lists |
| CKD stage 4 or 5 | Both | Read codes | [Our algorithm](https://github.com/Exeter-Diabetes/CPRD-Codelists/blob/main/readme.md#ckd-chronic-kidney-disease-stage) which uses creatinine measurements and CKD stage 5 codes. If no baseline CKD stage, set to 0 not missing |
| Atrial fibrillation | Both | Read codes | Medcode list (doesn't include all of their AF Read codes, but theirs add <0.1% to occurrences) and ICD10 (HES) codes |
| On BP meds | QRISK2 | QRISK (original) paper says they used beta-blockers, ACE-inhibitors, calcium channel blockers and thiazide-like diuretics.<br>QRISK2 paper says this variable is ‘treated hypertension’ =  diagnosis of hypertension (Read code list available) and ≥1 current prescription of at ≥1 antihypertensive agent.<br>QRISK3 paper says ‘Use of drugs at baseline was defined as at least two prescriptions, with the most recent one no more than 28 days before the date of entry to the cohort’ | Prodcode lists for the 4 types in QRISK paper, didn’t use hypertension codes as per QRISK2, used rule that had to have ≥2 prescriptions with last 1 within 28 days of index date as per QRISK3 paper |
| Rheumatoid arthritis | QRISK2 |  | Medcode and ICD10 (HES) lists |
| These only identify ~half of the people as QRISK2 codelist. 2% of drug starts have RA at drug start currently |
| HbA1c | QDHF |  | Medcode list, only including if valid/missing numunitid. Some missingness. |
| Cholesterol:HDL ratio | Both | Separate cholesterol and HDL Read code lists | Medcode lists for HDL and total cholesterol, only including if valid/missing numunitid. HDL codelist doesn’t include one of their seven codes – adds <0.1% to total occurrences. Not using specific cholesterol:HDL ratio codes as these add <1% additional patients |
| SBP | Both | Read code list (one code only) | Medcode list (includes their one SBP Read code), only including if valid/missing numunitid |
| BMI | Both | Read code list (one code only) | Medcode list (includes their one BMI Read code), only including if valid/missing numunitid |

&nbsp;

Additional: QRISK3 paper states that for SBP, BMI, and smoking status, the most recent values before baseline were used. For cholesterol:HDL ratio, the most recent before or after baseline was used as long as it was before CVD diagnosis/censoring/statins. QDHF paper states that for continuous variables, the closest values prior to baseline or within 6 months after were used. We have used the most recent value prior to baseline (maximum of two years earlier, or 6 months for HbA1c) or up to 7 days after if closer.

&nbsp;

### Algorithm coding

Code for implementing the QDiabetes-Heart Failure 2015 and QRISK2 2017 algorithms was adapted from the open source software available on their respective websites (QDHF: https://qdiabetes.org/index.php; QRISK2 no longer available) and are included in the [EHRBiomarkr R package](https://github.com/Exeter-Diabetes/EHRBiomarkr).

&nbsp;

### Implementation

Those with biomarker values outside the range used by the online versions of the QDiabetes-Heart Failure 2015 and QRISK2 2017 calculators were excluded. The ranges were:

* QDHF:
    *	Cholesterol:HDL ratio: missing or 1-11
    * HbA1c: 40-150 mmol/mol
    * SBP: missing or 70-210 mmHg
    * Age: 25-84 years
    * Also excluded if BMI<20 kg/m2 as very different from development cohort
* QRISK2:
  * Cholesterol:HDL ratio: missing or 1-12
  * SBP: missing or 70-210 mmHg
  * Age: 25-84 years
  * Also exclude if BMI<20 kg/m2 as very different from development cohort

&nbsp;

## Outcome variable definitions
