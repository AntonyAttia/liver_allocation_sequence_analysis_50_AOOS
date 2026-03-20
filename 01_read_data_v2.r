# 01_read_data_v2.r
# Loading raw data files into R
# v2: Study period starts at EP implementation (2021-03-25).

library(tidyverse)
library(haven)

# Where to save/load .rda files. Edit this to your project's data_processed folder.
DATA_DIR <- "/Users/antonyattia/Desktop/Stanford/Open Offers/Manuscript/liver_allocation_sequence_analysis_50_AOOS/data_processed"
if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)

# Study period lower bound = EP policy implementation date
date_lb <- as.Date("2021-03-25")

# Read raw data files
ptr <- read_sas("~/Downloads/ptr_hr_li/ptr_li.sas7bdat")
liver <- read_sas("~/Downloads/SAS Dataset 202512/Liver/liver_data.sas7bdat")
donor <- read_sas("~/Downloads/SAS Dataset 202512/Deceased Donor/deceased_donor_data.sas7bdat")
labpanels <- read_sas("~/Downloads/labpanels.sas7bdat")

# Reduce to required columns only
ptr_reduced <- ptr %>%
  select(
    DONOR_ID, MATCH_DATE, MATCH_ID_CODE, MATCH_OPO, WL_ID_CODE, OFFER_ACCEPT,
    ORGAN_PLACED, INITIAL_RESPONSE_DT, FINAL_RESPONSE_DT,
    PTR_SEQUENCE_NUM, PRIME_OPO_REFUSAL_ID, OPO_OSTXT_PRIMARY
  )

liver_reduced <- liver %>%
  select(
    TX_DATE, WL_ID_CODE, DONOR_ID,
    DON_TY, MULTIORG, AGE_GROUP,
    LISTING_CTR_CODE, AGE_DON, GENDER_DON, ETHCAT_DON, BMI_DON_CALC,
    COD_CAD_DON, NON_HRT_DON, ALCOHOL_HEAVY_DON, MACRO_FAT_LI_DON,
    SGOT_DON, SGPT_DON, TBILI_DON, HCV_NAT_DON,
    AGE, GENDER, ETHCAT, BMI_CALC, MELD_PELD_LAB_SCORE, FINAL_SERUM_SODIUM,
    TBILI_TX, INR_TX, ALBUMIN_TX, CREAT_TX, ASCITES_TX, ENCEPH_TX, EXC_AT_TX,
    DIAB, HCV_SEROSTATUS, DISTANCE, COLD_ISCH, SHARE_TY, DAYSWAIT_CHRON,
    WGT_KG_CALC, HGT_CM_CALC, MED_COND_TRR, BILIARY, DIAG, DGN_TCR, DGN2_TCR,
    PSTATUS, PTIME, GSTATUS, GTIME, REGION, TRR_ID_CODE
  )

donor_reduced <- donor %>%
  select(
    ABO_DON, AGE_DON, ALCOHOL_HEAVY_DON, BMI_DON_CALC, CLIN_INFECT_DON,
    COD_CAD_DON, CREAT_DON, DEATH_CIRCUM_DON, DEATH_MECH_DON, DONOR_ID,
    GENDER_DON, HCV_NAT, HEP_C_ANTI_DON, HGT_CM_DON_CALC, HISTORYOFDRUGUSE,
    LI_BIOPSY_DON, LI_DISCARD_CD, LI_FIBROSIS, MACRO_FAT_LI_DON, NON_HRT_DON, NUM_LI_RECOV,
    OPO_CTR_CODE, TBILI_DON, WGT_KG_DON_CALC,
    BUN_DON, CDC_RISK_HIV_DON, CONSENT_DON_CARD, HBV_CORE_DON,
    HIST_DIABETES_DON, HIST_HYPERTENS_DON, HIST_CANCER_DON, HIST_CAD,
    SODIUM170_VAL_DON, VESSELS_NUM_STEN_DON, BLOOD_INF_CONF_DON,
    BLOOD_INF_DON, LI_MACHINE_PERFUSION, WARM_ISCH_TM_DON,
    CLAMP_DATE, CLAMP_TIME
  )

save(ptr_reduced, liver_reduced, donor_reduced, labpanels, date_lb,
     file = file.path(DATA_DIR, "01_raw_data_v2.rda"))

cat("01_read_data_v2.r complete\n")
cat("  Saved to:", file.path(DATA_DIR, "01_raw_data_v2.rda"), "\n")
cat("  Study period starts:", as.character(date_lb), "(EP implementation)\n")
cat("  ptr_reduced:", nrow(ptr_reduced), "rows\n")
cat("  liver_reduced:", nrow(liver_reduced), "rows\n")
cat("  donor_reduced:", nrow(donor_reduced), "rows\n")
cat("  labpanels:", nrow(labpanels), "rows\n")
