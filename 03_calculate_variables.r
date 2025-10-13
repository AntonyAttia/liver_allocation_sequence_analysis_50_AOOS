# ==============================================================================
# 03_calculate_variables.r
# Calculate derived variables (DRI, DSRI, MELD, timing, etc.)
# ==============================================================================

library(tidyverse)
library(here)
library(lubridate)

# Load cohort from previous step
load(here("data_processed", "02_cohort_built.rda"))

cat("Starting variable calculations...\n")
cat("Input cohort size:", nrow(final_cohort), "\n\n")

# ==============================================================================
# 1. CENTERS CONTACTED (by distinct centers with specific N to offer, B/bypass excluded)
# ==============================================================================
cat("Calculating centers contacted...\n")

centres_per_run <- ptr_reduced %>%
  inner_join(final_cohort %>% select(DONOR_ID, MATCH_ID_CODE, PTR_SEQUENCE_NUM_tx),
             by = c("DONOR_ID", "MATCH_ID_CODE")) %>%
  filter(OFFER_ACCEPT == "N", PTR_SEQUENCE_NUM < PTR_SEQUENCE_NUM_tx) %>%
  select(DONOR_ID, MATCH_ID_CODE, WL_ID_CODE) %>%
  left_join(liver_reduced %>% distinct(WL_ID_CODE, LISTING_CTR_CODE),
            by = "WL_ID_CODE") %>%
  group_by(DONOR_ID, MATCH_ID_CODE) %>%
  summarise(
    n_unique_ctr_contacted = n_distinct(LISTING_CTR_CODE, na.rm = TRUE) + 1L,
    .groups = "drop"
  )

final_cohort <- final_cohort %>%
  left_join(centres_per_run, by = c("DONOR_ID", "MATCH_ID_CODE")) %>%
  mutate(n_unique_ctr_contacted = coalesce(n_unique_ctr_contacted, 1L))

# ==============================================================================
# 2. STANDARDIZE VARIABLE TYPES
# ==============================================================================
cat("Standardizing variable types...\n")

final_cohort <- final_cohort %>%
  mutate(
    # Segment flag from ORGAN_PLACED
    ORGAN_PLACED_uc = toupper(trimws(as.character(ORGAN_PLACED))),
    tx_split_flag_from_organ = case_when(
      ORGAN_PLACED_uc == "S" ~ 1L,
      ORGAN_PLACED_uc == "W" ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    # Numeric conversions
    AGE_DON_donor      = as.numeric(AGE_DON_donor),
    HGT_CM_DON_CALC    = as.numeric(HGT_CM_DON_CALC),
    COLD_ISCH          = as.numeric(COLD_ISCH),
    COLD_ISCH          = ifelse(!is.na(COLD_ISCH) & COLD_ISCH > 48, COLD_ISCH/60, COLD_ISCH),
    SHARE_TY           = as.numeric(SHARE_TY),
    BMI_DON_CALC_donor = as.numeric(BMI_DON_CALC_donor),
    SGOT_DON           = as.numeric(SGOT_DON),
    SGPT_DON           = as.numeric(SGPT_DON),
    SODIUM170_VAL_DON  = as.numeric(SODIUM170_VAL_DON),
    TBILI_DON_donor    = as.numeric(TBILI_DON_donor),
    
    # Character flags
    NON_HRT_DON_donor  = toupper(trimws(as.character(NON_HRT_DON_donor))),
    CDC_RISK_HIV_DON   = toupper(trimws(as.character(CDC_RISK_HIV_DON))),
    HBV_CORE_DON       = toupper(trimws(as.character(HBV_CORE_DON))),
    HEP_C_ANTI_DON     = toupper(trimws(as.character(HEP_C_ANTI_DON))),
    HIST_DIABETES_DON  = toupper(trimws(as.character(HIST_DIABETES_DON))),
    HIST_HYPERTENS_DON = toupper(trimws(as.character(HIST_HYPERTENS_DON))),
    GENDER_DON_donor   = toupper(trimws(as.character(GENDER_DON_donor)))
  )

# ==============================================================================
# 3. CALCULATE DRI (DONOR RISK INDEX)
# ==============================================================================
cat("Calculating DRI...\n")

final_cohort <- final_cohort %>%
  mutate(
    age_don    = AGE_DON_donor,
    eth_black  = ifelse(ETHCAT_DON == 2, 1L, 0L),
    eth_other  = ifelse(!is.na(ETHCAT_DON) & !(ETHCAT_DON %in% c(1,2)), 1L, 0L),
    cod_anoxia = as.integer(COD_CAD_DON_donor == 1),
    cod_cva    = as.integer(COD_CAD_DON_donor == 2),
    cod_other  = as.integer(!is.na(COD_CAD_DON_donor) & !(COD_CAD_DON_donor %in% c(1,2,3))),
    dcd_flag = ifelse(NON_HRT_DON_donor == "Y", 1L, 0L),    
    hgt_cm     = HGT_CM_DON_CALC,
    cit_hrs    = COLD_ISCH,
    share_ty   = SHARE_TY,
    tx_split_flag = ifelse(tx_split_flag_from_organ == 1L, 1L, 0L),
    
    DRI_base =
      case_when(age_don >= 40 & age_don < 50 ~ 0.154,
                age_don >= 50 & age_don < 60 ~ 0.274,
                age_don >= 60 & age_don < 70 ~ 0.424,
                age_don >= 70                ~ 0.501,
                .default = 0) +
      case_when(eth_black == 1 ~ 0.176,
                eth_other == 1 ~ 0.126,
                .default = 0) +
      case_when(cod_anoxia == 1 ~ 0.079,
                cod_cva    == 1 ~ 0.145,
                cod_other  == 1 ~ 0.184,
                TRUE ~ 0) +
      ifelse(dcd_flag == 1, 0.411, 0) +
      case_when(share_ty == 5 ~ 0.244,
                share_ty == 4 ~ 0.105,
                TRUE ~ 0) +
      ifelse(!is.na(cit_hrs), 0.010 * (cit_hrs - 8), 0) +
      ifelse(!is.na(hgt_cm),  0.066 * (170 - hgt_cm) * 0.1, 0) +
      ifelse(tx_split_flag == 1, 0.422, 0),
    
    DRI = exp(DRI_base)
  )

# ==============================================================================
# 4. CALCULATE DSRI (DISCARD RISK INDEX)
# ==============================================================================
cat("Calculating DSRI...\n")

.is_yes <- function(x) toupper(trimws(as.character(x))) %in% c("Y","YES","1","TRUE")
.is_pos <- function(x) toupper(trimws(as.character(x))) %in% c("P","POS","POSITIVE","1","Y","YES","TRUE")
.is_fem <- function(x) toupper(trimws(as.character(x))) %in% c("F","FEMALE")

final_cohort <- final_cohort %>%
  mutate(
    DSRI_base =
      case_when(.is_yes(NON_HRT_DON_donor) ~ 2.649, TRUE ~ 0) +
      case_when(!is.na(AGE_DON_donor) & AGE_DON_donor < 15 ~ -0.678,
                AGE_DON_donor >=15 & AGE_DON_donor < 20 ~ -0.995,
                AGE_DON_donor >=20 & AGE_DON_donor < 25 ~ -0.876,
                AGE_DON_donor >=25 & AGE_DON_donor < 30 ~ -0.576,
                AGE_DON_donor >=45 & AGE_DON_donor < 55 ~  0.454,
                AGE_DON_donor >=55 & AGE_DON_donor < 60 ~  0.612,
                AGE_DON_donor >=60 & AGE_DON_donor < 65 ~  0.699,
                AGE_DON_donor >=65 & AGE_DON_donor < 70 ~  0.755,
                AGE_DON_donor >=70 & AGE_DON_donor < 75 ~  0.793,
                AGE_DON_donor >=75                      ~  0.843,
                TRUE ~ 0) +
      case_when(!is.na(BMI_DON_CALC_donor) & BMI_DON_CALC_donor <= 20 ~ 0.162,
                BMI_DON_CALC_donor >= 30 & BMI_DON_CALC_donor < 35    ~ 0.371,
                BMI_DON_CALC_donor >= 35 & BMI_DON_CALC_donor < 40    ~ 0.666,
                BMI_DON_CALC_donor >= 40                               ~ 1.222,
                TRUE ~ 0) +
      case_when(.is_yes(CDC_RISK_HIV_DON) | .is_pos(CDC_RISK_HIV_DON) ~ -0.137, TRUE ~ 0) +
      case_when(COD_CAD_DON_donor == 2 ~ 0.160, TRUE ~ 0) +
      case_when(ETHCAT_DON == 2 ~ -0.513, TRUE ~ 0) +
      case_when(.is_fem(GENDER_DON_donor) ~ -0.085, TRUE ~ 0) +
      case_when(.is_pos(HBV_CORE_DON) ~ 0.169, TRUE ~ 0) +
      case_when(.is_pos(HEP_C_ANTI_DON) ~ 1.235, TRUE ~ 0) +
      case_when(.is_yes(HIST_DIABETES_DON) ~ 0.316, TRUE ~ 0) +
      case_when(.is_yes(HIST_HYPERTENS_DON) ~ 0.160, TRUE ~ 0) +
      case_when(!is.na(SGOT_DON) & SGOT_DON >= 100 & SGOT_DON < 200 ~ 0.208,
                !is.na(SGOT_DON) & SGOT_DON >= 200 & SGOT_DON < 500 ~ 0.487,
                !is.na(SGOT_DON) & SGOT_DON >= 500                  ~ 1.585,
                TRUE ~ 0) +
      case_when(!is.na(SGPT_DON) & SGPT_DON >=  50 & SGPT_DON < 100 ~ 0.255,
                !is.na(SGPT_DON) & SGPT_DON >= 100 & SGPT_DON < 200 ~ 0.208,
                !is.na(SGPT_DON) & SGPT_DON >= 200 & SGPT_DON < 500 ~ 0.485,
                !is.na(SGPT_DON) & SGPT_DON >= 500                  ~ 1.111,
                TRUE ~ 0) +
      case_when(!is.na(SODIUM170_VAL_DON) & SODIUM170_VAL_DON > 160 ~ 0.204, TRUE ~ 0) +
      case_when(!is.na(TBILI_DON_donor) & TBILI_DON_donor >=  1 & TBILI_DON_donor <  2 ~ 0.301,
                !is.na(TBILI_DON_donor) & TBILI_DON_donor >=  2 & TBILI_DON_donor <  3 ~ 0.797,
                !is.na(TBILI_DON_donor) & TBILI_DON_donor >=  3 & TBILI_DON_donor <  4 ~ 1.001,
                !is.na(TBILI_DON_donor) & TBILI_DON_donor >=  4 & TBILI_DON_donor <  5 ~ 1.420,
                !is.na(TBILI_DON_donor) & TBILI_DON_donor >=  5 & TBILI_DON_donor < 10 ~ 2.023,
                !is.na(TBILI_DON_donor) & TBILI_DON_donor >= 10                        ~ 3.228,
                TRUE ~ 0),
    
    DSRI = exp(DSRI_base)
  )

# ==============================================================================
# 5. CALCULATE BSA RATIO (BODY SURFACE AREA)
# ==============================================================================
cat("Calculating BSA ratio...\n")

calc_BSA <- function(weight_kg, height_cm) {
  sqrt((as.numeric(weight_kg) * as.numeric(height_cm)) / 3600)
}

final_cohort <- final_cohort %>%
  mutate(
    recip_wt_kg = as.numeric(WGT_KG_CALC),
    recip_ht_cm = as.numeric(HGT_CM_CALC),
    BSA_recipient_m2 = calc_BSA(recip_wt_kg, recip_ht_cm),
    
    donor_wt_kg = as.numeric(WGT_KG_DON_CALC),
    donor_ht_cm = as.numeric(HGT_CM_DON_CALC),
    BSA_donor_m2 = calc_BSA(donor_wt_kg, donor_ht_cm),
    
    BSA_ratio_donor_to_recipient = BSA_donor_m2 / BSA_recipient_m2,
    
    COLD_ISCH_num = as.numeric(COLD_ISCH),
    WAIT_DAYS = as.numeric(DAYSWAIT_CHRON)
  )

# ==============================================================================
# 6. CALCULATE ADJUSTED MELD SCORES including from lab values
# ==============================================================================
cat("Calculating adjusted MELD scores...\n")

final_cohort <- final_cohort %>%
  mutate(
    Na_tx   = as.numeric(FINAL_SERUM_SODIUM),
    Tb_tx   = as.numeric(TBILI_TX),
    INR_tx  = as.numeric(INR_TX),
    Alb_tx  = as.numeric(ALBUMIN_TX),
    Cr_tx   = as.numeric(CREAT_TX),
    
    # Clamp inputs
    Na_b    = pmin(pmax(Na_tx, 125), 137),
    Tb_b    = pmax(Tb_tx, 1),
    INR_b   = pmax(INR_tx, 1),
    Alb_b   = pmin(pmax(Alb_tx, 1.5), 3.5),
    Cr1_b   = pmin(pmax(Cr_tx, 1.0), 4.0),
    Cr2_b   = pmin(pmax(Cr_tx, 1.0), 3.0),
    
    fem_133 = if_else(toupper(trimws(as.character(GENDER))) %in% c("F","FEMALE"), 1.33, 0),
    
    # Raw scores
    MELD_raw    = 9.57*log(Cr1_b) + 3.78*log(Tb_b) + 11.20*log(INR_b) + 6.43,
    MELDNa_raw  = MELD_raw - 0.033*MELD_raw*(137 - Na_b) + 1.32*(137 - Na_b),
    MELD3_raw   = fem_133 + 4.56*log(Tb_b) + 0.82*(137 - Na_b) -
      0.24*(137 - Na_b)*log(Tb_b) + 9.09*log(INR_b) +
      11.14*log(Cr2_b) + 1.85*(3.5 - Alb_b) -
      1.83*(3.5 - Alb_b)*log(Cr2_b) + 6,
    
    # Rounded and capped
    MELD_adj      = pmin(40L, floor(MELD_raw   + 0.5)),
    MELDNa_adj    = pmin(40L, floor(MELDNa_raw + 0.5)),
    MELD3_adj     = pmin(40L, floor(MELD3_raw  + 0.5)),
    
    # Unified field
    MELD_PELD_ADJ = case_when(
      !is.na(AGE) & AGE >= 12 ~ MELD3_adj,
      TRUE                    ~ MELD_PELD_LAB_SCORE
    )
  )

# ==============================================================================
# 7. CALCULATE CLAMP TIMING VARIABLES
# ==============================================================================
cat("Calculating clamp timing variables...\n")

final_cohort <- final_cohort %>%
  mutate(
    CLAMP_DT = ymd_hms(paste(CLAMP_DATE, CLAMP_TIME), tz = "UTC", quiet = TRUE),
    diff_initial_clamp = as.numeric(difftime(INITIAL_RESPONSE_DT_tx, CLAMP_DT, units = "hours")),
    diff_final_clamp   = as.numeric(difftime(FINAL_RESPONSE_DT_tx, CLAMP_DT, units = "hours"))
  )

# ==============================================================================
# 8. CALCULATE RESPONSE TIMING FROM EARLIEST INITIAL RESPONSE IN MATCH-RUN
# ==============================================================================
cat("Calculating response timing variables...\n")

run_start <- ptr_reduced %>%
  filter(!is.na(DONOR_ID), !is.na(MATCH_ID_CODE)) %>%
  group_by(DONOR_ID, MATCH_ID_CODE) %>%
  summarise(
    earliest_init = if (all(is.na(INITIAL_RESPONSE_DT))) {
      as.POSIXct(NA)
    } else {
      min(INITIAL_RESPONSE_DT, na.rm = TRUE)
    },
    .groups = "drop"
  )

final_cohort <- final_cohort %>%
  left_join(run_start, by = c("DONOR_ID", "MATCH_ID_CODE")) %>%
  mutate(
    diff_initial_hours = case_when(
      is.na(earliest_init) ~ NA_real_,
      is.na(INITIAL_RESPONSE_DT_tx) ~ NA_real_,
      TRUE ~ as.numeric(difftime(INITIAL_RESPONSE_DT_tx, earliest_init, units = "hours"))
    ),
    diff_final_hours = case_when(
      is.na(earliest_init) ~ NA_real_,
      is.na(FINAL_RESPONSE_DT_tx) ~ NA_real_,
      TRUE ~ as.numeric(difftime(FINAL_RESPONSE_DT_tx, earliest_init, units = "hours"))
    )
  )

# ==============================================================================
# 9. CREATE CATEGORICAL VARIABLES
# ==============================================================================
cat("Creating categorical variables...\n")

final_cohort <- final_cohort %>%
  mutate(
    # Donation type
    donation_type = case_when(
      NON_HRT_DON_donor == "Y" ~ "DCD",
      NON_HRT_DON_donor == "N" ~ "DBD",
      TRUE ~ NA_character_
    ),
    
    # Region grouping
    region_type4 = case_when(
      REGION %in% c(1,2,9) ~ "Northeast",
      REGION %in% c(3,4,11) ~ "Southeast",
      REGION %in% c(7,8,10) ~ "Midwest",
      REGION %in% c(5,6) ~ "West",
      TRUE ~ NA_character_
    ),
    
    # Definition 1: Machine perfusion as-is (original variable)
    LI_MACHINE_PERFUSION_original = factor(
      LI_MACHINE_PERFUSION, 
      levels = c("Y", "N")
    ),
    
    # Definition 2: Post-FDA approval period (first of TMDX and Organox available)
    LI_MACHINE_PERFUSION_FDA_era = case_when(
      LI_MACHINE_PERFUSION == "Y" ~ "Yes",
      COLD_ISCH_num > 12 & TX_DATE >= as.Date("2021-09-29") ~ "Yes",
      LI_MACHINE_PERFUSION == "N" ~ "No",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Yes", "No")),
    
    # Definition 3: Cold ischemia >12
    LI_MACHINE_PERFUSION_cold_isch = case_when(
      LI_MACHINE_PERFUSION == "Y" ~ "Yes",
      COLD_ISCH_num > 12 ~ "Yes", 
      LI_MACHINE_PERFUSION == "N" ~ "No",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Yes", "No"))
  )
# ==============================================================================
# SAVE FINAL COHORT WITH ALL VARIABLES
# ==============================================================================
save(final_cohort, 
     file = here("data_processed", "03_cohort_with_variables.rda"))

cat("\n=== VARIABLE CALCULATION COMPLETE ===\n")
cat("Final cohort size:", nrow(final_cohort), "\n")
cat("Total variables:", ncol(final_cohort), "\n")
cat("Saved to: data_processed/03_cohort_with_variables.rda\n")