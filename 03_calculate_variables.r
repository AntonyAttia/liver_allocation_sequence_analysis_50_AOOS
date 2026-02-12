# 03_calculate_variables.r
# Calculate derived variables (DRI, DSRI, MELD, timing, etc.)

library(tidyverse)
library(here)
library(lubridate)

# Load cohort from previous step
load(here("data_processed", "02_cohort_built.rda"))
# Load raw data for refusal statistics calculation
load(here("data_processed", "01_raw_data.rda"))

cat("Input cohort size:", nrow(final_cohort), "\n\n")

# 1. CENTERS CONTACTED (by distinct centers with specific N to offer, B/bypass excluded)
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

# 2. STANDARDIZE VARIABLE TYPES
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

# 3. CALCULATE DRI (DONOR RISK INDEX)
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

# 4. CALCULATE DSRI (DISCARD RISK INDEX)
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

# 5. CALCULATE BSA RATIO (BODY SURFACE AREA)
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

# 6. CALCULATE ADJUSTED MELD SCORES including from lab values
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

# 7. CALCULATE CLAMP TIMING VARIABLES
final_cohort <- final_cohort %>%
  mutate(
    CLAMP_DT = ymd_hms(paste(CLAMP_DATE, CLAMP_TIME), tz = "UTC", quiet = TRUE),
    diff_initial_clamp = as.numeric(difftime(INITIAL_RESPONSE_DT_tx, CLAMP_DT, units = "hours")),
    diff_final_clamp   = as.numeric(difftime(FINAL_RESPONSE_DT_tx, CLAMP_DT, units = "hours"))
  )

# 8. CALCULATE RESPONSE TIMING FROM EARLIEST INITIAL RESPONSE IN MATCH-RUN
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

# 9. CREATE CATEGORICAL VARIABLES
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

# CALCULATE REFUSAL STATISTICS PER DONOR-RUN
# Open offer definition constants (matching 02_build_cohort.r)
OPEN_OFFER_TEXT_CODES <- c(798L, 799L, 898L)
OPEN_OFFER_TEXT_PATTERN <- regex("expedited|open offer|aggressive", ignore_case = TRUE)

# Create allocation type grouping
cohort_compare <- final_cohort %>%
  mutate(
    allocation_type = if_else(is_open_offer, "AOOS", "Sequential"),
    allocation_type = factor(allocation_type, levels = c("Sequential", "AOOS"))
  )

# Get allocation type for each donor-run (already filtered by date_lb via final_cohort)
cohort_keys <- cohort_compare %>%
  select(DONOR_ID, MATCH_ID_CODE, allocation_type) %>%
  distinct()

# Get center codes for patients
listing_center_codes <- liver_reduced %>%
  distinct(WL_ID_CODE, LISTING_CTR_CODE) %>%
  filter(!is.na(LISTING_CTR_CODE), !is.na(WL_ID_CODE))

# Prepare refusal data - filter to N refusals only
# Handle FINAL_RESPONSE_DT: if already POSIXct, use as-is; otherwise parse
refusals <- ptr_reduced %>%
  filter(OFFER_ACCEPT == "N", !is.na(FINAL_RESPONSE_DT)) %>%
  mutate(
    FINAL_RESPONSE_DT = case_when(
      inherits(FINAL_RESPONSE_DT, "POSIXct") ~ FINAL_RESPONSE_DT,
      inherits(FINAL_RESPONSE_DT, "Date") ~ as.POSIXct(FINAL_RESPONSE_DT),
      is.character(FINAL_RESPONSE_DT) ~ parse_date_time(FINAL_RESPONSE_DT, 
                                                          orders = c("ymd HMS", "ymd HM", "ymd H", "ymd", 
                                                                    "Ymd HMS", "Ymd HM", "Ymd H", "Ymd"),
                                                          quiet = TRUE),
      is.numeric(FINAL_RESPONSE_DT) ~ as.POSIXct(FINAL_RESPONSE_DT, origin = "1970-01-01"),
      TRUE ~ as.POSIXct(NA)
    )
  ) %>%
  filter(!is.na(FINAL_RESPONSE_DT)) %>%  # Remove rows where conversion failed
  left_join(listing_center_codes, by = "WL_ID_CODE") %>%
  filter(!is.na(LISTING_CTR_CODE)) %>%
  left_join(cohort_keys, by = c("DONOR_ID", "MATCH_ID_CODE")) %>%
  filter(!is.na(allocation_type)) %>%
  # EXCLUDE refusals that meet open offer definition
  # For refusals (OFFER_ACCEPT == "N"), open offer = PRIME_OPO_REFUSAL_ID in text codes AND text pattern matches
  mutate(
    is_open_offer_refusal = (PRIME_OPO_REFUSAL_ID %in% OPEN_OFFER_TEXT_CODES &
                             str_detect(tidyr::replace_na(OPO_OSTXT_PRIMARY, ""), OPEN_OFFER_TEXT_PATTERN))
  ) %>%
  filter(!is_open_offer_refusal)  # Exclude open offer refusals from multi-refusal calculation

cat("Total refusals for cohort donor-runs (excluding open offer refusals):", nrow(refusals), "\n")

# Identify multi-refusals: same DONOR_ID, MATCH_ID_CODE, LISTING_CTR_CODE, FINAL_RESPONSE_DT
multi_refusal_groups <- refusals %>%
  group_by(DONOR_ID, MATCH_ID_CODE, LISTING_CTR_CODE, FINAL_RESPONSE_DT, allocation_type) %>%
  summarise(n_recipients = n_distinct(WL_ID_CODE), .groups = "drop") %>%
  mutate(is_multi_refusal = n_recipients > 1)

# Count recipients bundled in multi-refusals per donor-run
recipients_in_multi_refusals <- refusals %>%
  inner_join(
    multi_refusal_groups %>%
      filter(is_multi_refusal) %>%
      select(DONOR_ID, MATCH_ID_CODE, LISTING_CTR_CODE, FINAL_RESPONSE_DT),
    by = c("DONOR_ID", "MATCH_ID_CODE", "LISTING_CTR_CODE", "FINAL_RESPONSE_DT")
  ) %>%
  group_by(DONOR_ID, MATCH_ID_CODE, allocation_type) %>%
  summarise(n_recipients_bundled = n_distinct(WL_ID_CODE), .groups = "drop")

# Count distinct refusals by centers that refused multiple times (back-and-forth) per donor-run
distinct_refusals_same_center <- multi_refusal_groups %>%
  filter(!is_multi_refusal) %>%
  group_by(DONOR_ID, MATCH_ID_CODE, LISTING_CTR_CODE, allocation_type) %>%
  summarise(n_refusals_by_center = n(), .groups = "drop") %>%
  filter(n_refusals_by_center > 1) %>%
  group_by(DONOR_ID, MATCH_ID_CODE, allocation_type) %>%
  summarise(n_distinct_refusals_same_center = sum(n_refusals_by_center), .groups = "drop")

# Count per donor-run
donor_run_refusal_stats <- multi_refusal_groups %>%
  group_by(DONOR_ID, MATCH_ID_CODE, allocation_type) %>%
  summarise(
    n_multi_refusal_decisions = sum(is_multi_refusal),
    n_distinct_refusals = sum(!is_multi_refusal),
    .groups = "drop"
  ) %>%
  left_join(recipients_in_multi_refusals, by = c("DONOR_ID", "MATCH_ID_CODE", "allocation_type")) %>%
  left_join(distinct_refusals_same_center, by = c("DONOR_ID", "MATCH_ID_CODE", "allocation_type")) %>%
  left_join(
    cohort_compare %>% 
      select(DONOR_ID, MATCH_ID_CODE, n_unique_ctr_contacted) %>%
      distinct(),
    by = c("DONOR_ID", "MATCH_ID_CODE")
  ) %>%
  mutate(
    n_recipients_bundled = coalesce(n_recipients_bundled, 0L),
    n_distinct_refusals_same_center = coalesce(n_distinct_refusals_same_center, 0L)
  )

cat("Donor-runs with refusal statistics:", nrow(donor_run_refusal_stats), "\n")


# CALCULATE MAX LAB VALUES FROM LABPANELS
# Get unique donors from final_cohort
final_cohort_donors <- final_cohort %>%
  distinct(DONOR_ID)

# Filter labpanels to only donors in final_cohort and calculate max values
labpanels_max <- labpanels %>%
  inner_join(final_cohort_donors, by = "DONOR_ID") %>%
  group_by(DONOR_ID) %>%
  summarise(
    max_SGOT = max(SGOT, na.rm = TRUE),
    max_SGPT = max(SGPT, na.rm = TRUE),
    max_BILIRUBIN = max(BILIRUBIN, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Replace -Inf with NA (occurs when all values are NA)
  mutate(
    max_SGOT = if_else(is.infinite(max_SGOT), NA_real_, max_SGOT),
    max_SGPT = if_else(is.infinite(max_SGPT), NA_real_, max_SGPT),
    max_BILIRUBIN = if_else(is.infinite(max_BILIRUBIN), NA_real_, max_BILIRUBIN)
  )

# Join max values to final_cohort (at DONOR_ID level, will join to all rows for that donor)
final_cohort <- final_cohort %>%
  left_join(labpanels_max, by = "DONOR_ID")

cat("- Donors with max_SGOT:", sum(!is.na(final_cohort$max_SGOT)), "\n")
cat("- Donors with max_SGPT:", sum(!is.na(final_cohort$max_SGPT)), "\n")
cat("- Donors with max_BILIRUBIN:", sum(!is.na(final_cohort$max_BILIRUBIN)), "\n\n")

# CREATE STANDARDIZED DISCARDED DATASET WITH COMPLETE DONOR TYPE DATA

# AOOS identification constants (matching 02_build_cohort.r)
OPEN_OFFER_BYPASS_CODES <- c(861L, 862L, 863L, 887L, 888L, 760L, 761L, 762L, 763L, 764L, 765L)
OPEN_OFFER_TEXT_CODES <- c(798L, 799L, 898L)
OPEN_OFFER_TEXT_PATTERN <- regex("expedited|open offer|aggressive", ignore_case = TRUE)

# Filter PTR to study period (matching date_lb from 01_read_data.r)
ptr_filtered <- ptr_reduced %>%
  filter(!is.na(MATCH_DATE), as.Date(MATCH_DATE) >= as.Date("2021-12-09"))

# Get all unique donors in filtered PTR
all_donors_ptr <- ptr_filtered %>%
  distinct(DONOR_ID)

# Get transplanted donors (in liver_data)
transplanted_donors <- liver_reduced %>%
  filter(!is.na(TX_DATE), !is.na(DONOR_ID)) %>%
  distinct(DONOR_ID)

# Discarded = donors NOT in liver_data
discarded_donors <- all_donors_ptr %>%
  filter(!DONOR_ID %in% transplanted_donors$DONOR_ID)

cat("Total PTR donors (study period):", nrow(all_donors_ptr), "\n")
cat("Transplanted donors:", nrow(transplanted_donors), "\n")
cat("Discarded donors (total):", nrow(discarded_donors), "\n")

# Join with donor_reduced to get donor type
discarded_with_donor_type <- discarded_donors %>%
  inner_join(donor_reduced %>% select(DONOR_ID, NON_HRT_DON), by = "DONOR_ID") %>%
  mutate(
    has_donor_type = !is.na(NON_HRT_DON),
    donor_type = case_when(
      NON_HRT_DON %in% c("Y", "YES", "1", "TRUE") ~ "DCD",
      NON_HRT_DON %in% c("N", "NO", "0", "FALSE") ~ "DBD",
      TRUE ~ NA_character_
    )
  )

cat("Discarded with donor type data:", sum(discarded_with_donor_type$has_donor_type, na.rm = TRUE), 
    "(", round(100 * sum(discarded_with_donor_type$has_donor_type, na.rm = TRUE) / nrow(discarded_donors), 1), "%)\n")

# Identify all AOOS match-runs (for donor-level classification)
aoos_keys_all <- ptr_filtered %>%
  filter(
    ((OFFER_ACCEPT == "B") &
       (PRIME_OPO_REFUSAL_ID %in% OPEN_OFFER_BYPASS_CODES |
          (PRIME_OPO_REFUSAL_ID %in% OPEN_OFFER_TEXT_CODES &
             str_detect(tidyr::replace_na(OPO_OSTXT_PRIMARY, ""), OPEN_OFFER_TEXT_PATTERN)))) |
      ((OFFER_ACCEPT == "N") &
         PRIME_OPO_REFUSAL_ID %in% OPEN_OFFER_TEXT_CODES &
         str_detect(tidyr::replace_na(OPO_OSTXT_PRIMARY, ""), OPEN_OFFER_TEXT_PATTERN))
  ) %>%
  distinct(DONOR_ID, MATCH_ID_CODE)

# For each discarded donor, check if ANY match-run was AOOS (donor-level)
donor_aoos_status <- ptr_filtered %>%
  distinct(DONOR_ID, MATCH_ID_CODE) %>%
  left_join(aoos_keys_all %>% mutate(is_aoos = TRUE), by = c("DONOR_ID", "MATCH_ID_CODE")) %>%
  group_by(DONOR_ID) %>%
  summarise(
    has_aoos_match_run = any(!is.na(is_aoos) & is_aoos == TRUE, na.rm = TRUE),
    .groups = "drop"
  )

# Final discarded with complete data (has donor type) 
discarded_complete <- discarded_with_donor_type %>%
  filter(has_donor_type == TRUE) %>%
  left_join(donor_aoos_status, by = "DONOR_ID") %>%
  mutate(
    allocation_type = if_else(coalesce(has_aoos_match_run, FALSE), "AOOS", "Sequential"),
    status = "Discarded"
  ) %>%
  select(DONOR_ID, allocation_type, status, donor_type)

cat("Discarded with complete data (final):", nrow(discarded_complete), "\n")
cat("  DBD:", sum(discarded_complete$donor_type == "DBD", na.rm = TRUE), "\n")
cat("  DCD:", sum(discarded_complete$donor_type == "DCD", na.rm = TRUE), "\n")
cat("  AOOS:", sum(discarded_complete$allocation_type == "AOOS", na.rm = TRUE), "\n")
cat("  In-Sequence:", sum(discarded_complete$allocation_type == "Sequential", na.rm = TRUE), "\n\n")

# ==============================================================================
# SAVE FINAL COHORT WITH ALL VARIABLES, REFUSAL STATISTICS, AND DISCARDED COMPLETE
# ==============================================================================
save(final_cohort, donor_run_refusal_stats, discarded_complete,
     file = here("data_processed", "03_cohort_with_variables.rda"))

cat("Final cohort size:", nrow(final_cohort), "\n")
cat("Total variables:", ncol(final_cohort), "\n")
cat("Donor-run refusal statistics:", nrow(donor_run_refusal_stats), "\n")
cat("Saved to: data_processed/03_cohort_with_variables.rda\n")