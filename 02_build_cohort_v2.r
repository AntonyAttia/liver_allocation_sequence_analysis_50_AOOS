# 02_build_cohort_v2.r
# Applying exclusion criteria and creation of final cohort
# v2: 4-group classification — <50, ≥50, AOOS (not EP), EP (887)
#     EP = bypass code 887 (policy-compliant expedited placement, Policy 9.10)
#     AOOS = 861, 862, 863, 888 (bypass) + 798/799/898 with text match — excluding 887
#     760-765 excluded (not in OPTN analytical definition, only included in Wehrle paper)

library(tidyverse)

# Same path as in 01_read_data_v2.r (where the .rda files are saved/loaded)
DATA_DIR <- "/Users/antonyattia/Desktop/Stanford/Open Offers/Manuscript/liver_allocation_sequence_analysis_50_AOOS/data_processed"

load(file.path(DATA_DIR, "01_raw_data_v2.rda"))

# ============================================================
# CONFIGURATION
# ============================================================

DIRECTED_DONATION_REFUSAL_CODE <- 851L

# EP: policy-compliant expedited placement (Policy 9.10)
EP_BYPASS_CODE <- 887L

# AOOS bypass codes (NOT including 887)
AOOS_BYPASS_CODES <- c(861L, 862L, 863L, 888L)

# Text-based codes (need pattern match)
AOOS_TEXT_CODES   <- c(798L, 799L, 898L)
AOOS_TEXT_PATTERN <- regex("expedited|open offer|aggressive", ignore_case = TRUE)

# Combined: all codes that flag a match run as "any out-of-sequence" (EP + AOOS)
ALL_OOS_BYPASS_CODES <- c(EP_BYPASS_CODE, AOOS_BYPASS_CODES)

# ============================================================
# STEP 1: ALL TRANSPLANTS IN STUDY PERIOD
# ============================================================
all_transplants <- liver_reduced %>%
  filter(!is.na(WL_ID_CODE), !is.na(TX_DATE), !is.na(DONOR_ID)) %>%
  filter(as.Date(TX_DATE) >= date_lb) %>%
  select(WL_ID_CODE, DONOR_ID, TX_DATE, DON_TY, MULTIORG)

cat("Total transplants in study period:", nrow(all_transplants), "\n\n")

# ============================================================
# STEP 2: JOIN WITH ACCEPTED OFFERS
# ============================================================
accepted_offers_with_tx <- ptr_reduced %>%
  filter(OFFER_ACCEPT == "Y") %>%
  inner_join(all_transplants, by = c("WL_ID_CODE", "DONOR_ID")) %>%
  transmute(
    WL_ID_CODE, DONOR_ID, MATCH_ID_CODE,
    PTR_SEQUENCE_NUM_tx = PTR_SEQUENCE_NUM,
    MATCH_DATE_tx = MATCH_DATE,
    INITIAL_RESPONSE_DT_tx = INITIAL_RESPONSE_DT,
    FINAL_RESPONSE_DT_tx = FINAL_RESPONSE_DT,
    TX_DATE, DON_TY, MULTIORG
  ) %>%
  distinct()

cat("Transplants with accepted offer records:", nrow(accepted_offers_with_tx), "\n\n")

# ============================================================
# STEP 3: EXCLUDE MULTI-ORGAN TRANSPLANTS
# ============================================================
after_multiorg_excl <- accepted_offers_with_tx %>%
  filter(MULTIORG != "Y")

cat("After multi-organ exclusions:", nrow(after_multiorg_excl), "\n\n")

# ============================================================
# STEP 4: EXCLUDE LIVING DONOR TRANSPLANTS
# ============================================================
after_living_excl <- after_multiorg_excl %>%
  filter(DON_TY != "L")

cat("After living donor exclusions:", nrow(after_living_excl), "\n\n")

# ============================================================
# STEP 5: EXCLUDE SPLIT LIVER CASES
# ============================================================
split_donors <- all_transplants %>%
  count(DONOR_ID) %>%
  filter(n > 1) %>%
  pull(DONOR_ID)

after_split_excl <- after_living_excl %>%
  filter(!DONOR_ID %in% split_donors)

cat("After split liver exclusions:", nrow(after_split_excl), "\n\n")

# ============================================================
# STEP 6: RESOLVE PTR-LEVEL DUPLICATES
# ============================================================
# After split exclusion, liver_reduced has exactly 1 transplant per donor.
# But PTR can have >1 acceptance row for the same (WL_ID_CODE, DONOR_ID)
# when the same patient accepted on 2 different match runs (re-runs).
# Diagnostic confirmed all 39 such cases are "different match runs, same recipient."
# Resolution: keep the row with the lowest PTR_SEQUENCE_NUM (primary acceptance).

n_before_dedup <- nrow(after_split_excl)

after_dedup <- after_split_excl %>%
  group_by(DONOR_ID) %>%
  arrange(PTR_SEQUENCE_NUM_tx) %>%
  slice(1) %>%
  ungroup()

n_ptr_dups_resolved <- n_before_dedup - nrow(after_dedup)
cat("After resolving PTR duplicates:", nrow(after_dedup),
    "(resolved", n_ptr_dups_resolved, "duplicate rows from",
    n_ptr_dups_resolved, "donors by keeping lowest seq#)\n\n")

# ============================================================
# STEP 7: EXCLUDE DIRECTED DONATIONS
# ============================================================
used_seq_by_match <- after_dedup %>%
  transmute(
    DONOR_ID, MATCH_ID_CODE,
    WL_ID_CODE_used = WL_ID_CODE,
    used_seq = PTR_SEQUENCE_NUM_tx
  ) %>%
  distinct()

ptr_directed <- ptr_reduced %>%
  filter(PRIME_OPO_REFUSAL_ID == DIRECTED_DONATION_REFUSAL_CODE) %>%
  select(DONOR_ID, MATCH_ID_CODE, WL_ID_CODE, PTR_SEQUENCE_NUM)

directed_exclude_keys <- ptr_directed %>%
  inner_join(used_seq_by_match, by = c("DONOR_ID", "MATCH_ID_CODE")) %>%
  filter(PTR_SEQUENCE_NUM < used_seq, WL_ID_CODE != WL_ID_CODE_used) %>%
  distinct(DONOR_ID, MATCH_ID_CODE)

after_directed_excl <- after_dedup %>%
  anti_join(directed_exclude_keys, by = c("DONOR_ID", "MATCH_ID_CODE"))

cat("After directed donation exclusions:", nrow(after_directed_excl), "\n\n")

# ============================================================
# STEP 8: CLASSIFY — EP vs AOOS vs IN-SEQUENCE
# ============================================================

# 8a. Identify EP match-runs (has 887 bypass)
ep_keys <- ptr_reduced %>%
  filter(OFFER_ACCEPT == "B", PRIME_OPO_REFUSAL_ID == EP_BYPASS_CODE) %>%
  distinct(DONOR_ID, MATCH_ID_CODE) %>%
  mutate(has_ep = TRUE)

# 8b. Identify AOOS match-runs (has AOOS bypass codes OR text-matched codes, but NOT EP)
aoos_bypass_keys <- ptr_reduced %>%
  filter(OFFER_ACCEPT == "B", PRIME_OPO_REFUSAL_ID %in% AOOS_BYPASS_CODES) %>%
  distinct(DONOR_ID, MATCH_ID_CODE)

aoos_text_keys <- ptr_reduced %>%
  filter(
    OFFER_ACCEPT %in% c("B", "N"),
    PRIME_OPO_REFUSAL_ID %in% AOOS_TEXT_CODES,
    str_detect(replace_na(OPO_OSTXT_PRIMARY, ""), AOOS_TEXT_PATTERN)
  ) %>%
  distinct(DONOR_ID, MATCH_ID_CODE)

aoos_keys <- bind_rows(aoos_bypass_keys, aoos_text_keys) %>%
  distinct() %>%
  mutate(has_aoos = TRUE)

# 8c. MATCH-RUN-LEVEL classification (not donor-level)
#     Each transplant has a specific (DONOR_ID, MATCH_ID_CODE) — the match run
#     where the organ was accepted. We check if THAT match run had EP or AOOS codes.
#     This avoids misclassifying cases where EP was triggered on a separate match run
#     but the organ was actually placed in-sequence on the original run.
cohort_with_seq <- after_directed_excl %>%
  left_join(ep_keys,   by = c("DONOR_ID", "MATCH_ID_CODE")) %>%
  left_join(aoos_keys, by = c("DONOR_ID", "MATCH_ID_CODE")) %>%
  mutate(
    has_ep   = coalesce(has_ep, FALSE),
    has_aoos = coalesce(has_aoos, FALSE),
    is_ep      = has_ep,
    is_aoos    = has_aoos & !has_ep,
    is_any_oos = has_ep | has_aoos,
    seq_flag = case_when(
      is_ep                      ~ "EP (Policy 9.10)",
      is_aoos                    ~ "AOOS (out-of-sequence)",
      PTR_SEQUENCE_NUM_tx < 50   ~ "<50 sequence",
      PTR_SEQUENCE_NUM_tx >= 50  ~ "\u226550 sequence"
    ),
    seq_flag = factor(seq_flag, levels = c(
      "<50 sequence", "\u226550 sequence", "AOOS (out-of-sequence)", "EP (Policy 9.10)"
    ))
  ) %>%
  select(-has_ep, -has_aoos)

cat("After classification:", nrow(cohort_with_seq), "\n")
print(table(cohort_with_seq$seq_flag))
cat("\n")

# ============================================================
# JOIN WITH FULL RECIPIENT AND DONOR DATA
# ============================================================
cohort_with_data <- cohort_with_seq %>%
  left_join(ptr_reduced,
            by = c("DONOR_ID", "MATCH_ID_CODE", "WL_ID_CODE",
                   "PTR_SEQUENCE_NUM_tx" = "PTR_SEQUENCE_NUM"),
            suffix = c("", "_offer")) %>%
  left_join(liver_reduced %>% distinct(DONOR_ID, WL_ID_CODE, .keep_all = TRUE),
            by = c("DONOR_ID", "WL_ID_CODE"),
            suffix = c("", "_recip")) %>%
  left_join(donor_reduced %>% distinct(DONOR_ID, .keep_all = TRUE),
            by = "DONOR_ID",
            suffix = c("", "_donor"))

drop_cols <- c(
  "AGE_DON", "ALCOHOL_HEAVY_DON", "BMI_DON_CALC", "COD_CAD_DON",
  "GENDER_DON", "MACRO_FAT_LI_DON", "NON_HRT_DON", "TBILI_DON",
  "MATCH_DATE", "INITIAL_RESPONSE_DT", "FINAL_RESPONSE_DT",
  "TX_DATE_recip"
)

cohort_with_data <- cohort_with_data %>%
  select(-any_of(drop_cols))

cat("After joining recipient and donor data:", nrow(cohort_with_data), "\n\n")

# ============================================================
# STEP 9: EXCLUDE MISSING CRITICAL VARIABLES
# ============================================================
n_before_step9 <- nrow(cohort_with_data)

n_excluded_meld <- sum(is.na(cohort_with_data$MELD_PELD_LAB_SCORE))
n_excluded_donor <- sum(
  is.na(cohort_with_data$AGE_DON_donor) &
    is.na(cohort_with_data$COD_CAD_DON_donor) &
    is.na(cohort_with_data$NON_HRT_DON_donor) &
    (is.na(cohort_with_data$GENDER_DON_donor) |
       trimws(cohort_with_data$GENDER_DON_donor) == "")
)

final_cohort <- cohort_with_data %>%
  filter(
    !(is.na(AGE_DON_donor) & is.na(COD_CAD_DON_donor) &
        is.na(NON_HRT_DON_donor) &
        (is.na(GENDER_DON_donor) | trimws(GENDER_DON_donor) == ""))
  ) %>%
  filter(!is.na(MELD_PELD_LAB_SCORE))

cat("After exclusions:", nrow(final_cohort), "\n")
cat("  Excluded for missing MELD:", n_excluded_meld, "\n")
cat("  Excluded for missing donor data:", n_excluded_donor, "\n\n")

# ============================================================
# SAVE COUNTS FOR CONSORT DIAGRAM
# ============================================================
counts <- list(
  n_initial = nrow(all_transplants),
  n_after_join = nrow(accepted_offers_with_tx),
  n_after_multiorg = nrow(after_multiorg_excl),
  n_after_living = nrow(after_living_excl),
  n_after_split = nrow(after_split_excl),
  n_after_dedup = nrow(after_dedup),
  n_ptr_dups_resolved = n_ptr_dups_resolved,
  n_after_directed = nrow(after_directed_excl),
  n_before_step9 = n_before_step9,
  n_final = nrow(final_cohort),
  n_excluded_meld = n_excluded_meld,
  n_excluded_donor = n_excluded_donor,
  final_cohort_breakdown = table(final_cohort$seq_flag)
)

# ============================================================
# SAVE
# ============================================================
save(final_cohort, file = file.path(DATA_DIR, "02_cohort_built_v2.rda"))
saveRDS(counts, file.path(DATA_DIR, "02_cohort_counts_v2.rds"))

cat("Final cohort:", nrow(final_cohort), "cases\n")
cat("Breakdown:\n")
print(table(final_cohort$seq_flag))
cat("\nSaved to: data_processed/02_cohort_built_v2.rda\n")
