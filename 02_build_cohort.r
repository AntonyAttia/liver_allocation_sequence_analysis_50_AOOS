# ==============================================================================
# 02_build_cohort.r
# Applying exclusion criteria and creation of final cohort
# ==============================================================================

library(tidyverse)
library(here)
library(DiagrammeR)
library(glue)

# Load data from previous step
load(here("data_processed", "01_raw_data.rda"))

# Configuration
DIRECTED_DONATION_REFUSAL_CODE <- 851L
OPEN_OFFER_BYPASS_CODES <- c(861L, 862L, 863L, 887L, 888L, 760L, 761L, 762L, 763L, 764L, 765L)
OPEN_OFFER_TEXT_CODES <- c(798L, 799L, 898L)
OPEN_OFFER_TEXT_PATTERN <- regex("expedited|open offer|aggressive", ignore_case = TRUE)

# ==============================================================================
# STEP 1: ALL TRANSPLANTS IN STUDY PERIOD
# ==============================================================================
all_transplants <- liver_reduced %>%
  filter(!is.na(WL_ID_CODE), !is.na(TX_DATE), !is.na(DONOR_ID)) %>%
  filter(as.Date(TX_DATE) >= date_lb) %>%
  select(WL_ID_CODE, DONOR_ID, TX_DATE, DON_TY, MULTIORG)

cat("=== STEP 1: ALL TRANSPLANTS ===\n")
cat("Total transplants in study period:", nrow(all_transplants), "\n\n")

# ==============================================================================
# STEP 2: JOIN WITH ACCEPTED OFFERS
# ==============================================================================
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

cat("=== STEP 2: JOIN WITH ACCEPTED OFFERS ===\n")
cat("Transplants with accepted offer records:", nrow(accepted_offers_with_tx), "\n\n")

# ==============================================================================
# STEP 3: EXCLUDE MULTI-ORGAN TRANSPLANTS
# ==============================================================================
after_multiorg_excl <- accepted_offers_with_tx %>%
  filter(MULTIORG != "Y")

cat("=== STEP 3: EXCLUDE MULTI-ORGAN ===\n")
cat("After multi-organ exclusions:", nrow(after_multiorg_excl), "\n\n")

# ==============================================================================
# STEP 4: EXCLUDE LIVING DONOR TRANSPLANTS
# ==============================================================================
after_living_excl <- after_multiorg_excl %>%
  filter(DON_TY != "L")

cat("=== STEP 4: EXCLUDE LIVING DONOR ===\n")
cat("After living donor exclusions:", nrow(after_living_excl), "\n\n")

# ==============================================================================
# STEP 5: EXCLUDE SPLIT LIVER CASES
# ==============================================================================
split_donors <- all_transplants %>%
  count(DONOR_ID) %>%
  filter(n > 1) %>%
  pull(DONOR_ID)

after_split_excl <- after_living_excl %>%
  filter(!DONOR_ID %in% split_donors)

cat("=== STEP 5: EXCLUDE SPLIT LIVER ===\n")
cat("After split liver exclusions:", nrow(after_split_excl), "\n\n")

# ==============================================================================
# STEP 6: REMOVE DUPLICATE DONORS
# ==============================================================================
donors_with_duplicates <- after_split_excl %>%
  group_by(DONOR_ID) %>%
  filter(n() > 1) %>%
  pull(DONOR_ID) %>%
  unique()

after_dedup <- after_split_excl %>%
  filter(!DONOR_ID %in% donors_with_duplicates)

cat("=== STEP 6: REMOVE DUPLICATE DONORS ===\n")
cat("After donor deduplication:", nrow(after_dedup), "\n\n")

# ==============================================================================
# STEP 7: EXCLUDE DIRECTED DONATIONS
# ==============================================================================
used_seq_by_match <- after_dedup %>%
  transmute(
    DONOR_ID, MATCH_ID_CODE,
    WL_ID_CODE_used = WL_ID_CODE,
    used_seq = PTR_SEQUENCE_NUM_tx
  ) %>%
  distinct()

ptr_directed <- ptr_reduced %>%
  filter(PRIME_OPO_REFUSAL_ID == 851) %>%
  select(DONOR_ID, MATCH_ID_CODE, WL_ID_CODE, PTR_SEQUENCE_NUM)

directed_exclude_keys <- ptr_directed %>%
  inner_join(used_seq_by_match, by = c("DONOR_ID", "MATCH_ID_CODE")) %>%
  filter(PTR_SEQUENCE_NUM < used_seq, WL_ID_CODE != WL_ID_CODE_used) %>%
  distinct(DONOR_ID, MATCH_ID_CODE)

after_directed_excl <- after_dedup %>%
  anti_join(directed_exclude_keys, by = c("DONOR_ID","MATCH_ID_CODE"))

cat("=== STEP 7: EXCLUDE DIRECTED DONATIONS ===\n")
cat("After directed donation exclusions:", nrow(after_directed_excl), "\n\n")

# ==============================================================================
# STEP 8: CLASSIFY OPEN OFFERS
# ==============================================================================
oo_keys <- ptr_reduced %>%
  filter(
    ((OFFER_ACCEPT == "B") &
       (PRIME_OPO_REFUSAL_ID %in% OPEN_OFFER_BYPASS_CODES |
          (PRIME_OPO_REFUSAL_ID %in% OPEN_OFFER_TEXT_CODES &
             str_detect(tidyr::replace_na(OPO_OSTXT_PRIMARY, ""), OPEN_OFFER_TEXT_PATTERN))))
    |
      ((OFFER_ACCEPT == "N") &
         PRIME_OPO_REFUSAL_ID %in% OPEN_OFFER_TEXT_CODES &
         str_detect(tidyr::replace_na(OPO_OSTXT_PRIMARY, ""), OPEN_OFFER_TEXT_PATTERN))
  ) %>%
  distinct(DONOR_ID, MATCH_ID_CODE)

cohort_with_seq <- after_directed_excl %>%
  left_join(oo_keys %>% mutate(is_open_offer = TRUE),
            by = c("DONOR_ID","MATCH_ID_CODE")) %>%
  mutate(
    is_open_offer = coalesce(is_open_offer, FALSE),
    seq_flag = case_when(
      is_open_offer ~ "AOOS (out-of-sequence)",
      PTR_SEQUENCE_NUM_tx < 50 ~ "<50 sequence",
      PTR_SEQUENCE_NUM_tx >= 50 ~ "≥50 sequence"
    )
  )

cat("=== STEP 8: OPEN OFFER CLASSIFICATION ===\n")
cat("After classification:", nrow(cohort_with_seq), "\n")
print(table(cohort_with_seq$seq_flag))
cat("\n")

# ==============================================================================
# JOIN WITH FULL RECIPIENT AND DONOR DATA
# ==============================================================================
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

# Drop duplicate columns
drop_cols <- c(
  "AGE_DON","ALCOHOL_HEAVY_DON","BMI_DON_CALC","COD_CAD_DON",
  "GENDER_DON","MACRO_FAT_LI_DON","NON_HRT_DON","TBILI_DON",
  "MATCH_DATE","INITIAL_RESPONSE_DT","FINAL_RESPONSE_DT",
  "TX_DATE_recip"
)

cohort_with_data <- cohort_with_data %>%
  select(-any_of(drop_cols))

cat("After joining recipient and donor data:", nrow(cohort_with_data), "\n\n")

# ==============================================================================
# STEP 9: EXCLUDE MISSING CRITICAL VARIABLES
# ==============================================================================

n_before_step9 <- nrow(cohort_with_data)

# Count exclusions before filtering
n_excluded_meld <- sum(is.na(cohort_with_data$MELD_PELD_LAB_SCORE))
n_excluded_donor <- sum(
  is.na(cohort_with_data$AGE_DON_donor) & 
    is.na(cohort_with_data$COD_CAD_DON_donor) & 
    is.na(cohort_with_data$NON_HRT_DON_donor) & 
    (is.na(cohort_with_data$GENDER_DON_donor) | 
       trimws(cohort_with_data$GENDER_DON_donor) == "")
)

# Apply exclusions
final_cohort <- cohort_with_data %>%
  filter(
    !(is.na(AGE_DON_donor) & is.na(COD_CAD_DON_donor) & 
        is.na(NON_HRT_DON_donor) & 
        (is.na(GENDER_DON_donor) | trimws(GENDER_DON_donor) == ""))
  ) %>%
  filter(!is.na(MELD_PELD_LAB_SCORE))

cat("=== STEP 9: EXCLUDE MISSING CRITICAL VARIABLES ===\n")
cat("After exclusions:", nrow(final_cohort), "\n")
cat("  - Excluded for missing MELD:", n_excluded_meld, "\n")
cat("  - Excluded for missing donor data:", n_excluded_donor, "\n\n")

# ==============================================================================
# CREATE COHORT FLOW DIAGRAM
# ==============================================================================
n_initial <- nrow(all_transplants)
n_after_join <- nrow(accepted_offers_with_tx)
n_after_multiorg <- nrow(after_multiorg_excl)
n_after_living <- nrow(after_living_excl)
n_after_split <- nrow(after_split_excl)
n_after_dedup <- nrow(after_dedup)
n_after_directed <- nrow(after_directed_excl)
n_final <- nrow(final_cohort)

final_cohort_breakdown <- table(final_cohort$seq_flag)
breakdown_text <- paste0(
  "<50 sequence: ", final_cohort_breakdown["<50 sequence"], "\\n",
  "≥50 sequence: ", final_cohort_breakdown["≥50 sequence"], "\\n",
  "AOOS: ", final_cohort_breakdown["AOOS (out-of-sequence)"]
)

#consort_diagram <- grViz(glue("
#digraph consort {{
#  graph [rankdir = TB, splines = ortho, nodesep=0.35, ranksep=0.45]
#  node [shape=box, style='rounded,filled', fillcolor='white', color='#666666', fontname='Helvetica', fontsize=10]
#  edge [color='#666666', arrowsize=0.6]
#  
 # A [label='All transplanted patients (STAR/PTR)\\nJan 2018 - April 2025\\nN = {format(n_initial, big.mark=',')}', fillcolor='#F8F8F8']
 # B [label='With accepted offer records\\nN = {format(n_after_join, big.mark=',')}']
  #C [label='After multi-organ exclusions\\nN = {format(n_after_multiorg, big.mark=',')}']
  #D [label='After living donor exclusions\\nN = {format(n_after_living, big.mark=',')}']
  #E [label='After split liver exclusions\\nN = {format(n_after_split, big.mark=',')}']
  #F [label='After donor deduplication\\nN = {format(n_after_dedup, big.mark=',')}']
#  G [label='After directed donation exclusions\\nN = {format(n_after_directed, big.mark=',')}']
 # H [label='Final analysis cohort\\nN = {format(n_final, big.mark=',')}\\n\\n{breakdown_text}', fillcolor='#EFEFEF']
#  
#  X1 [label='Excluded:\\nOffers not in PTR dataset\\n(large portion are living donors)\\nN = {format(n_initial - n_after_join, big.mark=',')}', 
#      shape=note, fillcolor='#FFE6E6', color='#CC0000']
#  X2 [label='Excluded:\\nMulti-organ transplants\\nN = {format(n_after_join - n_after_multiorg, big.mark=',')}', 
#      shape=note, fillcolor='white', color='#AAAAAA']
#  X3 [label='Excluded:\\nLiving donor transplants\\nN = {format(n_after_multiorg - n_after_living, big.mark=',')}', 
#      shape=note, fillcolor='white', color='#AAAAAA']
#  X4 [label='Excluded:\\nSplit liver cases\\nN = {format(n_after_living - n_after_split, big.mark=',')}', 
#      shape=note, fillcolor='white', color='#AAAAAA']
#  X5 [label='Excluded:\\nDuplicate donors\\nN = {format(n_after_split - n_after_dedup, big.mark=',')}', 
#      shape=note, fillcolor='white', color='#AAAAAA']
#  X6 [label='Excluded:\\nDirected donations\\nN = {format(n_after_dedup - n_after_directed, big.mark=',')}', 
#      shape=note, fillcolor='white', color='#AAAAAA']
#  X7 [label='Excluded:\\nMissing MELD: {n_excluded_meld}\\nMissing donor data: {n_excluded_donor}\\nTotal N = {format(n_before_step9 - n_final, big.mark=',')}', 
#      shape=note, fillcolor='white', color='#AAAAAA']
  
#  A -> B
#  B -> C
#  C -> D
#  D -> E
#  E -> F
#  F -> G
#  G -> H
#  
#  A -> X1 [arrowhead=none, color='#CC0000', style=bold]
#  B -> X2 [arrowhead=none, color='#BBBBBB']
#  C -> X3 [arrowhead=none, color='#BBBBBB']
#  D -> X4 [arrowhead=none, color='#BBBBBB']
#  E -> X5 [arrowhead=none, color='#BBBBBB']
#  F -> X6 [arrowhead=none, color='#BBBBBB']
#  G -> X7 [arrowhead=none, color='#BBBBBB']
#}}
#"))

#print(consort_diagram)
#cat("\nCONSORT diagram created\n\n")

counts <- list(
  n_initial = n_initial,
  n_after_join = n_after_join,
  n_after_multiorg = n_after_multiorg,
  n_after_living = n_after_living,
  n_after_split = n_after_split,
  n_after_dedup = n_after_dedup,
  n_after_directed = n_after_directed,
  n_before_step9 = n_before_step9,
  n_final = n_final,
  n_excluded_meld = n_excluded_meld,
  n_excluded_donor = n_excluded_donor,
  final_cohort_breakdown = table(final_cohort$seq_flag)
)

save(final_cohort, file = here("data_processed", "02_cohort_built.rda"))
saveRDS(counts, here("data_processed", "02_cohort_counts.rds"))

# ==============================================================================
# SAVE FINAL COHORT
# ==============================================================================
save(final_cohort, 
     file = here("data_processed", "02_cohort_built.rda"))

cat("Cohort building complete!\n")
cat("Final cohort:", nrow(final_cohort), "cases\n")
cat("Saved to: data_processed/02_cohort_built.rda\n")