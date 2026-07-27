# table1.r: Table 1 (4 groups) + Supplemental Tables S2 (DBD) and S3 (DCD)
dir.create(here("figures", "main"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("figures", "supplement"), showWarnings = FALSE, recursive = TRUE)
# Table 1: 4-group comparison
tbl1 <- timing_cohort %>%
  select(
    seq_flag,
    AGE, MELD_PELD_ADJ, GENDER,
    AGE_DON_donor, donation_type, BMI_DON_CALC_donor,
    MACRO_FAT_CAT, max_BILIRUBIN, max_SGOT, max_SGPT,
    DRI, DSRI,
    COLD_ISCH_cat, DISTANCE,
    LI_MACHINE_PERFUSION_original_display, LI_MACHINE_PERFUSION_FDA_era
  ) %>%
  tbl_summary(
    by = seq_flag,
    label = list(
      AGE ~ "Recipient age (years)",
      MELD_PELD_ADJ ~ "Recipient MELD/PELD",
      GENDER ~ "Recipient gender",
      AGE_DON_donor ~ "Donor age (years)",
      donation_type ~ "Donation type",
      BMI_DON_CALC_donor ~ "Donor BMI (kg/m\u00B2)",
      MACRO_FAT_CAT ~ "Donor macrosteatosis",
      max_BILIRUBIN ~ "Donor maximum bilirubin (mg/dL)",
      max_SGOT ~ "Donor maximum AST (U/L)",
      max_SGPT ~ "Donor maximum ALT (U/L)",
      DRI ~ "Donor Risk Index (DRI)",
      DSRI ~ "Discard Risk Index (DSRI)",
      COLD_ISCH_cat ~ "Cold ischemia time",
      DISTANCE ~ "Distance (nautical miles)",
      LI_MACHINE_PERFUSION_original_display ~ "Machine perfusion (OPTN coded)",
      LI_MACHINE_PERFUSION_FDA_era ~ "Machine perfusion (extrapolated)"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1, DRI ~ 2, DSRI ~ 2),
    missing = "no"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "kruskal.test",
    all_categorical() ~ "chisq.test"
  )) %>%
  modify_caption("**Table 1. Baseline characteristics of transplanted recipients by allocation pathway**") %>%
  modify_header(label = "**Characteristic**")

# The Rmd displays tbl1, tblS2 and tblS3; source() does not auto-print them.

# Save
tbl1_gt <- tbl1 %>% as_gt()
html_file <- tempfile(fileext = ".html")
gt::gtsave(tbl1_gt, filename = html_file)
tryCatch({
  webshot_to_tiff(
    url = paste0("file://", normalizePath(html_file)),
    path_tiff = fig_path("figures", "main", "Table1_Baseline_Characteristics.png"),
    vwidth = 1600, vheight = 2400, delay = 2, zoom = 2
  )
}, error = function(e) cat("webshot failed:", e$message, "\n"))

# Supplemental Table S2: DBD only
tblS2 <- timing_cohort %>%
  filter(donation_type == "DBD") %>%
  select(
    seq_flag, AGE, MELD_PELD_ADJ, GENDER,
    AGE_DON_donor, BMI_DON_CALC_donor, MACRO_FAT_CAT,
    max_BILIRUBIN, max_SGOT, max_SGPT, DRI, DSRI,
    COLD_ISCH_cat, DISTANCE,
    LI_MACHINE_PERFUSION_original_display, LI_MACHINE_PERFUSION_FDA_era
  ) %>%
  tbl_summary(
    by = seq_flag,
    label = list(
      AGE ~ "Recipient age (years)", MELD_PELD_ADJ ~ "Recipient MELD/PELD",
      GENDER ~ "Recipient gender", AGE_DON_donor ~ "Donor age (years)",
      BMI_DON_CALC_donor ~ "Donor BMI (kg/m\u00B2)",
      MACRO_FAT_CAT ~ "Donor macrosteatosis",
      max_BILIRUBIN ~ "Donor maximum bilirubin (mg/dL)",
      max_SGOT ~ "Donor maximum AST (U/L)", max_SGPT ~ "Donor maximum ALT (U/L)",
      DRI ~ "Donor Risk Index (DRI)", DSRI ~ "Discard Risk Index (DSRI)",
      COLD_ISCH_cat ~ "Cold ischemia time", DISTANCE ~ "Distance (nautical miles)",
      LI_MACHINE_PERFUSION_original_display ~ "Machine perfusion (OPTN coded)",
      LI_MACHINE_PERFUSION_FDA_era ~ "Machine perfusion (extrapolated)"
    ),
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})", all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1, DRI ~ 2, DSRI ~ 2),
    missing = "no"
  ) %>%
  add_p(test = list(all_continuous() ~ "kruskal.test", all_categorical() ~ "chisq.test")) %>%
  modify_caption("**Supplemental Table S2. Baseline characteristics of DBD transplanted recipients by allocation pathway**")

# Save Supplemental Table S2
tblS2_gt <- tblS2 %>% as_gt()
html_file_s2 <- tempfile(fileext = ".html")
gt::gtsave(tblS2_gt, filename = html_file_s2)
tryCatch({
  webshot_to_tiff(
    url = paste0("file://", normalizePath(html_file_s2)),
    path_tiff = fig_path("figures", "supplement", "TableS2_DBD_Baseline_Characteristics.png"),
    vwidth = 1600, vheight = 2400, delay = 2, zoom = 2
  )
}, error = function(e) cat("webshot S2 failed:", e$message, "\n"))

# Supplemental Table S3: DCD only
tblS3 <- timing_cohort %>%
  filter(donation_type == "DCD") %>%
  select(
    seq_flag, AGE, MELD_PELD_ADJ, GENDER,
    AGE_DON_donor, BMI_DON_CALC_donor, MACRO_FAT_CAT,
    max_BILIRUBIN, max_SGOT, max_SGPT, DRI, DSRI,
    COLD_ISCH_cat, DISTANCE,
    LI_MACHINE_PERFUSION_original_display, LI_MACHINE_PERFUSION_FDA_era
  ) %>%
  tbl_summary(
    by = seq_flag,
    label = list(
      AGE ~ "Recipient age (years)", MELD_PELD_ADJ ~ "Recipient MELD/PELD",
      GENDER ~ "Recipient gender", AGE_DON_donor ~ "Donor age (years)",
      BMI_DON_CALC_donor ~ "Donor BMI (kg/m\u00B2)",
      MACRO_FAT_CAT ~ "Donor macrosteatosis",
      max_BILIRUBIN ~ "Donor maximum bilirubin (mg/dL)",
      max_SGOT ~ "Donor maximum AST (U/L)", max_SGPT ~ "Donor maximum ALT (U/L)",
      DRI ~ "Donor Risk Index (DRI)", DSRI ~ "Discard Risk Index (DSRI)",
      COLD_ISCH_cat ~ "Cold ischemia time", DISTANCE ~ "Distance (nautical miles)",
      LI_MACHINE_PERFUSION_original_display ~ "Machine perfusion (OPTN coded)",
      LI_MACHINE_PERFUSION_FDA_era ~ "Machine perfusion (extrapolated)"
    ),
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})", all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1, DRI ~ 2, DSRI ~ 2),
    missing = "no"
  ) %>%
  add_p(test = list(all_continuous() ~ "kruskal.test", all_categorical() ~ "chisq.test")) %>%
  modify_caption("**Supplemental Table S3. Baseline characteristics of DCD transplanted recipients by allocation pathway**")

# Save Supplemental Table S3
tblS3_gt <- tblS3 %>% as_gt()
html_file_s3 <- tempfile(fileext = ".html")
gt::gtsave(tblS3_gt, filename = html_file_s3)
tryCatch({
  webshot_to_tiff(
    url = paste0("file://", normalizePath(html_file_s3)),
    path_tiff = fig_path("figures", "supplement", "TableS3_DCD_Baseline_Characteristics.png"),
    vwidth = 1600, vheight = 2400, delay = 2, zoom = 2
  )
}, error = function(e) cat("webshot S3 failed:", e$message, "\n"))
