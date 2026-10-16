# 00_global_config.r
# Shared theme, colours, anonymisation, timing cohort, and helpers.
# Sourced once by the main Rmd; figure scripts depend on objects created here.

# Figure export (publication TIFF)
FIG_EXT <- ".tiff"
FIG_DPI <- 300L

# Build a here() path and force the publication TIFF extension.
fig_path <- function(...) {
  args <- list(...)
  n <- length(args)
  args[[n]] <- paste0(tools::file_path_sans_ext(as.character(args[[n]])), FIG_EXT)
  do.call(here::here, args)
}

# Convert a PNG file to TIFF (used after webshot / rsvg / gtsave).
png_file_to_tiff <- function(path_png, path_tiff) {
  if (!requireNamespace("png", quietly = TRUE) || !requireNamespace("tiff", quietly = TRUE)) {
    stop("Install packages 'png' and 'tiff' for PNG to TIFF conversion (tables, cohort flow).")
  }
  arr <- png::readPNG(path_png, native = FALSE)
  tiff::writeTIFF(arr, path_tiff, compression = "LZW")
}

# webshot2 captures to PNG; convert to TIFF for journal submission.
webshot_to_tiff <- function(url, path_tiff, ...) {
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  webshot2::webshot(url = url, file = tmp, ...)
  png_file_to_tiff(tmp, path_tiff)
}

# Colour palette (4 groups)
SEQ_COLORS_4 <- c(
  "<50 sequence"           = "#0066CC",
  "\u226550 sequence"      = "#00468B",
  "AOOS (out-of-sequence)" = "#FF8C00",
  "EP (Policy 9.10)"       = "#8A3FFC"
)

GLOBAL_SEQUENTIAL_COLOR <- "#00468B"
GLOBAL_AOOS_COLOR       <- "#FF8C00"
GLOBAL_EP_COLOR         <- "#8A3FFC"

GLOBAL_SEQUENCE_COLORS <- c(
  "Sequential"              = GLOBAL_SEQUENTIAL_COLOR,
  "AOOS (out-of-sequence)"  = GLOBAL_AOOS_COLOR,
  "EP (Policy 9.10)"        = GLOBAL_EP_COLOR
)

EP_POLICY_DATE <- as.Date("2021-03-25")

# Plot theme
theme_main <- function(base_size = 12) {
  theme_minimal(base_family = "Helvetica", base_size = base_size) +
    theme(
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.line         = element_line(color = "black", linewidth = 1),
      axis.title        = element_text(face = "bold", color = "black"),
      axis.title.x      = element_text(face = "bold", color = "black", margin = margin(t = 10)),
      axis.title.y      = element_text(face = "bold", color = "black", margin = margin(r = 10)),
      axis.text         = element_text(color = "black", family = "Helvetica"),
      axis.ticks        = element_line(color = "black", linewidth = 0.5),
      panel.border      = element_blank(),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      plot.title        = element_text(face = "bold", family = "Helvetica", hjust = 0.5),
      legend.text       = element_text(family = "Helvetica"),
      legend.title      = element_text(face = "bold", family = "Helvetica")
    )
}
theme_set(theme_main())

format_pvalue <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("P < .001")
  if (p < 0.01)  return("P < .01")
  if (p < 0.05)  return("P < .05")
  paste0("P = ", sub("^0\\.", ".", sprintf("%.3f", p)))
}

# Anonymise OPO and center codes
set.seed(42)

all_opo_codes <- final_cohort %>%
  mutate(OPO_CTR_CODE_clean = case_when(
    OPO_CTR_CODE == "00527" ~ "11036",
    OPO_CTR_CODE == "05239" ~ "14353",
    OPO_CTR_CODE == "00093" ~ "19282",
    TRUE ~ OPO_CTR_CODE
  )) %>%
  filter(!is.na(OPO_CTR_CODE_clean)) %>%
  pull(OPO_CTR_CODE_clean) %>%
  unique()

all_center_codes <- final_cohort %>%
  filter(!is.na(LISTING_CTR_CODE)) %>%
  pull(LISTING_CTR_CODE) %>%
  unique()

opo_anon_mapping <- tibble(OPO_CTR_CODE_real = all_opo_codes) %>%
  mutate(random_order = sample(n()), OPO_anon = paste0("OPO ", random_order)) %>%
  arrange(random_order) %>%
  select(OPO_CTR_CODE_real, OPO_anon)

center_anon_mapping <- tibble(LISTING_CTR_CODE_real = all_center_codes) %>%
  mutate(random_order = sample(n()), Center_anon = paste0("Center ", random_order)) %>%
  arrange(random_order) %>%
  select(LISTING_CTR_CODE_real, Center_anon)

anonymize_opo <- function(code) {
  if (is.na(code)) return(NA_character_)
  code_clean <- case_when(
    code == "00527" ~ "11036", code == "05239" ~ "14353",
    code == "00093" ~ "19282", TRUE ~ code
  )
  idx <- match(code_clean, opo_anon_mapping$OPO_CTR_CODE_real)
  if (is.na(idx)) return(paste0("OPO_", code))
  opo_anon_mapping$OPO_anon[idx]
}
anonymize_center <- function(code) {
  if (is.na(code)) return(NA_character_)
  idx <- match(code, center_anon_mapping$LISTING_CTR_CODE_real)
  if (is.na(idx)) return(paste0("Center_", code))
  center_anon_mapping$Center_anon[idx]
}
anonymize_opo_vec    <- Vectorize(anonymize_opo)
anonymize_center_vec <- Vectorize(anonymize_center)

# Timing cohort (plausible timing subset used for most figures)
n_before_timing <- nrow(final_cohort)

timing_cohort_temp <- final_cohort %>%
  mutate(placement_time_hours = case_when(
    is.na(earliest_init)          ~ NA_real_,
    is.na(FINAL_RESPONSE_DT_tx)  ~ NA_real_,
    TRUE ~ as.numeric(difftime(FINAL_RESPONSE_DT_tx, earliest_init, units = "hours"))
  ))

timing_cohort <- timing_cohort_temp %>%
  filter(
    !is.na(FINAL_RESPONSE_DT_tx),
    as.Date(FINAL_RESPONSE_DT_tx) <= TX_DATE,
    !is.na(placement_time_hours),
    !is.na(diff_final_clamp),
    earliest_init <= FINAL_RESPONSE_DT_tx
  ) %>%
  mutate(
    seq_flag = factor(seq_flag, levels = c(
      "<50 sequence", "\u226550 sequence", "AOOS (out-of-sequence)", "EP (Policy 9.10)"
    )),
    placement_time_hours = as.numeric(difftime(FINAL_RESPONSE_DT_tx, earliest_init, units = "hours"))
  )

cat("Timing cohort: ", comma(nrow(timing_cohort)),
    " (from ", comma(n_before_timing), ")\n", sep = "")
print(table(timing_cohort$seq_flag))

timing_cohort <- timing_cohort %>%
  mutate(
    donation_type = case_when(
      NON_HRT_DON_donor %in% c("Y","YES","1","TRUE") ~ "DCD",
      NON_HRT_DON_donor %in% c("N","NO","0","FALSE") ~ "DBD",
      TRUE ~ NA_character_
    ),
    MACRO_FAT_CAT = {
      v <- suppressWarnings(as.numeric(MACRO_FAT_LI_DON_donor))
      case_when(
        is.na(v) ~ NA_character_,
        v == 0 ~ "0%", v > 0 & v < 30 ~ "1-29%",
        v >= 30 & v < 60 ~ "30-59%", v >= 60 ~ "\u226560%"
      ) %>% factor(levels = c("0%", "1-29%", "30-59%", "\u226560%"))
    },
    COLD_ISCH_cat = case_when(
      is.na(COLD_ISCH_num) ~ NA_character_,
      COLD_ISCH_num < 8  ~ "<8 hours",
      COLD_ISCH_num >= 8 ~ "\u22658 hours",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("<8 hours", "\u22658 hours")),
    LI_MACHINE_PERFUSION_original_display = case_when(
      LI_MACHINE_PERFUSION_original == "Y" ~ "Yes",
      LI_MACHINE_PERFUSION_original == "N" ~ "No",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Yes", "No"))
  )

# OPO codes for maps (mapping file not included)
opo_mapping_file <- Sys.getenv("AOOS_OPO_MAPPING_FILE")
opo_direct_mapping <- if (nzchar(opo_mapping_file)) read.csv(opo_mapping_file, colClasses = "character") else NULL

cat("Global config loaded.\n")
