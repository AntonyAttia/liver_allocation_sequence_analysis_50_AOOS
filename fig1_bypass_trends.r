# fig1_bypass_trends.r: Figure 1. Bypass code usage over time (full PTR)

library(lubridate)

EP_BYPASS_CODE    <- 887L
AOOS_BYPASS_CODES <- c(861L, 862L, 863L, 888L)
AOOS_TEXT_CODES   <- c(798L, 799L, 898L)
AOOS_TEXT_PATTERN <- regex("expedited|open offer|aggressive", ignore_case = TRUE)

# All bypass rows from full PTR (not limited to study period)
bypass_rows <- ptr_reduced %>%
  filter(OFFER_ACCEPT == "B", !is.na(MATCH_DATE)) %>%
  mutate(match_date = as.Date(MATCH_DATE),
         month = floor_date(match_date, "month"))

# Text-matched rows (N or B with text codes + pattern)
text_rows <- ptr_reduced %>%
  filter(OFFER_ACCEPT %in% c("B", "N"),
         PRIME_OPO_REFUSAL_ID %in% AOOS_TEXT_CODES,
         str_detect(replace_na(OPO_OSTXT_PRIMARY, ""), AOOS_TEXT_PATTERN),
         !is.na(MATCH_DATE)) %>%
  mutate(match_date = as.Date(MATCH_DATE),
         month = floor_date(match_date, "month"))

# Bypass codes (distinct match-runs per month per code)
bypass_coded <- bypass_rows %>%
  filter(PRIME_OPO_REFUSAL_ID %in% c(EP_BYPASS_CODE, AOOS_BYPASS_CODES)) %>%
  mutate(code = as.character(PRIME_OPO_REFUSAL_ID)) %>%
  distinct(DONOR_ID, MATCH_ID_CODE, code, month)

# Text-matched codes
text_coded <- text_rows %>%
  mutate(code = as.character(PRIME_OPO_REFUSAL_ID)) %>%
  distinct(DONOR_ID, MATCH_ID_CODE, code, month)

all_coded_raw <- bind_rows(bypass_coded, text_coded) %>%
  distinct(DONOR_ID, MATCH_ID_CODE, code, month)

# Legend labels, ordering, and palette (exclude zero-count codes)
expected_codes <- c("887", "861", "862", "863", "799", "888", "798", "898")

code_totals <- all_coded_raw %>%
  count(code, name = "total") %>%
  arrange(desc(total))

codes_present <- code_totals %>%
  filter(total > 0) %>%
  pull(code)

codes_missing <- setdiff(expected_codes, codes_present)
if (length(codes_missing) > 0) {
  cat("Figure 1: codes with 0 match-runs (excluded from plot): ",
      paste(codes_missing, collapse = ", "), "\n", sep = "")
}

ordered_present_codes <- expected_codes[expected_codes %in% codes_present]

code_labels_real <- c(
  "887" = "[EP] 887 Not offered \u2013 expedited placement",
  "861" = "[OPTN] 861 Operational \u2013 OPO",
  "862" = "[OPTN] 862 Donor medical urgency",
  "863" = "[OPTN] 863 Expedited placement attempt",
  "799" = "[OPTN] 799 Other + keyword text match",
  "888" = "[Lit] 888 Other bypass",
  "798" = "[Lit] 798 + keyword text match",
  "898" = "[Lit] 898 + keyword text match"
)
code_labels_real <- code_labels_real[ordered_present_codes]

optn_codes <- c("861", "862", "863", "799")
lit_codes  <- c("888", "798", "898")

optn_present <- ordered_present_codes[ordered_present_codes %in% optn_codes]
lit_present  <- ordered_present_codes[ordered_present_codes %in% lit_codes]

legend_rows <- max(1L, length(optn_present), length(lit_present))

.spacer_labels <- function(prefix, n, space_char) {
  if (n <= 0) return(character(0))
  ids <- paste0(prefix, seq_len(n))
  lbl <- vapply(seq_len(n), function(i) strrep(space_char, i), character(1))
  stats::setNames(lbl, ids)
}

code_labels <- c(
  stats::setNames(code_labels_real["887"], "887"),
  .spacer_labels("_spacer_ep_",   legend_rows - 1L,                  "\u2007"),
  code_labels_real[optn_present],
  .spacer_labels("_spacer_optn_", legend_rows - length(optn_present), "\u2008"),
  code_labels_real[lit_present],
  .spacer_labels("_spacer_lit_",  legend_rows - length(lit_present),  "\u2009")
)

color_real <- c(
  "887" = "#D62728",
  "861" = "#1F77B4",
  "862" = "#2CA02C",
  "863" = "#9467BD",
  "799" = "#17BECF",
  "888" = "#FF7F0E",
  "798" = "#8C564B",
  "898" = "#E377C2"
)

color_real <- color_real[names(code_labels_real)]

real_label_strings   <- unname(code_labels[!grepl("^_spacer_", names(code_labels))])
spacer_label_strings <- unname(code_labels[grepl("^_spacer_", names(code_labels))])

code_palette <- c(
  stats::setNames(unname(color_real), real_label_strings),
  stats::setNames(rep("transparent", length(spacer_label_strings)),
                  spacer_label_strings)
)

EP_POLICY_COLOR <- "#D62728"
OPTN_COMPLIANCE_DATE <- as.Date("2025-08-01")
OPTN_COMPLIANCE_COLOR <- "#333333"

all_coded <- all_coded_raw %>%
  filter(code %in% names(code_labels_real)) %>%
  mutate(code_label = factor(code_labels_real[code],
                             levels = unname(code_labels)))

# Monthly counts (full granularity for tables and reporting)
timeline_data <- all_coded %>%
  count(month, code, code_label, name = "n_match_runs")

# Figure 1 only: aggregate into 3 series (EP, OPTN AOOS definition, Lit)
ep_codes_fig   <- "887"
optn_codes_fig <- c("861", "862", "863", "799")
lit_codes_fig  <- c("798", "898", "888")
timeline_3line <- all_coded %>%
  mutate(series = case_when(
    code %in% ep_codes_fig   ~ "EP",
    code %in% optn_codes_fig ~ "OPTN AOOS",
    code %in% lit_codes_fig  ~ "Lit",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(series)) %>%
  group_by(month, series) %>%
  summarise(n_match_runs = n_distinct(DONOR_ID, MATCH_ID_CODE), .groups = "drop")
month_seq <- seq(min(timeline_3line$month, na.rm = TRUE),
                 max(timeline_3line$month, na.rm = TRUE), by = "month")
timeline_3line <- timeline_3line %>%
  complete(month = month_seq, series = c("EP", "OPTN AOOS", "Lit"), fill = list(n_match_runs = 0)) %>%
  mutate(series = factor(series, levels = c("EP", "OPTN AOOS", "Lit")))
fig1_series_colors <- c("EP" = "#D62728", "OPTN AOOS" = "#1F77B4", "Lit" = "#FF7F0E")

# Printed output for results
monthly_wide <- all_coded %>%
  count(month, code, name = "n_match_runs") %>%
  arrange(month, code) %>%
  tidyr::pivot_wider(names_from = code, values_from = n_match_runs, values_fill = 0)
print(monthly_wide, n = Inf)

all_time_by_code <- all_coded %>%
  count(code, name = "total") %>%
  mutate(prop_pct = round(100 * total / sum(total), 1))
total_code_occurrences <- sum(all_time_by_code$total)
distinct_runs_all <- all_coded_raw %>%
  filter(code %in% names(code_labels_real)) %>%
  distinct(DONOR_ID, MATCH_ID_CODE, month) %>%
  nrow()
cat("Distinct match runs (any of these codes) over all time: ", distinct_runs_all, "\n", sep = "")
cat("Total code occurrences (match run × code): ", total_code_occurrences, "\n", sep = "")
print(all_time_by_code %>% mutate(label = unname(code_labels_real[code])) %>% select(code, label, total, prop_pct), n = Inf)

ep_policy_month <- floor_date(EP_POLICY_DATE, "month")
all_coded_post <- all_coded %>% filter(month >= ep_policy_month)
post_by_code <- all_coded_post %>%
  count(code, name = "total") %>%
  mutate(prop_pct = round(100 * total / sum(total), 1))
total_post <- sum(post_by_code$total)
distinct_runs_post <- all_coded_raw %>%
  filter(code %in% names(code_labels_real), month >= ep_policy_month) %>%
  distinct(DONOR_ID, MATCH_ID_CODE, month) %>%
  nrow()
cat("Distinct match runs (any of these codes) after policy: ", distinct_runs_post, "\n", sep = "")
cat("Total code occurrences after policy: ", total_post, "\n", sep = "")
print(post_by_code %>% mutate(label = unname(code_labels_real[code])) %>% select(code, label, total, prop_pct), n = Inf)

optn_compliance_month <- floor_date(OPTN_COMPLIANCE_DATE, "month")
all_coded_optn <- all_coded %>% filter(month >= optn_compliance_month)
optn_by_code <- all_coded_optn %>%
  count(code, name = "total") %>%
  mutate(prop_pct = if (sum(total) > 0) round(100 * total / sum(total), 1) else NA_real_)
total_optn <- sum(optn_by_code$total, na.rm = TRUE)
distinct_runs_optn <- all_coded_raw %>%
  filter(code %in% names(code_labels_real), month >= optn_compliance_month) %>%
  distinct(DONOR_ID, MATCH_ID_CODE, month) %>%
  nrow()

# Complete month-by-code grid for granular timeline (tables/reporting)
month_seq <- seq(min(timeline_data$month, na.rm = TRUE),
                 max(timeline_data$month, na.rm = TRUE), by = "month")
timeline_data <- timeline_data %>%
  select(month, code_label, n_match_runs) %>%
  complete(month = month_seq, code_label, fill = list(n_match_runs = 0))

# Figure 1: 3-line plot (EP, OPTN AOOS definition, Lit)
fig1_x_range <- range(timeline_3line$month, na.rm = TRUE)
fig1_ep_line_x   <- EP_POLICY_DATE - 21L
fig1_optn_line_x <- OPTN_COMPLIANCE_DATE - 21L
p_fig1 <- timeline_3line %>%
  ggplot(aes(x = month, y = n_match_runs, color = series)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = fig1_ep_line_x, linetype = "dotted", color = EP_POLICY_COLOR, linewidth = 0.8) +
  annotate("text", x = fig1_ep_line_x, y = Inf,
           label = "EP Policy 9.10\n03/2021", vjust = 1.5, hjust = -0.1, color = EP_POLICY_COLOR, size = 3) +
  geom_vline(xintercept = fig1_optn_line_x, linetype = "dotted", color = OPTN_COMPLIANCE_COLOR, linewidth = 0.8) +
  annotate("text", x = fig1_optn_line_x, y = Inf,
           label = "OPTN warning\n08/2025", vjust = 1.5, hjust = -0.1, color = OPTN_COMPLIANCE_COLOR, size = 3) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",
              expand = expansion(mult = c(0.02, 0)), limits = c(fig1_x_range[1], as.Date("2026-12-01"))) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(
    values = fig1_series_colors,
    labels = c(
      "EP" = "Expedited placement (EP, code 887)",
      "OPTN AOOS" = "OPTN AOOS definition (codes 861, 862, 863, 799)",
      "Lit" = "Literature-based AOOS definition (798, 898, 888)"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Bypass/Refusal Code Usage Over Time",
    subtitle = "Distinct match runs per month by definitional group",
    x = "Month", y = "Distinct match runs per month", color = NULL
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(linewidth = 1.2))) +
  coord_cartesian(clip = "off") +
  theme_main(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.key.width = grid::unit(18, "pt"),
    plot.margin = margin(t = 6, r = 8, b = 4, l = 4, unit = "mm")
  )

print(p_fig1)
dir.create(here("figures", "main"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("figures", "supplement"), showWarnings = FALSE, recursive = TRUE)
ggsave(fig_path("figures", "main", "Figure1_Bypass_Trends.png"), p_fig1, width = 14, height = 8, dpi = FIG_DPI, compression = "lzw")

# Definition text for figure caption and supplement (OPTN vs Lit)
fig1_definitions_text <- paste(
  "EP: OPTN code 887 (not offered \u2013 expedited placement).",
  "OPTN AOOS: OPTN-defined codes 861, 862, 863, and 799 (with keyword text).",
  "Lit: literature-based codes 798 and 898 (with keyword text) and 888; see supplement for code-level detail."
)
# For use in Rmd: fig1_definitions_text

# Supplemental Table S1: bypass/refusal code summary (all-time and after EP policy only)
# Row order: summary rows, then EP (887), then all OPTN (861, 862, 863, 799), then all Lit (798, 898, 888)
tab_s1_code_order <- c("887", "861", "862", "863", "799", "798", "898", "888")
tab_s1_bypass_summary <- all_time_by_code %>%
  mutate(label = unname(code_labels_real[code])) %>%
  left_join(post_by_code %>% select(code, total_post = total, prop_pct_post = prop_pct), by = "code") %>%
  select(code, label, total_all = total, prop_pct_all = prop_pct, total_post, prop_pct_post) %>%
  filter(code %in% tab_s1_code_order) %>%
  arrange(match(code, tab_s1_code_order)) %>%
  select(-code)

tab_s1_summary_rows <- tibble(
  label = c("Distinct match runs (any code)", "Total code occurrences (match run × code)"),
  total_all = c(distinct_runs_all, total_code_occurrences),
  prop_pct_all = c(NA_real_, NA_real_),
  total_post = c(distinct_runs_post, total_post),
  prop_pct_post = c(NA_real_, NA_real_)
)

tab_s1_table_df <- bind_rows(tab_s1_summary_rows, tab_s1_bypass_summary)

tab_s1_gt <- tab_s1_table_df %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Supplemental Table S1. Bypass/refusal code summary",
    subtitle = "All-time vs after EP policy (≥ Mar 2021)"
  ) %>%
  gt::tab_spanner(label = "All-time", columns = c(total_all, prop_pct_all)) %>%
  gt::tab_spanner(label = "After EP policy (≥ Mar 2021)", columns = c(total_post, prop_pct_post)) %>%
  gt::cols_label(
    label = "Metric / Bypass code",
    total_all = "n",
    prop_pct_all = "%",
    total_post = "n",
    prop_pct_post = "%"
  ) %>%
  gt::fmt_number(columns = c(total_all, total_post), decimals = 0) %>%
  gt::fmt_number(columns = c(prop_pct_all, prop_pct_post), decimals = 1) %>%
  gt::sub_missing(columns = c(prop_pct_all, prop_pct_post), missing_text = "-") %>%
  gt::tab_style(style = gt::cell_text(weight = "bold"), locations = gt::cells_body(rows = 1:2, columns = 1)) %>%
  gt::tab_options(table.font.names = "Helvetica", table.font.size = 14) %>%
  gt::tab_source_note(
    gt::md("**Definitions:** EP = code 887. OPTN AOOS = OPTN-defined codes 861, 862, 863, 799. Lit = literature-based 798, 898, 888. **+ text:** codes 798, 799, 898 require matching refusal/offer text for keywords: \"expedited\", \"open offer\", \"aggressive\".")
  )
tab_s1_html <- tempfile(fileext = ".html")
gt::gtsave(tab_s1_gt, tab_s1_html)
tryCatch({
  webshot_to_tiff(
    url = paste0("file://", normalizePath(tab_s1_html)),
    path_tiff = fig_path("figures", "supplement", "TableS1_Bypass_Summary.png"),
    vwidth = 1600, vheight = 900, delay = 2, zoom = 2
  )
}, error = function(e) cat("webshot Table S1 failed:", e$message, "\n"))
readr::write_csv(tab_s1_table_df, here("figures", "supplement", "TableS1_Bypass_Summary.csv"))
