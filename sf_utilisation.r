# sf_utilisation.r — Supplemental Utilisation Figures (separate panels per pathway)
# Stratified by DCD/DBD. Each pathway shown individually.
# Moved to supplement with undercount caveat.

EP_BYPASS_CODE    <- 887L
AOOS_BYPASS_CODES <- c(861L, 862L, 863L, 888L)
AOOS_TEXT_CODES   <- c(798L, 799L, 898L)
AOOS_TEXT_PATTERN <- regex("expedited|open offer|aggressive", ignore_case = TRUE)

# Accepted organs with timing and pathway
accepted_util <- timing_cohort %>%
  distinct(DONOR_ID, MATCH_ID_CODE, diff_final_clamp, seq_flag, donation_type) %>%
  mutate(hours_from_clamp = diff_final_clamp, status = "Utilized") %>%
  filter(!is.na(donation_type))

# Discarded organs (from 03)
discarded_util <- discarded_complete %>%
  mutate(
    seq_flag_disc = case_when(
      allocation_type == "AOOS (out-of-sequence)" ~ "AOOS (out-of-sequence)",
      allocation_type == "EP (Policy 9.10)"       ~ "EP (Policy 9.10)",
      TRUE ~ "Sequential"
    )
  )

cutoff_hours <- seq(-24, 24, by = 2)

# For each pathway, calculate utilization at each cutoff
# In-sequence: all transplants placed at <50 or >=50 sequence + all Sequential discards (one denominator)
calculate_pathway_util <- function(pathway_label, disc_label, dt_filter = NULL) {
  accepted_sub <- accepted_util %>%
    filter(
      pathway_label == "In-sequence" & seq_flag %in% c("<50 sequence", "\u226550 sequence") |
        as.character(seq_flag) == pathway_label
    )
  if (!is.null(dt_filter)) accepted_sub <- accepted_sub %>% filter(donation_type == dt_filter)

  disc_sub <- discarded_util %>% filter(seq_flag_disc == disc_label)
  if (!is.null(dt_filter)) disc_sub <- disc_sub %>% filter(donor_type == dt_filter)

  n_disc <- nrow(disc_sub)

  lapply(cutoff_hours, function(hr) {
    n_util <- sum(accepted_sub$hours_from_clamp >= hr, na.rm = TRUE)
    total <- n_util + n_disc
    tibble(cutoff_hours = hr, pathway = pathway_label,
           donation_type = ifelse(is.null(dt_filter), "All", dt_filter),
           n_utilized = n_util, n_discarded = n_disc, total = total,
           pct_utilized = if_else(total > 0, round(100 * n_util / total, 1), 0))
  }) %>% bind_rows()
}

# Three pathways: In-sequence (pooled <50 + >=50), AOOS, EP. Discards: Sequential / AOOS / EP.
pathways <- c("In-sequence", "AOOS (out-of-sequence)", "EP (Policy 9.10)")
disc_labels <- c("Sequential", "AOOS (out-of-sequence)", "EP (Policy 9.10)")
donor_types <- c("DBD", "DCD")

all_util_data <- list()
for (i in seq_along(pathways)) {
  for (dt in donor_types) {
    all_util_data[[length(all_util_data) + 1]] <- calculate_pathway_util(pathways[i], disc_labels[i], dt)
  }
}
util_data <- bind_rows(all_util_data)

# One plot per pathway (DBD vs DCD lines)
make_util_plot <- function(data, pathway_name) {
  pd <- data %>% filter(pathway == pathway_name)
  ggplot(pd, aes(x = cutoff_hours, y = pct_utilized, color = donation_type, linetype = donation_type)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_x_continuous(breaks = seq(-24, 24, by = 4),
                       labels = function(x) ifelse(x == 0, "0", ifelse(x < 0, as.character(x), paste0("+", x)))) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_color_manual(values = c("DBD" = "#003366", "DCD" = "#8B4513")) +
    labs(title = pathway_name,
         x = "Hours from cross-clamp", y = "Utilization rate (%)",
         color = "Donor Type", linetype = "Donor Type") +
    theme_main(base_size = 11) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
}

p_inseq <- make_util_plot(util_data, "In-sequence")
p_aoos  <- make_util_plot(util_data, "AOOS (out-of-sequence)")
p_ep    <- make_util_plot(util_data, "EP (Policy 9.10)")

# Combined as one figure with panels A B C (three pathways); y-axis "Utilization rate (%)" is self-explanatory
combined_util <- (p_inseq | p_aoos | p_ep) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
combined_util <- combined_util + plot_annotation(tag_levels = "A")

print(combined_util)
# Wide aspect ratio; slightly shorter height so figure is less tall
ggsave(here("figures", "supplement", "SF_Utilisation_Combined.png"), combined_util,
       width = 14, height = 4.5, dpi = 300)

# Utilisation summary table (selected time points); order In-sequence, AOOS, EP
pathway_order <- c("In-sequence", "AOOS (out-of-sequence)", "EP (Policy 9.10)")
util_summary <- util_data %>%
  filter(cutoff_hours %in% c(-12, -8, -4, 0, 4)) %>%
  select(pathway, donation_type, cutoff_hours, n_utilized, pct_utilized, total) %>%
  mutate(pathway = factor(pathway, levels = pathway_order)) %>%
  arrange(pathway, donation_type, cutoff_hours)

cat("\n=== Supplemental Utilisation Table (selected time points) ===\n")
cat("CAVEAT: AOOS attempts on non-transplanted livers are significantly undercounted.\n")
cat("These numbers should be interpreted with extreme caution.\n\n")
print(knitr::kable(util_summary, format = "simple"))

# Overall utilization rate and total counts by pathway (for reporting)
# (A) At 24-hour window: utilized = acceptances with hours_from_clamp >= -24 (used for curves)
overall_util <- util_data %>%
  filter(cutoff_hours == -24) %>%
  group_by(pathway) %>%
  summarise(n_utilized = sum(n_utilized), n_discarded = sum(n_discarded), .groups = "drop") %>%
  mutate(total = n_utilized + n_discarded,
         utilization_pct = if_else(total > 0, round(100 * n_utilized / total, 1), NA_real_),
         pathway = factor(pathway, levels = pathway_order)) %>%
  arrange(pathway)

cat("\n=== Overall utilization by pathway (at 24h window: utilized = acceptance >= -24h from clamp) ===\n")
cat("Pathway          | Utilized (transplanted) | Discarded | Total | Utilization (%)\n")
cat("-----------------|-------------------------|-----------|-------|------------------\n")
pathway_display <- c("In-sequence" = "In-sequence",
                     "AOOS (out-of-sequence)" = "AOOS", "EP (Policy 9.10)" = "EP")
for (i in seq_len(nrow(overall_util))) {
  r <- overall_util[i, ]
  lab <- pathway_display[as.character(r$pathway)]
  if (is.na(lab)) lab <- as.character(r$pathway)
  cat(sprintf("%-16s | %21s | %9s | %5s | %s%%\n",
              lab, format(r$n_utilized, big.mark = ","),
              format(r$n_discarded, big.mark = ","),
              format(r$total, big.mark = ","),
              format(r$utilization_pct, nsmall = 1)))
}
cat("\n")
print(overall_util)

# (B) Whole period: utilized = all accepted transplants (no time cutoff); aligns with cohort n
accepted_by_pathway <- accepted_util %>%
  mutate(pathway = case_when(
    seq_flag %in% c("<50 sequence", "\u226550 sequence") ~ "In-sequence",
    TRUE ~ as.character(seq_flag)
  )) %>%
  group_by(pathway) %>%
  summarise(n_utilized = n(), .groups = "drop")
discarded_by_pathway <- discarded_util %>%
  mutate(pathway = case_when(
    seq_flag_disc == "Sequential" ~ "In-sequence",
    TRUE ~ seq_flag_disc
  )) %>%
  group_by(pathway) %>%
  summarise(n_discarded = n(), .groups = "drop")
overall_util_full <- accepted_by_pathway %>%
  left_join(discarded_by_pathway, by = "pathway") %>%
  mutate(
    n_discarded = coalesce(n_discarded, 0L),
    total = n_utilized + n_discarded,
    utilization_pct = if_else(total > 0, round(100 * n_utilized / total, 1), NA_real_),
    pathway = factor(pathway, levels = pathway_order)
  ) %>%
  arrange(pathway)

cat("\n=== Overall utilization by pathway (whole period: all utilized = all transplanted) ===\n")
cat("Pathway          | Utilized (transplanted) | Discarded | Total | Utilization (%)\n")
cat("-----------------|-------------------------|-----------|-------|------------------\n")
for (i in seq_len(nrow(overall_util_full))) {
  r <- overall_util_full[i, ]
  lab <- pathway_display[as.character(r$pathway)]
  if (is.na(lab)) lab <- as.character(r$pathway)
  cat(sprintf("%-16s | %21s | %9s | %5s | %s%%\n",
              lab, format(r$n_utilized, big.mark = ","),
              format(r$n_discarded, big.mark = ","),
              format(r$total, big.mark = ","),
              format(r$utilization_pct, nsmall = 1)))
}
cat("\n")
print(overall_util_full)

# Diagnostic: why overall_util shows 5,430 AOOS utilized vs 5,486 in timing_cohort
# overall_util is computed at cutoff_hours == -24, so n_utilized = count with hours_from_clamp >= -24.
# Acceptances *more than* 24h before clamp (hours_from_clamp < -24) are therefore excluded.
n_aoos_timing   <- sum(timing_cohort$seq_flag == "AOOS (out-of-sequence)", na.rm = TRUE)
n_aoos_util     <- accepted_util %>% filter(seq_flag == "AOOS (out-of-sequence)") %>% nrow()
aoos_accepted   <- accepted_util %>% filter(seq_flag == "AOOS (out-of-sequence)")
n_aoos_at_cutoff <- sum(aoos_accepted$hours_from_clamp >= -24, na.rm = TRUE)
n_aoos_before_24 <- sum(aoos_accepted$hours_from_clamp < -24, na.rm = TRUE)
n_aoos_na_hr    <- sum(is.na(aoos_accepted$hours_from_clamp))
cat("\n=== AOOS count diagnostic (timing cohort vs overall_util table) ===\n")
cat(sprintf("AOOS in timing_cohort: %d; in accepted_util: %d\n", n_aoos_timing, n_aoos_util))
cat(sprintf("Overall util table uses cutoff_hours = -24 => n_utilized = count with hours_from_clamp >= -24.\n"))
cat(sprintf("  AOOS with hours_from_clamp >= -24 (counted): %d\n", n_aoos_at_cutoff))
cat(sprintf("  AOOS with hours_from_clamp <  -24 (excluded): %d (accepted >24h before clamp)\n", n_aoos_before_24))
cat(sprintf("  AOOS with NA hours_from_clamp: %d\n", n_aoos_na_hr))
cat(sprintf("  => Table 'Utilized' = %d; cohort n AOOS = %d; difference = %d\n",
            n_aoos_at_cutoff, n_aoos_timing, n_aoos_timing - n_aoos_at_cutoff))
