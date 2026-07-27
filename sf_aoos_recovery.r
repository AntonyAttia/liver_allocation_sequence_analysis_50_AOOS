# sf_aoos_recovery.r
# Monthly AOOS and EP transplants alongside national liver recovery and
# transplantation. This is a descriptive figure; no pre/post test is performed.

library(tidyverse)
library(lubridate)
library(patchwork)

needed_donor_cols <- c("DON_DATE", "NUM_LI_RECOV", "NUM_LI_TX")
missing_cols <- setdiff(needed_donor_cols, names(donor_reduced))
if (length(missing_cols) > 0) {
  stop("donor_reduced missing columns: ", paste(missing_cols, collapse = ", "),
       ". Re-run 01_read_data.r to rebuild 01_raw_data_v2.rda.")
}

if (exists("timing_cohort") && is.data.frame(timing_cohort) && nrow(timing_cohort) > 0) {
  cohort_tx <- timing_cohort
  cohort_source <- "timing_cohort"
} else if (exists("final_cohort") && is.data.frame(final_cohort) && nrow(final_cohort) > 0) {
  cohort_tx <- final_cohort
  cohort_source <- "final_cohort"
} else {
  stop("Need timing_cohort or final_cohort in the environment.")
}

WINDOW_START <- as.Date("2024-01-01")
FEB_DIRECTIVE_DATE <- as.Date("2025-02-21")
OPTN_MEMO_DATE     <- as.Date("2025-08-05")

aoos_col  <- if (exists("GLOBAL_AOOS_COLOR")) GLOBAL_AOOS_COLOR else "#FF8C00"
ep_col    <- if (exists("GLOBAL_EP_COLOR")) GLOBAL_EP_COLOR else "#8A3FFC"
# Panel B uses two blues for print contrast.
recov_col <- "#00A6D6"
tx_col    <- "#00205B"

message("Preparing AOOS recovery figure from ", cohort_source)

donor_liver <- donor_reduced %>%
  mutate(DON_DATE = as.Date(DON_DATE)) %>%
  filter(!is.na(DON_DATE), NUM_LI_RECOV >= 1L, DON_DATE >= WINDOW_START)

latest_date <- max(donor_liver$DON_DATE, na.rm = TRUE)
message(
  "Donor span: ",
  format(min(donor_liver$DON_DATE)),
  " to ",
  format(latest_date)
)

monthly_recov <- donor_liver %>%
  mutate(period = floor_date(DON_DATE, "month")) %>%
  group_by(period) %>%
  summarise(
    livers_recovered = sum(NUM_LI_RECOV, na.rm = TRUE),
    livers_tx        = sum(NUM_LI_TX, na.rm = TRUE),
    .groups = "drop"
  )

monthly_pathway <- cohort_tx %>%
  mutate(
    tx_date = as.Date(TX_DATE),
    is_aoos = as.character(seq_flag) == "AOOS (out-of-sequence)",
    is_ep   = as.character(seq_flag) == "EP (Policy 9.10)"
  ) %>%
  filter(!is.na(tx_date), tx_date >= WINDOW_START) %>%
  mutate(period = floor_date(tx_date, "month")) %>%
  group_by(period) %>%
  summarise(
    n_aoos_tx = sum(is_aoos),
    n_ep_tx   = sum(is_ep),
    .groups = "drop"
  )

monthly <- full_join(monthly_recov, monthly_pathway, by = "period") %>%
  mutate(
    across(c(livers_recovered, livers_tx, n_aoos_tx, n_ep_tx), ~ replace_na(.x, 0L))
  ) %>%
  arrange(period)

x_lim <- c(min(monthly$period), latest_date)
x_scale <- scale_x_date(
  limits = x_lim,
  expand = expansion(mult = c(0.02, 0)),
  date_breaks = "3 months",
  date_labels = "%Y-%m"
)

vline_policy <- list(
  geom_vline(xintercept = FEB_DIRECTIVE_DATE, linetype = "dotted",
             color = "grey40", linewidth = 0.6),
  geom_vline(xintercept = OPTN_MEMO_DATE, linetype = "dashed",
             color = "grey30", linewidth = 0.6)
)

top_long <- monthly %>%
  pivot_longer(c(n_aoos_tx, n_ep_tx), names_to = "metric", values_to = "n") %>%
  mutate(metric = factor(
    dplyr::case_match(
      metric,
      "n_aoos_tx" ~ "AOOS transplants",
      "n_ep_tx" ~ "EP transplants"
    ),
    levels = c("AOOS transplants", "EP transplants")
  ))

p_top <- ggplot(top_long, aes(period, n, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  vline_policy +
  coord_cartesian(ylim = c(0, NA)) +
  x_scale +
  scale_color_manual(values = c(
    "AOOS transplants" = aoos_col,
    "EP transplants" = ep_col
  )) +
  labs(
    title = "A. AOOS and EP transplants / month (study cohort)",
    x = NULL, y = "Transplants", color = NULL
  ) +
  theme_main(11) +
  theme(legend.position = "bottom")

recov_long <- monthly %>%
  pivot_longer(c(livers_recovered, livers_tx), names_to = "metric", values_to = "n") %>%
  mutate(metric = dplyr::case_match(
    metric,
    "livers_recovered" ~ "Livers recovered",
    "livers_tx" ~ "Livers transplanted"
  ))

p_recov <- ggplot(recov_long, aes(period, n, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  vline_policy +
  coord_cartesian(ylim = c(0, NA)) +
  x_scale +
  scale_color_manual(values = c(
    "Livers recovered" = recov_col,
    "Livers transplanted" = tx_col
  )) +
  labs(
    title = "B. National liver recovery and transplantation / month",
    x = NULL, y = "Count per month", color = NULL
  ) +
  theme_main(11) +
  theme(legend.position = "bottom")

aoos_recovery_figure <- (p_top / p_recov) +
  plot_annotation(
    caption = "Dotted line = 21 Feb 2025 HRSA directive; dashed line = 5 Aug 2025 OPTN memo."
  ) &
  theme(plot.caption = element_text(hjust = 0, size = 10))
print(aoos_recovery_figure)

out_path <- fig_path("figures", "supplement", "S3_AOOS_Recovery.png")
ggsave(out_path, aoos_recovery_figure,
       width = 12, height = 8, dpi = FIG_DPI, compression = "lzw")
message("Saved: ", out_path)
