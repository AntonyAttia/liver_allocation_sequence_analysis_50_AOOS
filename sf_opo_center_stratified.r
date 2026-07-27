# sf_opo_center_stratified.r
# Supplement Figure 7. OPO-center pair analysis stratified by DBD and DCD.
# Four panels: A (DBD counts), B (DCD counts), C (DBD timing), D (DCD timing)

dir.create(here("figures", "supplement"), showWarnings = FALSE, recursive = TRUE)

ref_label <- "Median"

# Helper: for one donor type, return list of plot data (counts + timing) and top-10 pair labels
stratified_pair_data <- function(data, dtype) {
  sub <- data %>%
    filter(donation_type == dtype,
           seq_flag %in% c("<50 sequence", "\u226550 sequence", "AOOS (out-of-sequence)")) %>%
    mutate(
      OPO_CTR_CODE_current = case_when(
        OPO_CTR_CODE == "00527" ~ "11036", OPO_CTR_CODE == "05239" ~ "14353",
        OPO_CTR_CODE == "00093" ~ "19282", TRUE ~ OPO_CTR_CODE
      ),
      offer_type = ifelse(is_aoos, "AOOS", "Sequential")
    ) %>%
    filter(!is.na(OPO_CTR_CODE_current), !is.na(LISTING_CTR_CODE))

  pair_counts <- sub %>%
    group_by(OPO_CTR_CODE = OPO_CTR_CODE_current, LISTING_CTR_CODE, offer_type) %>%
    summarise(n_offers = n(), .groups = "drop") %>%
    pivot_wider(names_from = offer_type, values_from = n_offers, values_fill = 0) %>%
    mutate(
      OPO_anon = anonymize_opo_vec(OPO_CTR_CODE),
      Center_anon = anonymize_center_vec(LISTING_CTR_CODE),
      pair_label = paste0(OPO_anon, " \u2192 ", Center_anon),
      total = AOOS + Sequential,
      aoos_rate = if_else(total > 0, 100 * AOOS / total, 0)
    )

  one_per_center <- pair_counts %>%
    filter(AOOS > 0) %>%
    group_by(LISTING_CTR_CODE) %>%
    slice_max(AOOS, n = 1, with_ties = FALSE) %>%
    ungroup()

  med_counts <- pair_counts %>%
    inner_join(one_per_center %>% select(OPO_CTR_CODE, LISTING_CTR_CODE),
               by = c("OPO_CTR_CODE", "LISTING_CTR_CODE")) %>%
    summarise(AOOS = median(AOOS, na.rm = TRUE), Sequential = median(Sequential, na.rm = TRUE), .groups = "drop")

  above_med <- pair_counts %>%
    filter(AOOS > med_counts$AOOS[1], Sequential > med_counts$Sequential[1]) %>%
    arrange(desc(aoos_rate))

  top10 <- above_med %>% slice_head(n = 10) %>% pull(pair_label)
  median_pair_rate <- median(one_per_center$aoos_rate, na.rm = TRUE)
  median_pair_counts <- one_per_center %>%
    summarise(median_aoos = median(AOOS, na.rm = TRUE), median_sequential = median(Sequential, na.rm = TRUE), .groups = "drop")

  # Counts plot data (top 10 + median reference)
  pair_counts_plot <- above_med %>%
    filter(pair_label %in% top10) %>%
    select(pair_label, AOOS, Sequential) %>%
    pivot_longer(c(AOOS, Sequential), names_to = "offer_type", values_to = "n_offers") %>%
    mutate(pair_label_fct = factor(pair_label, levels = rev(top10))) %>%
    bind_rows(
      tibble(pair_label = "MED", pair_label_fct = factor(ref_label, levels = c(ref_label, rev(top10))),
             offer_type = "AOOS", n_offers = median_pair_counts$median_aoos),
      tibble(pair_label = "MED", pair_label_fct = factor(ref_label, levels = c(ref_label, rev(top10))),
             offer_type = "Sequential", n_offers = median_pair_counts$median_sequential)
    ) %>%
    mutate(pair_label_fct = factor(pair_label_fct, levels = c(ref_label, rev(top10))))

  pair_rates_labels <- above_med %>%
    filter(pair_label %in% top10) %>%
    select(pair_label, aoos_rate) %>%
    mutate(pair_label_fct = factor(pair_label, levels = rev(top10)),
           rate_label = paste0(round(aoos_rate, 1), "%")) %>%
    bind_rows(tibble(pair_label = "MED", pair_label_fct = factor(ref_label, levels = c(ref_label, rev(top10))),
                    aoos_rate = median_pair_rate, rate_label = paste0(round(median_pair_rate, 0), "%"))) %>%
    mutate(pair_label_fct = factor(pair_label_fct, levels = c(ref_label, rev(top10)))) %>%
    left_join(pair_counts_plot %>% filter(offer_type == "AOOS") %>% select(pair_label_fct, n_offers) %>% rename(label_y = n_offers),
              by = "pair_label_fct")

  # Timing: median time to clamp for top 10 pairs (from sub = this donor type only)
  pair_timing <- sub %>%
    mutate(
      offer_type = ifelse(is_aoos, "AOOS", "Sequential"),
      OPO_anon = anonymize_opo_vec(OPO_CTR_CODE_current),
      Center_anon = anonymize_center_vec(LISTING_CTR_CODE),
      pair_label = paste0(OPO_anon, " \u2192 ", Center_anon)
    ) %>%
    filter(!is.na(diff_final_clamp), !is.na(FINAL_RESPONSE_DT_tx),
           as.Date(FINAL_RESPONSE_DT_tx) <= TX_DATE,
           pair_label %in% top10) %>%
    group_by(pair_label, offer_type) %>%
    summarise(median_time_to_clamp = median(diff_final_clamp, na.rm = TRUE), .groups = "drop") %>%
    mutate(pair_label_fct = factor(pair_label, levels = c(ref_label, rev(top10))))

  med_timing <- sub %>%
    mutate(
      offer_type = ifelse(is_aoos, "AOOS", "Sequential"),
      OPO_anon = anonymize_opo_vec(OPO_CTR_CODE_current),
      Center_anon = anonymize_center_vec(LISTING_CTR_CODE),
      pair_label = paste0(OPO_anon, " \u2192 ", Center_anon)
    ) %>%
    filter(!is.na(diff_final_clamp)) %>%
    inner_join(one_per_center %>% select(OPO_CTR_CODE, LISTING_CTR_CODE),
               by = c("OPO_CTR_CODE_current" = "OPO_CTR_CODE", "LISTING_CTR_CODE")) %>%
    group_by(pair_label, offer_type) %>%
    summarise(median_time_to_clamp = median(diff_final_clamp, na.rm = TRUE), .groups = "drop") %>%
    group_by(offer_type) %>%
    summarise(median_time_to_clamp = median(median_time_to_clamp, na.rm = TRUE), .groups = "drop")

  pair_timing <- pair_timing %>%
    bind_rows(
      tibble(pair_label = "MED", pair_label_fct = factor(ref_label, levels = c(ref_label, rev(top10))),
             offer_type = "AOOS", median_time_to_clamp = med_timing$median_time_to_clamp[med_timing$offer_type == "AOOS"]),
      tibble(pair_label = "MED", pair_label_fct = factor(ref_label, levels = c(ref_label, rev(top10))),
             offer_type = "Sequential", median_time_to_clamp = med_timing$median_time_to_clamp[med_timing$offer_type == "Sequential"])
    ) %>%
    mutate(pair_label_fct = factor(pair_label_fct, levels = c(ref_label, rev(top10))))

  timing_wide <- pair_timing %>%
    pivot_wider(names_from = offer_type, values_from = median_time_to_clamp, id_cols = pair_label_fct)

  list(
    counts_plot = pair_counts_plot,
    rates_labels = pair_rates_labels,
    timing_plot = pair_timing,
    timing_wide = timing_wide,
    dtype = dtype
  )
}

dbd_dat <- stratified_pair_data(timing_cohort, "DBD")
dcd_dat <- stratified_pair_data(timing_cohort, "DCD")

# Summary stats for results (bottom line of Panels C and D)
med_dbd <- dbd_dat$timing_wide %>% filter(pair_label_fct == ref_label)
med_dcd <- dcd_dat$timing_wide %>% filter(pair_label_fct == ref_label)
if (nrow(med_dbd) > 0) {
  aoos_dbd <- med_dbd$AOOS[1]
  seq_dbd  <- med_dbd$Sequential[1]
  diff_dbd <- aoos_dbd - seq_dbd
  cat("DBD:  Median hours to clamp \u2014 AOOS: ", round(aoos_dbd, 2), " h, In-sequence: ", round(seq_dbd, 2), " h; difference (AOOS \u2212 In-sequence): ", round(diff_dbd, 2), " h\n", sep = "")
}
if (nrow(med_dcd) > 0) {
  aoos_dcd <- med_dcd$AOOS[1]
  seq_dcd  <- med_dcd$Sequential[1]
  diff_dcd <- aoos_dcd - seq_dcd
  cat("DCD:  Median hours to clamp \u2014 AOOS: ", round(aoos_dcd, 2), " h, In-sequence: ", round(seq_dcd, 2), " h; difference (AOOS \u2212 In-sequence): ", round(diff_dcd, 2), " h\n", sep = "")
}
cat("(Reference = median across one-pair-per-center pairs for each donor type.)\n\n")

# Top 10 only: median counts and median timing (for results)
top10_dbd_counts <- dbd_dat$counts_plot %>% filter(pair_label_fct != ref_label) %>% group_by(offer_type) %>% summarise(med = median(n_offers), .groups = "drop")
top10_dbd_timing <- dbd_dat$timing_wide %>% filter(pair_label_fct != ref_label)
top10_dcd_counts <- dcd_dat$counts_plot %>% filter(pair_label_fct != ref_label) %>% group_by(offer_type) %>% summarise(med = median(n_offers), .groups = "drop")
top10_dcd_timing <- dcd_dat$timing_wide %>% filter(pair_label_fct != ref_label)

# DBD
if (nrow(top10_dbd_counts) > 0) {
  med_aoos_dbd   <- top10_dbd_counts %>% filter(offer_type == "AOOS") %>% pull(med)
  med_seq_dbd    <- top10_dbd_counts %>% filter(offer_type == "Sequential") %>% pull(med)
  med_aoos_dbd   <- if (length(med_aoos_dbd) == 0) NA_real_ else med_aoos_dbd[1]
  med_seq_dbd    <- if (length(med_seq_dbd) == 0) NA_real_ else med_seq_dbd[1]
  cat("DBD \u2014 Counts (median of top 10): AOOS n = ", round(med_aoos_dbd, 0), ", In-sequence n = ", round(med_seq_dbd, 0), "\n", sep = "")
}
if (nrow(top10_dbd_timing) > 0) {
  t_aoos <- median(top10_dbd_timing$AOOS, na.rm = TRUE)
  t_seq  <- median(top10_dbd_timing$Sequential, na.rm = TRUE)
  cat("DBD \u2014 Proximity (median of top 10): AOOS ", round(t_aoos, 2), " h, In-sequence ", round(t_seq, 2), " h; difference (AOOS \u2212 In-sequence): ", round(t_aoos - t_seq, 2), " h\n", sep = "")
}
# DCD
if (nrow(top10_dcd_counts) > 0) {
  med_aoos_dcd   <- top10_dcd_counts %>% filter(offer_type == "AOOS") %>% pull(med)
  med_seq_dcd    <- top10_dcd_counts %>% filter(offer_type == "Sequential") %>% pull(med)
  med_aoos_dcd   <- if (length(med_aoos_dcd) == 0) NA_real_ else med_aoos_dcd[1]
  med_seq_dcd    <- if (length(med_seq_dcd) == 0) NA_real_ else med_seq_dcd[1]
  cat("DCD \u2014 Counts (median of top 10): AOOS n = ", round(med_aoos_dcd, 0), ", In-sequence n = ", round(med_seq_dcd, 0), "\n", sep = "")
}
if (nrow(top10_dcd_timing) > 0) {
  t_aoos <- median(top10_dcd_timing$AOOS, na.rm = TRUE)
  t_seq  <- median(top10_dcd_timing$Sequential, na.rm = TRUE)
  cat("DCD \u2014 Proximity (median of top 10): AOOS ", round(t_aoos, 2), " h, In-sequence ", round(t_seq, 2), " h; difference (AOOS \u2212 In-sequence): ", round(t_aoos - t_seq, 2), " h\n", sep = "")
}
cat("\n")

# Panel A: DBD counts
p_dbd_counts <- ggplot(dbd_dat$counts_plot, aes(x = pair_label_fct, y = n_offers, fill = offer_type)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(data = dbd_dat$rates_labels,
            aes(x = pair_label_fct, y = label_y, label = rate_label),
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 2.8, fontface = "bold", color = "white") +
  scale_fill_manual(values = c("AOOS" = GLOBAL_AOOS_COLOR, "Sequential" = GLOBAL_SEQUENTIAL_COLOR),
                    labels = c("AOOS" = "AOOS", "Sequential" = "In-Sequence"), name = "Allocation Type") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(tag = "A \u2013 DBD", x = "OPO \u2192 Listing Center", y = "Number of Transplants (n)") +
  theme_main(base_size = 11) +
  theme(plot.tag = element_text(size = 11, face = "bold", hjust = 0, vjust = 1), plot.tag.position = c(0, 1),
        axis.text.y = element_text(size = 8), legend.position = "right", panel.grid.major.y = element_blank())

# Panel B: DCD counts
p_dcd_counts <- ggplot(dcd_dat$counts_plot, aes(x = pair_label_fct, y = n_offers, fill = offer_type)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(data = dcd_dat$rates_labels,
            aes(x = pair_label_fct, y = label_y, label = rate_label),
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 2.8, fontface = "bold", color = "white") +
  scale_fill_manual(values = c("AOOS" = GLOBAL_AOOS_COLOR, "Sequential" = GLOBAL_SEQUENTIAL_COLOR),
                    labels = c("AOOS" = "AOOS", "Sequential" = "In-Sequence"), name = "Allocation Type") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(tag = "B \u2013 DCD", x = "OPO \u2192 Listing Center", y = "Number of Transplants (n)") +
  theme_main(base_size = 11) +
  theme(plot.tag = element_text(size = 11, face = "bold", hjust = 0, vjust = 1), plot.tag.position = c(0, 1),
        axis.text.y = element_text(size = 8), legend.position = "right", panel.grid.major.y = element_blank())

# Panel C: DBD timing
p_dbd_timing <- ggplot(dbd_dat$timing_plot, aes(x = pair_label_fct, y = median_time_to_clamp, fill = offer_type, color = offer_type)) +
  geom_segment(data = dbd_dat$timing_wide,
               aes(x = pair_label_fct, xend = pair_label_fct, y = AOOS, yend = Sequential),
               color = "gray70", linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(shape = 21, size = 3.5, stroke = 1.5, position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  scale_fill_manual(values = c("AOOS" = GLOBAL_AOOS_COLOR, "Sequential" = GLOBAL_SEQUENTIAL_COLOR),
                    labels = c("AOOS" = "AOOS", "Sequential" = "In-Sequence"), name = "Allocation Type") +
  scale_color_manual(values = c("AOOS" = GLOBAL_AOOS_COLOR, "Sequential" = GLOBAL_SEQUENTIAL_COLOR), guide = "none") +
  coord_flip() +
  labs(tag = "C \u2013 DBD", x = "OPO \u2192 Listing Center",
       y = "Hours from Clamp\n(negative = pre-clamp, 0 = at clamp, positive = post-clamp)") +
  theme_main(base_size = 11) +
  theme(plot.tag = element_text(size = 11, face = "bold", hjust = 0, vjust = 1), plot.tag.position = c(0, 1),
        axis.text.y = element_text(size = 8), legend.position = "right", panel.grid.major.y = element_blank())

# Panel D: DCD timing
p_dcd_timing <- ggplot(dcd_dat$timing_plot, aes(x = pair_label_fct, y = median_time_to_clamp, fill = offer_type, color = offer_type)) +
  geom_segment(data = dcd_dat$timing_wide,
               aes(x = pair_label_fct, xend = pair_label_fct, y = AOOS, yend = Sequential),
               color = "gray70", linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(shape = 21, size = 3.5, stroke = 1.5, position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  scale_fill_manual(values = c("AOOS" = GLOBAL_AOOS_COLOR, "Sequential" = GLOBAL_SEQUENTIAL_COLOR),
                    labels = c("AOOS" = "AOOS", "Sequential" = "In-Sequence"), name = "Allocation Type") +
  scale_color_manual(values = c("AOOS" = GLOBAL_AOOS_COLOR, "Sequential" = GLOBAL_SEQUENTIAL_COLOR), guide = "none") +
  coord_flip() +
  labs(tag = "D \u2013 DCD", x = "OPO \u2192 Listing Center",
       y = "Hours from Clamp\n(negative = pre-clamp, 0 = at clamp, positive = post-clamp)") +
  theme_main(base_size = 11) +
  theme(plot.tag = element_text(size = 11, face = "bold", hjust = 0, vjust = 1), plot.tag.position = c(0, 1),
        axis.text.y = element_text(size = 8), legend.position = "right", panel.grid.major.y = element_blank())

title_counts <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Number of Transplants\n(Top 10 OPO\u2013Center pairs by AOOS rate + median)",
           size = 4, fontface = "bold", hjust = 0.5, vjust = 0.5) +
  xlim(0, 1) + ylim(0, 1) + theme_void() + theme(plot.margin = margin(5, 0, 5, 0))

title_timing <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Proximity of Acceptance to Cross-Clamp\n(Median hours from acceptance to clamp; 0 = at clamp)",
           size = 4, fontface = "bold", hjust = 0.5, vjust = 0.5) +
  xlim(0, 1) + ylim(0, 1) + theme_void() + theme(plot.margin = margin(5, 0, 5, 0))

# Layout: title, (A, B), title, (C, D)
sf7_opo_center <- title_counts /
  (p_dbd_counts + p_dcd_counts) /
  title_timing /
  (p_dbd_timing + p_dcd_timing) +
  plot_layout(nrow = 4, heights = c(0.12, 1, 0.12, 1), guides = "collect") &
  theme(plot.margin = margin(5, 5, 5, 5))

print(sf7_opo_center)
ggsave(fig_path("figures", "supplement", "S7_OPO_Center_Stratified.png"), sf7_opo_center, width = 14, height = 12, dpi = FIG_DPI, compression = "lzw")
