# fig2_refusal_timing.r: Figure 2. Refusal patterns (Panel A) + Placement timing (Panel B)
# 4-group: <50, ≥50, AOOS, EP. EP gets its own bin on placement timing x-axis

# Panel A: refusal bar chart (4 groups)
timing_cohort_keys <- timing_cohort %>% distinct(DONOR_ID, MATCH_ID_CODE)

donor_run_refusal_stats_filtered <- donor_run_refusal_stats %>%
  inner_join(timing_cohort_keys, by = c("DONOR_ID", "MATCH_ID_CODE"))

comparison_summary <- donor_run_refusal_stats_filtered %>%
  mutate(total_distinct_refusals = n_multi_refusal_decisions + n_distinct_refusals) %>%
  mutate(alloc_label = case_when(
    seq_flag %in% c("<50 sequence", "\u226550 sequence") ~ "In-Sequence",
    seq_flag == "AOOS (out-of-sequence)" ~ "AOOS",
    seq_flag == "EP (Policy 9.10)" ~ "EP",
    TRUE ~ "Other"
  )) %>%
  group_by(alloc_label) %>%
  summarise(
    median_multi_refusals       = median(n_multi_refusal_decisions),
    median_distinct_refusals    = median(n_distinct_refusals),
    median_total_refusals       = median(total_distinct_refusals),
    median_unique_centers       = median(n_unique_ctr_contacted, na.rm = TRUE),
    median_recipients_bundled   = median(n_recipients_bundled),
    .groups = "drop"
  )

plot_comparison <- comparison_summary %>%
  select(alloc_label, median_unique_centers, median_distinct_refusals,
         median_multi_refusals, median_recipients_bundled) %>%
  pivot_longer(-alloc_label, names_to = "metric", values_to = "value") %>%
  mutate(
    metric_clean = case_when(
      metric == "median_unique_centers"       ~ "Unique centers contacted",
      metric == "median_distinct_refusals"    ~ "Single patient refusals",
      metric == "median_multi_refusals"       ~ "Multi-patient refusal events",
      metric == "median_recipients_bundled"   ~ "Patients in multi-refusals"
    ),
    metric_clean = factor(metric_clean, levels = c(
      "Patients in multi-refusals", "Multi-patient refusal events",
      "Single patient refusals", "Unique centers contacted"
    )),
    alloc_label = factor(alloc_label, levels = c("In-Sequence", "AOOS", "EP"))
  )

bar_colors <- c("In-Sequence" = GLOBAL_SEQUENTIAL_COLOR,
                "AOOS" = GLOBAL_AOOS_COLOR,
                "EP" = GLOBAL_EP_COLOR)

p_refusals <- ggplot(plot_comparison, aes(x = metric_clean, y = value, fill = alloc_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.9) +
  geom_text(aes(label = ifelse(value == 0, "0",
                 ifelse(value < 1, sprintf("%.1f", value), sprintf("%.0f", value)))),
            position = position_dodge(width = 0.7),
            hjust = -0.1, size = 4, fontface = "bold", color = "black", family = "Helvetica") +
  scale_fill_manual(values = bar_colors, name = "Allocation Type") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "A", x = NULL, y = "Median count per donor-run") +
  theme_main(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0, margin = margin(b = 5)),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    axis.text.y = element_text(size = 11)
  )

# Panel B: placement timing (dual axis)
# Sequence bins: 1..5, 6-10, 11-20, 21-30, 31-50, 51-100, 101-500, 500+, AOOS, EP
timing_summary <- timing_cohort %>%
  mutate(
    sequence_bin = case_when(
      seq_flag == "AOOS (out-of-sequence)" ~ "AOOS",
      seq_flag == "EP (Policy 9.10)"       ~ "EP",
      PTR_SEQUENCE_NUM_tx == 1  ~ "1",
      PTR_SEQUENCE_NUM_tx == 2  ~ "2",
      PTR_SEQUENCE_NUM_tx == 3  ~ "3",
      PTR_SEQUENCE_NUM_tx == 4  ~ "4",
      PTR_SEQUENCE_NUM_tx == 5  ~ "5",
      PTR_SEQUENCE_NUM_tx >= 6   & PTR_SEQUENCE_NUM_tx <= 10  ~ "6-10",
      PTR_SEQUENCE_NUM_tx >= 11  & PTR_SEQUENCE_NUM_tx <= 20  ~ "11-20",
      PTR_SEQUENCE_NUM_tx >= 21  & PTR_SEQUENCE_NUM_tx <= 30  ~ "21-30",
      PTR_SEQUENCE_NUM_tx >= 31  & PTR_SEQUENCE_NUM_tx <= 50  ~ "31-50",
      PTR_SEQUENCE_NUM_tx >= 51  & PTR_SEQUENCE_NUM_tx <= 100 ~ "51-100",
      PTR_SEQUENCE_NUM_tx >= 101 & PTR_SEQUENCE_NUM_tx <= 500 ~ "101-500",
      PTR_SEQUENCE_NUM_tx > 500  ~ "500+",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(sequence_bin) %>%
  summarise(
    n_obs = n(),
    placement_time_median = median(placement_time_hours, na.rm = TRUE),
    clamp_timing_median   = median(diff_final_clamp, na.rm = TRUE),
    .groups = "drop"
  )

seq_order <- c("1","2","3","4","5","6-10","11-20","21-30",
               "31-50","51-100","101-500","500+","AOOS","EP")
timing_summary <- timing_summary %>%
  mutate(sequence_bin = factor(sequence_bin, levels = seq_order)) %>%
  arrange(sequence_bin)

y1_min <- min(timing_summary$placement_time_median, na.rm = TRUE)
y1_max <- max(timing_summary$placement_time_median, na.rm = TRUE)
y2_min <- min(timing_summary$clamp_timing_median, na.rm = TRUE)
y2_max <- max(timing_summary$clamp_timing_median, na.rm = TRUE)
scale_factor <- (y1_max - y1_min) / (y2_max - y2_min)
to_primary   <- function(y) (y - y2_min) * scale_factor + y1_min
to_secondary <- function(y) (y - y1_min) / scale_factor + y2_min

plot_data <- timing_summary %>%
  mutate(clamp_on_primary = to_primary(clamp_timing_median),
         is_oos = sequence_bin %in% c("AOOS", "EP"))

p_dual <- ggplot(plot_data, aes(x = sequence_bin, group = 1)) +
  geom_line(aes(y = placement_time_median, colour = "Placement time"), linewidth = 1.5) +
  geom_point(aes(y = placement_time_median, colour = "Placement time"), size = 3) +
  geom_line(aes(y = clamp_on_primary, colour = "Acceptance relative to cross-clamp"),
            linewidth = 1.5, linetype = "dashed") +
  geom_point(aes(y = clamp_on_primary, colour = "Acceptance relative to cross-clamp"),
             shape = 15, size = 3) +
  geom_point(data = subset(plot_data, sequence_bin == "AOOS"),
             aes(y = placement_time_median), shape = 8, size = 5, colour = GLOBAL_AOOS_COLOR) +
  geom_point(data = subset(plot_data, sequence_bin == "AOOS"),
             aes(y = clamp_on_primary), shape = 8, size = 5, colour = GLOBAL_AOOS_COLOR) +
  geom_point(data = subset(plot_data, sequence_bin == "EP"),
             aes(y = placement_time_median), shape = 8, size = 5, colour = GLOBAL_EP_COLOR) +
  geom_point(data = subset(plot_data, sequence_bin == "EP"),
             aes(y = clamp_on_primary), shape = 8, size = 5, colour = GLOBAL_EP_COLOR) +
  scale_y_continuous(
    name = "Placement time (hours, median)",
    limits = c(0, 32),
    breaks = seq(0, 32, 2),
    sec.axis = sec_axis(~ to_secondary(.), name = "Acceptance relative to cross-clamp (hours, median)",
                        breaks = seq(floor(y2_min), 4, 2))
  ) +
  scale_colour_manual(name = "Metric",
    values = c("Placement time" = "#003366",
               "Acceptance relative to cross-clamp" = "#8B4513")) +
  labs(x = "Sequence bin", title = "B") +
  theme_main(base_size = 12) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    plot.title = element_text(size = 16, face = "bold", hjust = 0, margin = margin(b = 5)),
    plot.margin = margin(t = 15, r = 10, b = 15, l = 10)
  )

# Combined
fig2_combined <- p_refusals | p_dual
print(fig2_combined)
ggsave(fig_path("figures", "main", "Figure2_Refusal_Timing.png"), fig2_combined, width = 22, height = 8, dpi = FIG_DPI, compression = "lzw")

# Print placement timing table
print(timing_summary, n = Inf)
