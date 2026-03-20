# fig3_response_timing.r — Figure 3: Response timing relative to cross-clamp (4 groups)
# A/B: CDF by DBD/DCD, C: boxplot, D: % after clamp

final_plot_data <- timing_cohort %>%
  filter(!is.na(donation_type), !is.na(seq_flag))

sequence_colors <- SEQ_COLORS_4

box_medians <- final_plot_data %>%
  group_by(donation_type, seq_flag) %>%
  summarise(
    med = median(diff_final_clamp, na.rm = TRUE),
    q25 = quantile(diff_final_clamp, 0.25, na.rm = TRUE),
    q75 = quantile(diff_final_clamp, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
cat("\n=== Acceptance Timing Medians ===\n")
print(box_medians)

make_cumulative_plot <- function(data, donor_type_filter, plot_label) {
  pd <- data %>% filter(donation_type == donor_type_filter)
  max_val <- max(pd$diff_final_clamp, na.rm = TRUE)
  ggplot(pd, aes(x = diff_final_clamp, color = seq_flag)) +
    stat_ecdf(linewidth = 1.5, alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1, alpha = 0.7) +
    coord_cartesian(xlim = c(-50, max_val)) +
    scale_color_manual(values = sequence_colors, name = "Allocation Pathway") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), labels = percent_format(accuracy = 1)) +
    labs(title = plot_label, x = "Hours from Cross-Clamp", y = "Cumulative Proportion of Acceptance") +
    theme_main(base_size = 11) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0, margin = margin(b = 8)),
          legend.position = "none", plot.margin = margin(8, 8, 8, 8)) +
    annotate("text", x = -3, y = 0.5, label = "Cross-clamp",
             angle = 90, hjust = 0.5, vjust = 0.5, size = 3.5, fontface = "bold")
}

make_boxplot <- function(data, plot_label) {
  ggplot(data, aes(x = donation_type, y = diff_final_clamp, fill = seq_flag)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7, linewidth = 1) +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.9,
                 outlier.alpha = 0.5, outlier.size = 1.5, linewidth = 0.8, color = "black") +
    scale_fill_manual(values = sequence_colors, name = "Allocation Pathway") +
    coord_cartesian(ylim = c(-50, max(data$diff_final_clamp, na.rm = TRUE))) +
    labs(title = plot_label, x = "Donation Type",
         y = "Timing of Acceptance (hours from cross-clamp)") +
    theme_main(base_size = 11) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0, margin = margin(b = 8)),
          legend.position = "none", plot.margin = margin(8, 8, 8, 8)) +
    annotate("text", x = 0.4, y = 0, label = "Cross-clamp",
             hjust = 0, vjust = -0.3, size = 3.35, fontface = "bold")
}

make_proportion_plot <- function(data) {
  prop_data <- data %>%
    group_by(donation_type, seq_flag) %>%
    summarise(n_post = sum(diff_final_clamp > 0, na.rm = TRUE), total = n(), .groups = "drop") %>%
    mutate(prop = (n_post / total) * 100,
           label = paste0(n_post, "/", total, "\n(", round(prop, 1), "%)"))
  ggplot(prop_data, aes(x = donation_type, y = prop, fill = seq_flag)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.72, alpha = 0.9, color = "black", linewidth = 0.8) +
    geom_text(aes(label = label), position = position_dodge(width = 0.75),
              vjust = -0.3, hjust = 0.5, size = 2.7, fontface = "bold") +
    scale_fill_manual(values = sequence_colors, name = "Allocation Pathway") +
    labs(title = "D", x = "Donation Type", y = "% of Acceptances After Cross-Clamp") +
    ylim(0, max(prop_data$prop) * 1.2) +
    theme_main(base_size = 11) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0, margin = margin(b = 8)),
          legend.position = "none", plot.margin = margin(8, 8, 8, 8))
}

p_dbd <- make_cumulative_plot(final_plot_data, "DBD", "A: DBD")
p_dcd <- make_cumulative_plot(final_plot_data, "DCD", "B: DCD")
p_box <- make_boxplot(final_plot_data, "C")
p_prop <- make_proportion_plot(final_plot_data)

# Shared legend
legend_data <- data.frame(x = 1:4, y = 1, alloc = names(sequence_colors))
color_key <- ggplot(legend_data, aes(x = x, y = y, color = alloc)) +
  geom_point(size = 0) +
  scale_color_manual(values = sequence_colors, name = "Allocation Pathway") +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 11),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 5)))

fig3 <- ((p_dbd | p_dcd) / (p_box | p_prop)) / color_key +
  plot_layout(heights = c(1, 1, 0.12))

print(fig3)
ggsave(here("figures", "main", "Figure3_Response_Timing.png"), fig3, width = 14, height = 12, dpi = 300)

# Cumulative table
time_points <- c(-30, -20, -10, 0, 2, 4, 6, 8, 10, 20)
cum_list <- lapply(time_points, function(tp) {
  final_plot_data %>%
    group_by(donation_type, seq_flag) %>%
    summarise(time_point = tp, n_accepted = sum(diff_final_clamp <= tp, na.rm = TRUE),
              n_total = n(), cumulative_prop = (n_accepted / n_total) * 100, .groups = "drop")
})
cumulative_table <- bind_rows(cum_list) %>%
  mutate(seq_flag_clean = case_when(
    seq_flag == "<50 sequence" ~ "<50",
    seq_flag == "\u226550 sequence" ~ "\u226550",
    seq_flag == "AOOS (out-of-sequence)" ~ "AOOS",
    seq_flag == "EP (Policy 9.10)" ~ "EP"
  )) %>%
  select(donation_type, seq_flag_clean, time_point, cumulative_prop) %>%
  pivot_wider(names_from = c(donation_type, seq_flag_clean), values_from = cumulative_prop, names_sep = "_") %>%
  arrange(time_point) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.1f%%", .x)))

cat("\n=== Cumulative Acceptance Table ===\n")
print(knitr::kable(cumulative_table, format = "simple"))
