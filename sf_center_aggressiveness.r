# sf_center_aggressiveness.r: Supplement Figure 6. Center aggressiveness and distance differential
# Panels A and B: scatter plots, panels C and D: AOOS by tertile, panel E: distance difference by OPO

dir.create(here("figures", "supplement"), showWarnings = FALSE, recursive = TRUE)

# Sequential offers only (for aggressiveness)
sequential_offers_aggr <- timing_cohort %>%
  filter(!is_aoos, !is_ep, !is.na(LISTING_CTR_CODE), !is.na(n_unique_ctr_contacted)) %>%
  mutate(n_rejections_before_acceptance = n_unique_ctr_contacted - 1)

center_summary_aggr <- sequential_offers_aggr %>%
  group_by(LISTING_CTR_CODE) %>%
  summarise(
    n_offers = n(),
    max_rejections = max(n_rejections_before_acceptance, na.rm = TRUE),
    n_after_10_rejections = sum(n_rejections_before_acceptance >= 10),
    pct_after_10_rejections = round(100 * sum(n_rejections_before_acceptance >= 10) / n(), 1),
    .groups = "drop"
  )

centers_robust <- center_summary_aggr %>%
  arrange(n_after_10_rejections) %>%
  mutate(
    tertile = ntile(n_after_10_rejections, 3),
    category_by_count = case_when(
      tertile == 3 ~ "High: Frequently Accept After \u226510 Rejections",
      tertile == 2 ~ "Medium: Occasionally Accept After \u226510 Rejections",
      tertile == 1 ~ "Low: Rarely Accept After \u226510 Rejections"
    ),
    category_by_count = factor(category_by_count, levels = c(
      "Low: Rarely Accept After \u226510 Rejections",
      "Medium: Occasionally Accept After \u226510 Rejections",
      "High: Frequently Accept After \u226510 Rejections"
    ), ordered = TRUE)
  )

n_high   <- sum(centers_robust$tertile == 3)
n_medium <- sum(centers_robust$tertile == 2)
n_low    <- sum(centers_robust$tertile == 1)
legend_label_high   <- paste0("High: Frequently Accept (", n_high, " centers)")
legend_label_medium <- paste0("Medium: Occasionally Accept (", n_medium, " centers)")
legend_label_low    <- paste0("Low: Rarely Accept (", n_low, " centers)")

tertile_fill <- c(
  "High: Frequently Accept After \u226510 Rejections"   = "#D02670",
  "Medium: Occasionally Accept After \u226510 Rejections" = "#8A3FFC",
  "Low: Rarely Accept After \u226510 Rejections"      = "#24A148"
)

# Panel A: scatter (sequential after ≥10 vs total sequential)
p1_scatter <- centers_robust %>%
  ggplot(aes(x = n_after_10_rejections, y = n_offers, fill = category_by_count)) +
  geom_point(alpha = 0.8, stroke = 1, shape = 21, color = "black", size = 4) +
  scale_fill_manual(values = tertile_fill, guide = "none") +
  labs(title = "A", x = "No. In-Sequence Offers Accepted\nAfter \u226510 Centers Rejected", y = "Total No. In-Sequence Offers Accepted") +
  theme_main(base_size = 11) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0, margin = margin(b = 5)))

# Panel B: scatter (sequential after ≥10 vs AOOS per center)
aoos_by_center <- timing_cohort %>%
  filter(is_aoos, !is.na(LISTING_CTR_CODE)) %>%
  group_by(LISTING_CTR_CODE) %>%
  summarise(n_aoos_accepted = n(), .groups = "drop")

centers_with_aoos <- centers_robust %>%
  left_join(aoos_by_center, by = "LISTING_CTR_CODE") %>%
  mutate(n_aoos_accepted = coalesce(n_aoos_accepted, 0L))

# Spearman correlation (in-sequence after ≥10 rejections vs AOOS accepted) for Panel B
cor_test_result <- cor.test(centers_with_aoos$n_after_10_rejections,
                            centers_with_aoos$n_aoos_accepted,
                            method = "spearman", exact = FALSE)
spearman_rho <- round(cor_test_result$estimate, 3)
spearman_p <- cor_test_result$p.value
if (spearman_p < 0.001) {
  p_text <- "P < 0.001"
} else if (spearman_p < 0.01) {
  p_text <- "P < 0.01"
} else if (spearman_p < 0.05) {
  p_text <- "P < 0.05"
} else {
  p_text <- paste0("P = ", sprintf("%.3f", spearman_p))
}
cor_text <- paste0("\u03c1 = ", spearman_rho, ", ", p_text)
x_max <- max(centers_with_aoos$n_after_10_rejections, na.rm = TRUE)
y_max <- max(centers_with_aoos$n_aoos_accepted, na.rm = TRUE)
cat("Spearman \u03c1:", spearman_rho, "|", p_text, "\n")

p2_scatter <- centers_with_aoos %>%
  ggplot(aes(x = n_after_10_rejections, y = n_aoos_accepted, fill = category_by_count)) +
  geom_point(alpha = 0.8, stroke = 1, shape = 21, color = "black", size = 4) +
  scale_fill_manual(values = tertile_fill, guide = "none") +
  annotate("text",
           x = x_max * 0.02,
           y = y_max * 1.07,
           label = cor_text,
           hjust = 0,
           vjust = 1,
           size = 4,
           fontface = "bold",
           family = "Helvetica",
           color = "black") +
  labs(title = "B", x = "No. In-Sequence Offers Accepted\nAfter \u226510 Centers Rejected", y = "Total No. AOOS Offers Accepted") +
  theme_main(base_size = 11) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0, margin = margin(b = 5)))

# Panel C: single stacked bar (AOOS volume by tertile)
# Canonical order so stack and label positions match (Low = bottom, Medium = middle, High = top)
panel_c_levels <- c(
  "Low: Rarely Accept After \u226510 Rejections",
  "Medium: Occasionally Accept After \u226510 Rejections",
  "High: Frequently Accept After \u226510 Rejections"
)
aoos_by_count_aggr <- timing_cohort %>%
  filter(is_aoos, !is.na(LISTING_CTR_CODE)) %>%
  select(DONOR_ID, MATCH_ID_CODE, LISTING_CTR_CODE) %>%
  left_join(centers_robust %>% select(LISTING_CTR_CODE, category_by_count), by = "LISTING_CTR_CODE") %>%
  filter(!is.na(category_by_count)) %>%
  count(category_by_count, name = "n_aoos") %>%
  mutate(
    category_by_count = factor(category_by_count, levels = panel_c_levels, ordered = TRUE),
    pct_aoos = round(100 * n_aoos / sum(n_aoos), 1)
  ) %>%
  arrange(match(category_by_count, panel_c_levels))

# Label positions: use same order as ggplot stacks (with reverse=FALSE, last level = bottom)
# So stack is High (bottom), Medium (middle), Low (top). Build labels in that order.
panel_c_stack_order <- rev(panel_c_levels)  # High, Medium, Low
aoos_stack_labels <- aoos_by_count_aggr %>%
  arrange(match(category_by_count, panel_c_stack_order)) %>%
  mutate(
    ymin = c(0, cumsum(head(n_aoos, -1))),
    ymax = cumsum(n_aoos),
    ymid = (ymin + ymax) / 2,
    label_text = paste0(n_aoos, " (", pct_aoos, "%)"),
    label_color = unname(tertile_fill[as.character(category_by_count)])
  )

p3_stacked <- aoos_by_count_aggr %>%
  ggplot(aes(x = 1, y = n_aoos, fill = category_by_count)) +
  geom_col(position = position_stack(reverse = FALSE), alpha = 0.9, width = 0.5, color = "black", linewidth = 0.8) +
  geom_text(data = aoos_stack_labels, aes(x = 1.3, y = ymid, label = label_text, color = label_color),
            fontface = "bold", size = 4, hjust = 0, vjust = 0.5, inherit.aes = FALSE) +
  scale_color_identity() +
  scale_fill_manual(
    values = tertile_fill,
    labels = c("High: Frequently Accept After \u226510 Rejections" = legend_label_high,
               "Medium: Occasionally Accept After \u226510 Rejections" = legend_label_medium,
               "Low: Rarely Accept After \u226510 Rejections" = legend_label_low),
    name = "Center Category (Tertiles):\nAcceptance After \u226510 Prior Rejections"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_x_continuous(limits = c(0.5, 1.7), expand = c(0, 0)) +
  labs(title = "C", x = NULL, y = "AOOS transplant volume") +
  theme_main(base_size = 11) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0, margin = margin(b = 5)),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.text = element_text(size = 9))

# Panel D: median AOOS rate (%) by tertile
centers_with_transplants <- timing_cohort %>%
  filter(!is.na(LISTING_CTR_CODE)) %>%
  left_join(centers_robust %>% select(LISTING_CTR_CODE, category_by_count), by = "LISTING_CTR_CODE") %>%
  filter(!is.na(category_by_count)) %>%
  mutate(allocation_type = ifelse(is_aoos, "AOOS", "Sequential")) %>%
  group_by(LISTING_CTR_CODE, category_by_count, allocation_type) %>%
  summarise(n_transplants = n(), .groups = "drop") %>%
  pivot_wider(names_from = allocation_type, values_from = n_transplants, values_fill = 0) %>%
  mutate(total_transplants = AOOS + Sequential,
         aoos_rate = if_else(total_transplants > 0, 100 * AOOS / total_transplants, 0))

normalized_rates <- centers_with_transplants %>%
  group_by(category_by_count) %>%
  summarise(n_centers = n(), median_aoos_rate = median(aoos_rate, na.rm = TRUE), .groups = "drop") %>%
  mutate(category_by_count = factor(category_by_count, levels = c(
    "Low: Rarely Accept After \u226510 Rejections",
    "Medium: Occasionally Accept After \u226510 Rejections",
    "High: Frequently Accept After \u226510 Rejections"
  ), ordered = TRUE)) %>%
  arrange(category_by_count)

p4_median_rate <- normalized_rates %>%
  ggplot(aes(x = category_by_count, y = median_aoos_rate, fill = category_by_count)) +
  geom_col(alpha = 0.9, color = "black", linewidth = 0.8) +
  geom_text(aes(label = paste0(sprintf("%.1f", median_aoos_rate), "%"), y = median_aoos_rate / 2),
            vjust = 0.5, color = "white", size = 4, fontface = "bold") +
  scale_fill_manual(values = tertile_fill, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "D", x = NULL, y = "Median AOOS Rate (%)") +
  theme_main(base_size = 11) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0, margin = margin(b = 5)),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Panel E: distance difference by OPO (segmented by tertile)
center_acceptance_proportions <- timing_cohort %>%
  mutate(
    OPO_CTR_CODE_current = case_when(
      OPO_CTR_CODE == "00527" ~ "11036", OPO_CTR_CODE == "05239" ~ "14353",
      OPO_CTR_CODE == "00093" ~ "19282", TRUE ~ OPO_CTR_CODE
    )
  ) %>%
  filter(!is.na(DISTANCE), !is.na(LISTING_CTR_CODE)) %>%
  left_join(centers_robust %>% select(LISTING_CTR_CODE, category_by_count), by = "LISTING_CTR_CODE") %>%
  filter(!is.na(category_by_count)) %>%
  group_by(OPO_CTR_CODE = OPO_CTR_CODE_current, category_by_count) %>%
  summarise(n_offers = n(), .groups = "drop") %>%
  group_by(OPO_CTR_CODE) %>%
  mutate(total_offers = sum(n_offers), proportion = n_offers / total_offers) %>%
  ungroup()

distance_data <- timing_cohort %>%
  mutate(
    OPO_CTR_CODE_current = case_when(
      OPO_CTR_CODE == "00527" ~ "11036", OPO_CTR_CODE == "05239" ~ "14353",
      OPO_CTR_CODE == "00093" ~ "19282", TRUE ~ OPO_CTR_CODE
    ),
    offer_type = if_else(is_aoos, "AOOS", "Sequential")
  ) %>%
  filter(!is.na(DISTANCE)) %>%
  group_by(OPO_CTR_CODE = OPO_CTR_CODE_current, offer_type) %>%
  summarise(median_distance = median(DISTANCE, na.rm = TRUE), n = n(), .groups = "drop") %>%
  pivot_wider(names_from = offer_type, values_from = c(median_distance, n), names_sep = "_") %>%
  filter(!is.na(median_distance_AOOS), !is.na(median_distance_Sequential)) %>%
  mutate(
    distance_diff = median_distance_Sequential - median_distance_AOOS,
    OPO_short = anonymize_opo_vec(OPO_CTR_CODE),
    n_open_offers = n_AOOS
  ) %>%
  arrange(desc(distance_diff))

distance_plot_data <- distance_data %>% filter(abs(distance_diff) <= 300)
total_n_positive <- sum(distance_plot_data$n_open_offers[distance_plot_data$distance_diff > 0], na.rm = TRUE)
total_n_negative <- sum(distance_plot_data$n_open_offers[distance_plot_data$distance_diff < 0], na.rm = TRUE)

stacked_acceptance_data <- distance_plot_data %>%
  select(OPO_CTR_CODE, OPO_short, distance_diff, n_open_offers) %>%
  left_join(center_acceptance_proportions, by = "OPO_CTR_CODE") %>%
  arrange(OPO_CTR_CODE, match(category_by_count, c("High: Frequently Accept After \u226510 Rejections",
    "Medium: Occasionally Accept After \u226510 Rejections", "Low: Rarely Accept After \u226510 Rejections"))) %>%
  group_by(OPO_CTR_CODE) %>%
  mutate(
    category_by_count = factor(category_by_count, levels = c(
      "High: Frequently Accept After \u226510 Rejections",
      "Medium: Occasionally Accept After \u226510 Rejections",
      "Low: Rarely Accept After \u226510 Rejections"
    ), ordered = TRUE),
    cumsum_prop = cumsum(proportion),
    cumsum_prop_prev = lag(cumsum_prop, default = 0),
    ymin = if_else(distance_diff >= 0, distance_diff * cumsum_prop_prev, distance_diff * cumsum_prop),
    ymax = if_else(distance_diff >= 0, distance_diff * cumsum_prop, distance_diff * cumsum_prop_prev)
  ) %>%
  ungroup()

n_opos <- length(unique(distance_plot_data$OPO_short))
p5_distance <- stacked_acceptance_data %>%
  ggplot(aes(x = reorder(OPO_short, distance_diff))) +
  geom_rect(aes(ymin = ymin, ymax = ymax, fill = category_by_count), alpha = 0.9, width = 0.7, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = tertile_fill, guide = "none") +
  scale_y_continuous(limits = c(-300, 300), expand = expansion(mult = c(0.05, 0.1))) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  geom_text(data = distance_plot_data, aes(x = OPO_short, y = distance_diff, label = paste0("n=", n_open_offers)),
            hjust = ifelse(distance_plot_data$distance_diff >= 0, -0.1, 1.1), size = 2.8, fontface = "bold", inherit.aes = FALSE) +
  annotate("text", x = n_opos * 0.15, y = 80, label = paste0("AOOS more distant\nthan in-sequence\n(n=", total_n_negative, ")"),
           color = "#D02670", fontface = "bold", size = 4, hjust = 0, vjust = 0.5) +
  annotate("text", x = n_opos * 0.85, y = -200, label = paste0("AOOS more local\nthan in-sequence\n(n=", total_n_positive, ")"),
           color = "#24A148", fontface = "bold", size = 4, hjust = 0, vjust = 0.5) +
  coord_flip() +
  labs(title = "E", x = "OPO", y = "Distance Difference (nautical miles)\nIn-Sequence - AOOS") +
  theme_main(base_size = 11) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0, margin = margin(b = 5)),
        axis.text.y = element_text(size = 8))

# Layout: (A, B) / (C, D), then E, then legend below
p3_stacked_noleg <- p3_stacked + theme(legend.position = "none")
left_grid <- (p1_scatter | p2_scatter) / (p3_stacked_noleg | p4_median_rate) +
  plot_layout(heights = c(1, 1))

# Shared legend from Panel C
p3_legend <- p3_stacked +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.text = element_text(size = 9))
center_legend <- cowplot::get_legend(p3_legend)
legend_plot <- cowplot::ggdraw(center_legend) + theme(plot.margin = margin(8, 30, 8, 30))

main_figure <- left_grid | p5_distance + plot_layout(widths = c(1.6, 1))
sf6_center_aggr <- main_figure / legend_plot + plot_layout(heights = c(1, 0.08))

print(sf6_center_aggr)
ggsave(fig_path("figures", "supplement", "S6_Center_Aggressiveness.png"), sf6_center_aggr, width = 22, height = 12, dpi = FIG_DPI, compression = "lzw")

# Printed summary
tertile_aoos_summary <- centers_with_aoos %>%
  group_by(category_by_count) %>%
  summarise(n_centers = n(), total_aoos = sum(n_aoos_accepted), .groups = "drop")
cat("AOOS by tertile:\n")
print(tertile_aoos_summary)
cat("\nMedian AOOS rate by tertile (%):\n")
print(normalized_rates %>% select(category_by_count, n_centers, median_aoos_rate))
