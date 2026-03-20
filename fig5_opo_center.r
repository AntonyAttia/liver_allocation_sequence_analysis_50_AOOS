# fig5_opo_center.r — Figure 5: OPO-Center Pair Analysis (AOOS only)
# Panel A: counts, Panel B: proximity to clamp

opo_center_counts <- timing_cohort %>%
  filter(seq_flag %in% c("<50 sequence", "\u226550 sequence", "AOOS (out-of-sequence)")) %>%
  mutate(
    OPO_CTR_CODE_current = case_when(
      OPO_CTR_CODE == "00527" ~ "11036", OPO_CTR_CODE == "05239" ~ "14353",
      OPO_CTR_CODE == "00093" ~ "19282", TRUE ~ OPO_CTR_CODE
    ),
    offer_type = ifelse(is_aoos, "AOOS", "Sequential")
  ) %>%
  filter(!is.na(OPO_CTR_CODE_current), !is.na(LISTING_CTR_CODE)) %>%
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

# One pair per center (highest AOOS)
one_pair_per_center <- opo_center_counts %>%
  filter(AOOS > 0) %>%
  group_by(LISTING_CTR_CODE) %>%
  slice_max(AOOS, n = 1, with_ties = FALSE) %>%
  ungroup()

median_center_counts <- opo_center_counts %>%
  inner_join(one_pair_per_center %>% select(OPO_CTR_CODE, LISTING_CTR_CODE),
             by = c("OPO_CTR_CODE", "LISTING_CTR_CODE")) %>%
  summarise(AOOS = median(AOOS, na.rm = TRUE), Sequential = median(Sequential, na.rm = TRUE))

median_aoos_count <- median_center_counts$AOOS[1]
median_seq_count  <- median_center_counts$Sequential[1]
median_pair_rate  <- median(one_pair_per_center$aoos_rate, na.rm = TRUE)

# Above-median pairs, top 10
above_median_pairs <- opo_center_counts %>%
  filter(AOOS > median_aoos_count, Sequential > median_seq_count) %>%
  arrange(desc(aoos_rate))

top_10_pairs <- above_median_pairs %>% slice_head(n = 10) %>% pull(pair_label)
ref_label <- "Median"

# Panel A data
pair_counts_plot <- above_median_pairs %>%
  filter(pair_label %in% top_10_pairs) %>%
  select(pair_label, AOOS, Sequential) %>%
  pivot_longer(c(AOOS, Sequential), names_to = "offer_type", values_to = "n_offers") %>%
  mutate(pair_label_fct = factor(pair_label, levels = rev(top_10_pairs))) %>%
  bind_rows(
    tibble(pair_label = "MED", pair_label_fct = factor(ref_label, levels = c(ref_label, rev(top_10_pairs))),
           offer_type = "AOOS", n_offers = median_center_counts$AOOS),
    tibble(pair_label = "MED", pair_label_fct = factor(ref_label, levels = c(ref_label, rev(top_10_pairs))),
           offer_type = "Sequential", n_offers = median_center_counts$Sequential)
  ) %>%
  mutate(pair_label_fct = factor(pair_label_fct, levels = c(ref_label, rev(top_10_pairs))))

pair_rates_labels <- above_median_pairs %>%
  filter(pair_label %in% top_10_pairs) %>%
  select(pair_label, aoos_rate) %>%
  mutate(pair_label_fct = factor(pair_label, levels = rev(top_10_pairs)),
         rate_label = paste0(round(aoos_rate, 1), "%")) %>%
  bind_rows(tibble(pair_label = "MED",
                   pair_label_fct = factor(ref_label, levels = c(ref_label, rev(top_10_pairs))),
                   aoos_rate = median_pair_rate,
                   rate_label = paste0(round(median_pair_rate, 1), "%"))) %>%
  mutate(pair_label_fct = factor(pair_label_fct, levels = c(ref_label, rev(top_10_pairs))))

pair_rates_labels <- pair_rates_labels %>%
  left_join(pair_counts_plot %>% filter(offer_type == "AOOS") %>%
              select(pair_label_fct, n_offers) %>% rename(label_y = n_offers),
            by = "pair_label_fct")

p1 <- ggplot(pair_counts_plot, aes(x = pair_label_fct, y = n_offers, fill = offer_type)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(data = pair_rates_labels,
            aes(x = pair_label_fct, y = label_y, label = rate_label),
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 3,
            fontface = "bold", color = "white") +
  scale_fill_manual(values = c("AOOS" = GLOBAL_AOOS_COLOR, "Sequential" = GLOBAL_SEQUENTIAL_COLOR),
                    labels = c("AOOS" = "AOOS", "Sequential" = "In-Sequence"),
                    name = "Allocation Type") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Panel A: Number of Transplants",
       subtitle = "Top 10 OPO-Center pairs by AOOS rate (above-median) + median",
       x = "OPO \u2192 Listing Center", y = "Number of Transplants (n)") +
  theme_gastro(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
        axis.text.y = element_text(size = 8), legend.position = "right",
        panel.grid.major.y = element_blank())

# Panel B: proximity to clamp
pair_timing_plot <- timing_cohort %>%
  filter(seq_flag %in% c("<50 sequence", "\u226550 sequence", "AOOS (out-of-sequence)")) %>%
  mutate(
    OPO_CTR_CODE_current = case_when(
      OPO_CTR_CODE == "00527" ~ "11036", OPO_CTR_CODE == "05239" ~ "14353",
      OPO_CTR_CODE == "00093" ~ "19282", TRUE ~ OPO_CTR_CODE
    ),
    offer_type = ifelse(is_aoos, "AOOS", "Sequential"),
    OPO_anon = anonymize_opo_vec(OPO_CTR_CODE_current),
    Center_anon = anonymize_center_vec(LISTING_CTR_CODE),
    pair_label = paste0(OPO_anon, " \u2192 ", Center_anon)
  ) %>%
  filter(!is.na(diff_final_clamp), !is.na(FINAL_RESPONSE_DT_tx),
         as.Date(FINAL_RESPONSE_DT_tx) <= TX_DATE,
         pair_label %in% top_10_pairs) %>%
  group_by(pair_label, offer_type) %>%
  summarise(median_time_to_clamp = median(diff_final_clamp, na.rm = TRUE), .groups = "drop") %>%
  mutate(pair_label_fct = factor(pair_label, levels = c(ref_label, rev(top_10_pairs))))

# Median reference timing
med_timing <- timing_cohort %>%
  filter(seq_flag %in% c("<50 sequence", "\u226550 sequence", "AOOS (out-of-sequence)")) %>%
  mutate(
    OPO_CTR_CODE_current = case_when(
      OPO_CTR_CODE == "00527" ~ "11036", OPO_CTR_CODE == "05239" ~ "14353",
      OPO_CTR_CODE == "00093" ~ "19282", TRUE ~ OPO_CTR_CODE
    ),
    offer_type = ifelse(is_aoos, "AOOS", "Sequential"),
    OPO_anon = anonymize_opo_vec(OPO_CTR_CODE_current),
    Center_anon = anonymize_center_vec(LISTING_CTR_CODE),
    pair_label = paste0(OPO_anon, " \u2192 ", Center_anon)
  ) %>%
  filter(!is.na(diff_final_clamp)) %>%
  inner_join(one_pair_per_center %>% select(OPO_CTR_CODE, LISTING_CTR_CODE),
             by = c("OPO_CTR_CODE_current" = "OPO_CTR_CODE", "LISTING_CTR_CODE")) %>%
  group_by(pair_label, offer_type) %>%
  summarise(median_time_to_clamp = median(diff_final_clamp, na.rm = TRUE), .groups = "drop") %>%
  group_by(offer_type) %>%
  summarise(median_time_to_clamp = median(median_time_to_clamp, na.rm = TRUE), .groups = "drop")

pair_timing_plot <- pair_timing_plot %>%
  bind_rows(
    tibble(pair_label = "MED",
           pair_label_fct = factor(ref_label, levels = c(ref_label, rev(top_10_pairs))),
           offer_type = "AOOS",
           median_time_to_clamp = med_timing$median_time_to_clamp[med_timing$offer_type == "AOOS"]),
    tibble(pair_label = "MED",
           pair_label_fct = factor(ref_label, levels = c(ref_label, rev(top_10_pairs))),
           offer_type = "Sequential",
           median_time_to_clamp = med_timing$median_time_to_clamp[med_timing$offer_type == "Sequential"])
  ) %>%
  mutate(pair_label_fct = factor(pair_label_fct, levels = c(ref_label, rev(top_10_pairs))))

panel2_wide <- pair_timing_plot %>%
  pivot_wider(names_from = offer_type, values_from = median_time_to_clamp, id_cols = pair_label_fct)

p2 <- ggplot(pair_timing_plot, aes(x = pair_label_fct, y = median_time_to_clamp,
                                    fill = offer_type, color = offer_type)) +
  geom_segment(data = panel2_wide,
               aes(x = pair_label_fct, xend = pair_label_fct, y = AOOS, yend = Sequential),
               color = "gray70", linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(shape = 21, size = 4, stroke = 1.5, position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  scale_fill_manual(values = c("AOOS" = GLOBAL_AOOS_COLOR, "Sequential" = GLOBAL_SEQUENTIAL_COLOR),
                    labels = c("AOOS" = "AOOS", "Sequential" = "In-Sequence"), name = "Allocation Type") +
  scale_color_manual(values = c("AOOS" = GLOBAL_AOOS_COLOR, "Sequential" = GLOBAL_SEQUENTIAL_COLOR),
                     guide = "none") +
  coord_flip() +
  labs(title = "Panel B: Proximity of Acceptance to Cross-Clamp",
       subtitle = "Median hours from acceptance to cross-clamp (0 = at clamp)",
       x = "OPO \u2192 Listing Center",
       y = "Hours from Clamp\n(negative = pre-clamp, 0 = at clamp, positive = post-clamp)") +
  theme_gastro(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
        axis.text.y = element_text(size = 8), legend.position = "right",
        panel.grid.major.y = element_blank())

fig5 <- p1 / p2 + plot_layout(nrow = 2, heights = c(1, 1))
print(fig5)
ggsave(here("figures", "main", "Figure5_OPO_Center_Pairs.png"), fig5, width = 14, height = 14, dpi = 300)
