# sf_machine_perfusion.r — SF-MP: Machine Perfusion Analysis (4 groups)

# MP proportions by 4-group
mp_summary <- timing_cohort %>%
  filter(!is.na(seq_flag)) %>%
  group_by(seq_flag) %>%
  summarise(
    n_total = n(),
    n_mp_optn_yes = sum(LI_MACHINE_PERFUSION_original == "Y", na.rm = TRUE),
    pct_mp_optn = round(100 * n_mp_optn_yes / n_total, 1),
    n_mp_extrap_yes = sum(LI_MACHINE_PERFUSION_FDA_era == "Yes", na.rm = TRUE),
    pct_mp_extrap = round(100 * n_mp_extrap_yes / n_total, 1),
    .groups = "drop"
  )

cat("\n=== Machine Perfusion by Allocation Pathway ===\n")
print(mp_summary)

# MP proportion bar chart (OPTN coded)
p_mp_bar <- mp_summary %>%
  ggplot(aes(x = seq_flag, y = pct_mp_optn, fill = seq_flag)) +
  geom_col(alpha = 0.9, color = "black", linewidth = 0.5, width = 0.7) +
  geom_text(aes(label = paste0(pct_mp_optn, "%\n(n=", n_mp_optn_yes, ")")),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = SEQ_COLORS_4, guide = "none") +
  ylim(0, max(mp_summary$pct_mp_optn) * 1.3) +
  labs(title = "A: Machine Perfusion Rate (OPTN Coded)",
       x = "Allocation Pathway", y = "% Machine Perfusion") +
  theme_main(base_size = 12) +
  theme(axis.text.x = element_text(size = 9, angle = 15, hjust = 1))

# MP proportion bar chart (extrapolated)
p_mp_extrap <- mp_summary %>%
  ggplot(aes(x = seq_flag, y = pct_mp_extrap, fill = seq_flag)) +
  geom_col(alpha = 0.9, color = "black", linewidth = 0.5, width = 0.7) +
  geom_text(aes(label = paste0(pct_mp_extrap, "%\n(n=", n_mp_extrap_yes, ")")),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = SEQ_COLORS_4, guide = "none") +
  ylim(0, max(mp_summary$pct_mp_extrap) * 1.3) +
  labs(title = "B: Machine Perfusion Rate (Extrapolated)",
       x = "Allocation Pathway", y = "% Machine Perfusion") +
  theme_main(base_size = 12) +
  theme(axis.text.x = element_text(size = 9, angle = 15, hjust = 1))

# One panel per pathway: each has 3 lines (Without MP, OPTN coded, Extrapolated)
# MP_type: mutually exclusive
mp_timing <- timing_cohort %>%
  filter(!is.na(LI_MACHINE_PERFUSION_original), !is.na(LI_MACHINE_PERFUSION_FDA_era), !is.na(diff_final_clamp)) %>%
  mutate(
    MP_type = case_when(
      LI_MACHINE_PERFUSION_original == "Y" ~ "With MP (OPTN coded)",
      LI_MACHINE_PERFUSION_FDA_era == "Yes" ~ "With MP (extrapolated)",
      TRUE ~ "Without MP"
    ),
    MP_type = factor(MP_type, levels = c("Without MP", "With MP (OPTN coded)", "With MP (extrapolated)"))
  )

# Colours/linetypes for the 3 MP categories (same in every panel)
mp_type_colors <- c("Without MP" = "gray40", "With MP (OPTN coded)" = "#0066CC", "With MP (extrapolated)" = "#CC6600")
mp_type_linetypes <- c("Without MP" = "dashed", "With MP (OPTN coded)" = "solid", "With MP (extrapolated)" = "solid")

# Build one timing plot per pathway
pathway_order <- c("<50 sequence", "\u226550 sequence", "AOOS (out-of-sequence)", "EP (Policy 9.10)")
pathway_labels_short <- c("<50 sequence" = "<50", "\u226550 sequence" = "\u226550",
                          "AOOS (out-of-sequence)" = "AOOS", "EP (Policy 9.10)" = "EP")

p_mp_timing_list <- list()
for (i in seq_along(pathway_order)) {
  pathway_i <- pathway_order[i]
  lab <- pathway_labels_short[pathway_i]
  dat_i <- mp_timing %>% filter(seq_flag == pathway_i)
  p_mp_timing_list[[i]] <- ggplot(dat_i, aes(x = diff_final_clamp, color = MP_type, linetype = MP_type)) +
    geom_density(linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_color_manual(values = mp_type_colors, name = NULL) +
    scale_linetype_manual(values = mp_type_linetypes, name = NULL) +
    coord_cartesian(xlim = c(-30, 30)) +
    labs(title = paste0(lab, ": Acceptance Timing by MP"),
         x = "Hours from Cross-Clamp", y = "Density") +
    theme_main(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 11))
}

# 6-panel figure: row 1 = A, B; row 2 = C,D,E,F (one timing panel per pathway)
row1 <- p_mp_bar | p_mp_extrap
row2 <- (p_mp_timing_list[[1]] | p_mp_timing_list[[2]]) / (p_mp_timing_list[[3]] | p_mp_timing_list[[4]])
fig_mp <- row1 / row2 +
  plot_layout(heights = c(1, 1.2), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

print(fig_mp)
ggsave(here("figures", "supplement", "SF_Machine_Perfusion.png"), fig_mp, width = 14, height = 14, dpi = 300)

# Reportable stats: median (IQR) hours from clamp, % after 0h and 6h (for caption/results)
mp_report <- mp_timing %>%
  mutate(With_MP = MP_type != "Without MP") %>%
  group_by(seq_flag, MP_type, With_MP) %>%
  summarise(
    n = n(),
    median_h = round(median(diff_final_clamp), 1),
    q25 = round(quantile(diff_final_clamp, 0.25), 1),
    q75 = round(quantile(diff_final_clamp, 0.75), 1),
    pct_after_0h = round(100 * mean(diff_final_clamp > 0, na.rm = TRUE), 1),
    pct_after_6h = round(100 * mean(diff_final_clamp > 6, na.rm = TRUE), 1),
    .groups = "drop"
  )

pathway_short <- c("<50 sequence" = "<50", "\u226550 sequence" = "\u226550",
                   "AOOS (out-of-sequence)" = "AOOS", "EP (Policy 9.10)" = "EP")

cat("\n=== Acceptance timing by MP (Panels C–F) — reportable ===\n\n")
for (path in pathway_order) {
  lab <- pathway_short[path]
  r <- mp_report %>% filter(seq_flag == path)
  cat(lab, ":\n")
  for (i in seq_len(nrow(r))) {
    ri <- r[i, ]
    cat("  ", as.character(ri$MP_type), " (n=", ri$n, "): median ", ri$median_h, " h (IQR ", ri$q25, "–", ri$q75, "); ",
        ri$pct_after_0h, "% after 0 h, ", ri$pct_after_6h, "% after 6 h\n", sep = "")
  }
  cat("\n")
}

# One-line summary for AOOS and EP (with MP allows later acceptance)
cat("=== Summary (Panels E & F) ===\n")
mp_timing_with <- mp_timing %>% mutate(With_MP = MP_type != "Without MP")
for (path in c("AOOS (out-of-sequence)", "EP (Policy 9.10)")) {
  lab <- pathway_short[path]
  d <- mp_timing_with %>% filter(seq_flag == path)
  without <- d %>% filter(!With_MP) %>% summarise(median_h = median(diff_final_clamp), pct_6h = 100 * mean(diff_final_clamp > 6), .groups = "drop")
  with_mp  <- d %>% filter(With_MP)  %>% summarise(median_h = median(diff_final_clamp), pct_6h = 100 * mean(diff_final_clamp > 6), .groups = "drop")
  if (nrow(without) == 1 && nrow(with_mp) == 1) {
    cat(lab, ": With MP, median ", round(with_mp$median_h, 1), " h from clamp (", round(with_mp$pct_6h, 1), "% after 6 h) vs without MP median ", round(without$median_h, 1), " h (", round(without$pct_6h, 1), "% after 6 h).\n", sep = "")
  }
}
cat("\n")

# ---- QUICK FIGURE: 3 groups (EP+AOOS combined), 4 panels stacked — REMOVE AFTER ----
timing_3 <- timing_cohort %>%
  filter(!is.na(seq_flag)) %>%
  mutate(
    seq_flag_3 = case_when(
      seq_flag %in% c("AOOS (out-of-sequence)", "EP (Policy 9.10)") ~ "AOOS",
      TRUE ~ as.character(seq_flag)
    ),
    seq_flag_3 = factor(seq_flag_3, levels = c("<50 sequence", "\u226550 sequence", "AOOS"))
  )
mp3 <- timing_3 %>%
  group_by(seq_flag_3) %>%
  summarise(
    n_total = n(),
    n_mp_extrap = sum(LI_MACHINE_PERFUSION_FDA_era == "Yes", na.rm = TRUE),
    pct_mp = 100 * n_mp_extrap / n_total,
    .groups = "drop"
  )
p_bar3 <- mp3 %>%
  ggplot(aes(x = seq_flag_3, y = pct_mp, fill = seq_flag_3)) +
  geom_col(alpha = 0.9, color = "black", linewidth = 0.4, width = 0.6) +
  geom_text(aes(label = paste0(round(pct_mp, 1), "%")), vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c(SEQ_COLORS_4[1], SEQ_COLORS_4[2], SEQ_COLORS_4[3]), guide = "none") +
  ylim(0, max(mp3$pct_mp) * 1.15) +
  labs(x = "Allocation pathway", y = "Machine perfusion rate (%)") +
  theme_main(base_size = 10) +
  theme(axis.text.x = element_text(size = 9, angle = 15, hjust = 1))

# Same density curves as original panels C–F: 2 lines = Without / With machine perfusion (extrapolated only)
timing_3_lines <- timing_3 %>%
  filter(!is.na(diff_final_clamp), !is.na(LI_MACHINE_PERFUSION_FDA_era)) %>%
  mutate(
    MP_status = if_else(LI_MACHINE_PERFUSION_FDA_era == "Yes", "With machine perfusion", "Without machine perfusion"),
    MP_status = factor(MP_status, levels = c("Without machine perfusion", "With machine perfusion"))
  )
pathways_3 <- c("<50 sequence", "\u226550 sequence", "AOOS")
pathway_titles_3 <- c("<50 sequence" = "<50 sequence", "\u226550 sequence" = "\u226550 sequence", "AOOS" = "AOOS")
mp_line_colors <- c("Without machine perfusion" = "gray40", "With machine perfusion" = "#CC6600")
p_line_list <- list()
for (i in seq_along(pathways_3)) {
  path_i <- pathways_3[i]
  title_i <- pathway_titles_3[path_i]
  dat_i <- timing_3_lines %>% filter(seq_flag_3 == path_i)
  p_line_list[[i]] <- ggplot(dat_i, aes(x = diff_final_clamp, color = MP_status, linetype = MP_status)) +
    geom_density(linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_color_manual(values = mp_line_colors, name = NULL) +
    scale_linetype_manual(values = c("Without machine perfusion" = "dashed", "With machine perfusion" = "solid"), name = NULL) +
    coord_cartesian(xlim = c(-30, 30)) +
    labs(title = title_i, x = "Hours from cross-clamp at acceptance", y = "Density") +
    theme_main(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10))
}
fig_quick <- p_bar3 / p_line_list[[1]] / p_line_list[[2]] / p_line_list[[3]] +
  plot_layout(heights = c(1, 1, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")
print(fig_quick)
ggsave(here("figures", "supplement", "quick_MP_3group_4panel.png"), fig_quick, width = 8, height = 12, dpi = 300)
