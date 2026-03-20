# sf_dri_dsri.r — SF-DRI/DSRI: Distributions by 4 groups (draft6 style: side-by-side, filled)

dri_dsri_complete <- timing_cohort %>%
  filter(!is.na(DRI), !is.na(DSRI), !is.na(seq_flag))

n_dri_dsri_by_group <- dri_dsri_complete %>%
  count(seq_flag, name = "n_complete")

kw_dri  <- kruskal.test(DRI  ~ seq_flag, data = dri_dsri_complete)
kw_dsri <- kruskal.test(DSRI ~ seq_flag, data = dri_dsri_complete)
fmt_p <- function(p) if (p < 0.001) "p < 0.001" else paste0("p = ", format(round(p, 3), nsmall = 3))
p_lab_dri  <- fmt_p(kw_dri$p.value)
p_lab_dsri <- fmt_p(kw_dsri$p.value)

# Ridgeline plots — each group in its own band, shared x-axis, no overlap
if (!requireNamespace("ggridges", quietly = TRUE)) install.packages("ggridges")
library(ggridges)

# Reverse factor order so <50 is at top (reads naturally top-to-bottom)
dri_dsri_complete <- dri_dsri_complete %>%
  mutate(seq_flag = factor(seq_flag, levels = rev(levels(seq_flag))))

# Shared ridge theme: more vertical room, title not encroaching, less dead space at bottom
ridge_theme <- theme_main(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0, margin = margin(b = 10)),
        plot.margin = margin(t = 8, r = 10, b = 4, l = 6),
        axis.text.y = element_text(size = 10))

# Panel A: DRI (x-axis 0–5, ticks at each integer)
p_dri <- ggplot(dri_dsri_complete, aes(x = DRI, y = seq_flag, fill = seq_flag, color = seq_flag)) +
  geom_density_ridges(alpha = 0.4, linewidth = 0.8, scale = 1.0, rel_min_height = 0.005) +
  scale_x_continuous(limits = c(0, 5), breaks = 1:5) +
  scale_fill_manual(values = SEQ_COLORS_4) +
  scale_color_manual(values = SEQ_COLORS_4) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.25))) +
  labs(title = "A. Donor Risk Index (DRI)", x = "DRI", y = NULL) +
  annotate("text", x = 4.8, y = 4.8, label = p_lab_dri, hjust = 1, size = 3.5) +
  ridge_theme

# Panel B: DSRI (x limited to 25; scale_x so density is computed on clipped data — smoother curve)
p_dsri <- ggplot(dri_dsri_complete, aes(x = DSRI, y = seq_flag, fill = seq_flag, color = seq_flag)) +
  geom_density_ridges(alpha = 0.4, linewidth = 0.8, scale = 1.0, rel_min_height = 0.005) +
  scale_x_continuous(limits = c(0, 25)) +
  scale_fill_manual(values = SEQ_COLORS_4) +
  scale_color_manual(values = SEQ_COLORS_4) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.25))) +
  labs(title = "B. Discard Risk Index (DSRI)", x = "DSRI", y = NULL) +
  annotate("text", x = 24, y = 4.8, label = p_lab_dsri, hjust = 1, size = 3.5) +
  ridge_theme

fig_dri_dsri <- p_dri + p_dsri

print(fig_dri_dsri)
ggsave(here("figures", "supplement", "SF_DRI_DSRI.png"), fig_dri_dsri, width = 14, height = 7, dpi = 300)

# Summary
cat("\n=== DRI/DSRI — N with complete data ===\n")
print(n_dri_dsri_by_group)
cat("Total:", sum(n_dri_dsri_by_group$n_complete), "\n")
cat("Kruskal-Wallis DRI:", p_lab_dri, "\n")
cat("Kruskal-Wallis DSRI:", p_lab_dsri, "\n")

cat("\n=== DRI Summary ===\n")
dri_dsri_complete %>%
  group_by(seq_flag) %>%
  summarise(n = n(), median = median(DRI, na.rm = TRUE),
            q25 = quantile(DRI, 0.25, na.rm = TRUE),
            q75 = quantile(DRI, 0.75, na.rm = TRUE), .groups = "drop") %>%
  print()

cat("\n=== DSRI Summary ===\n")
dri_dsri_complete %>%
  group_by(seq_flag) %>%
  summarise(n = n(), median = median(DSRI, na.rm = TRUE),
            q25 = quantile(DSRI, 0.25, na.rm = TRUE),
            q75 = quantile(DSRI, 0.75, na.rm = TRUE), .groups = "drop") %>%
  print()
