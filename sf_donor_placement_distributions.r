# sf_donor_placement_distributions.r
# Donor and placement characteristics by allocation pathway.
# Panels: A DRI, B donor age, C CIT, D distance

library(tidyverse)
library(patchwork)
if (!requireNamespace("ggridges", quietly = TRUE)) {
  stop("Install the 'ggridges' package to generate this figure.")
}
library(ggridges)

ridge_theme <- theme_main(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0, margin = margin(b = 10)),
    plot.margin = margin(t = 10, r = 10, b = 4, l = 6),
    axis.text.y = element_text(size = 10)
  )

make_ridge <- function(df, col, title, xlab, xlim, p_lab, p_x,
                       x_breaks = NULL) {
  d <- df %>%
    filter(!is.na(.data[[col]]), !is.na(seq_flag)) %>%
    mutate(seq_flag = factor(seq_flag, levels = rev(levels(factor(seq_flag)))))

  # No left padding so ridges meet the y-axis.
  x_scale <- if (is.null(x_breaks)) {
    scale_x_continuous(expand = expansion(mult = c(0, 0.02)))
  } else {
    scale_x_continuous(breaks = x_breaks, expand = expansion(mult = c(0, 0.02)))
  }

  ggplot(d, aes(x = .data[[col]], y = seq_flag, fill = seq_flag, color = seq_flag)) +
    geom_density_ridges(
      alpha = 0.4, linewidth = 0.8, scale = 0.85,
      rel_min_height = 0.005
    ) +
    x_scale +
    scale_fill_manual(values = SEQ_COLORS_4) +
    scale_color_manual(values = SEQ_COLORS_4) +
    # Light top headroom so <50 peaks + p-labels are not clipped
    scale_y_discrete(expand = expansion(mult = c(0.05, 0.42))) +
    coord_cartesian(xlim = xlim) +
    labs(title = title, x = xlab, y = NULL) +
    annotate("text", x = p_x, y = 4.72, label = p_lab, hjust = 1, size = 3.5) +
    ridge_theme
}

fmt_p <- function(p) if (p < 0.001) "p < 0.001" else paste0("p = ", format(round(p, 3), nsmall = 3))

dat <- timing_cohort %>% filter(!is.na(seq_flag))

kw_dri <- kruskal.test(DRI ~ seq_flag, data = dat %>% filter(!is.na(DRI)))
kw_age <- kruskal.test(AGE_DON_donor ~ seq_flag, data = dat %>% filter(!is.na(AGE_DON_donor)))
kw_cit <- kruskal.test(COLD_ISCH ~ seq_flag, data = dat %>% filter(!is.na(COLD_ISCH)))
kw_dist <- kruskal.test(DISTANCE ~ seq_flag, data = dat %>% filter(!is.na(DISTANCE)))

p_dri <- make_ridge(
  dat, "DRI", "A. Donor Risk Index (DRI)", "DRI", c(0.8, 5),
  fmt_p(kw_dri$p.value), 4.8, x_breaks = 1:5
)

p_age <- make_ridge(
  dat, "AGE_DON_donor", "B. Donor age", "Donor age (years)", c(0, 85),
  fmt_p(kw_age$p.value), 78
)

p_cit <- make_ridge(
  dat, "COLD_ISCH", "C. Cold ischemia time", "CIT (hours)", c(0, 24),
  fmt_p(kw_cit$p.value), 22
)

p_dist <- make_ridge(
  dat, "DISTANCE", "D. Distance", "Distance (nautical miles)", c(0, 1500),
  fmt_p(kw_dist$p.value), 1400
)

fig_risk <- (p_dri | p_age) / (p_cit | p_dist)

print(fig_risk)
ggsave(
  fig_path("figures", "supplement", "S4_Donor_Placement_Distributions.png"),
  fig_risk, width = 14, height = 12, dpi = FIG_DPI, compression = "lzw"
)

message("DRI: ", fmt_p(kw_dri$p.value))
message("Donor age: ", fmt_p(kw_age$p.value))
message("CIT: ", fmt_p(kw_cit$p.value))
message("Distance: ", fmt_p(kw_dist$p.value))
