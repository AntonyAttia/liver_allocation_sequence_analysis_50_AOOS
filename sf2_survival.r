# sf2_survival.r — Figure S2: Kaplan-Meier survival + Supplemental Table S4 (4 groups)

surv_df <- timing_cohort %>%
  mutate(
    seq_flag = factor(seq_flag, levels = c("<50 sequence", "\u226550 sequence",
                                            "AOOS (out-of-sequence)", "EP (Policy 9.10)")),
    PTIME_yrs = PTIME / 365.25,
    GTIME_yrs = GTIME / 365.25,
    donation_type = factor(donation_type, levels = c("DBD", "DCD"))
  ) %>%
  filter(!is.na(donation_type), !is.na(PSTATUS), !is.na(GSTATUS)) %>%
  droplevels()

surv_dcd <- surv_df %>% filter(donation_type == "DCD") %>% droplevels()
surv_dbd <- surv_df %>% filter(donation_type == "DBD") %>% droplevels()

surv_dcd$seq_flag <- relevel(surv_dcd$seq_flag, ref = "<50 sequence")
surv_dbd$seq_flag <- relevel(surv_dbd$seq_flag, ref = "<50 sequence")

fit_patient_dcd <- survfit(Surv(PTIME_yrs, PSTATUS) ~ seq_flag, data = surv_dcd)
fit_graft_dcd   <- survfit(Surv(GTIME_yrs, GSTATUS) ~ seq_flag, data = surv_dcd)
fit_patient_dbd <- survfit(Surv(PTIME_yrs, PSTATUS) ~ seq_flag, data = surv_dbd)
fit_graft_dbd   <- survfit(Surv(GTIME_yrs, GSTATUS) ~ seq_flag, data = surv_dbd)

pval_patient_dcd <- surv_pvalue(fit_patient_dcd, data = surv_dcd)$pval
pval_graft_dcd   <- surv_pvalue(fit_graft_dcd, data = surv_dcd)$pval
pval_patient_dbd <- surv_pvalue(fit_patient_dbd, data = surv_dbd)$pval
pval_graft_dbd   <- surv_pvalue(fit_graft_dbd, data = surv_dbd)$pval

labs_short <- c("<50", "\u226550", "AOOS", "EP")
pal <- c("#0066CC", "#00468B", "#FF8C00", "#8A3FFC")

make_surv_panel <- function(fit, data, title_str) {
  ggsurvplot(
    fit, data = data, title = title_str,
    xlim = c(0, 2), ylim = c(0.6, 1), break.time.by = 0.5,
    xlab = "Years since transplant", ylab = "Survival probability",
    legend.title = "Pathway", legend.labs = labs_short, palette = pal,
    conf.int = FALSE, risk.table = TRUE, pval = TRUE,
    pval.coord = c(0.2, 0.65), pval.size = 5,
    risk.table.fontsize = 3, risk.table.y.text.col = TRUE,
    tables.theme = theme_void(),
    ggtheme = theme_gastro(base_size = 12)
  )
}

p1 <- make_surv_panel(fit_patient_dcd, surv_dcd, bquote(bold("A") ~ ". DCD - Patient Survival"))
p2 <- make_surv_panel(fit_graft_dcd, surv_dcd, bquote(bold("B") ~ ". DCD - Graft Survival"))
p3 <- make_surv_panel(fit_patient_dbd, surv_dbd, bquote(bold("C") ~ ". DBD - Patient Survival"))
p4 <- make_surv_panel(fit_graft_dbd, surv_dbd, bquote(bold("D") ~ ". DBD - Graft Survival"))

p1$plot <- p1$plot + theme(plot.title = element_text(size = 16, face = "plain", hjust = 0))
p2$plot <- p2$plot + theme(plot.title = element_text(size = 16, face = "plain", hjust = 0))
p3$plot <- p3$plot + theme(plot.title = element_text(size = 16, face = "plain", hjust = 0))
p4$plot <- p4$plot + theme(plot.title = element_text(size = 16, face = "plain", hjust = 0))

p1_stack <- p1$plot / p1$table + plot_layout(heights = c(4, 1))
p2_stack <- p2$plot / p2$table + plot_layout(heights = c(4, 1))
p3_stack <- p3$plot / p3$table + plot_layout(heights = c(4, 1))
p4_stack <- p4$plot / p4$table + plot_layout(heights = c(4, 1))

combined_survival <- (p1_stack | p2_stack) / (p3_stack | p4_stack) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        plot.margin = margin(6, 6, 6, 6))

print(combined_survival)
ggsave(here("figures", "supplement", "FigureSF2_Survival.png"), combined_survival, width = 14, height = 16, dpi = 300)

# Print p-values and 2-year estimates (follow plot xlim: 0–2 years since transplant)
cat("\n=== SURVIVAL P-VALUES (log-rank, 0–2 years since transplant) ===\n")
cat("DCD Patient:", sprintf("%.4f", pval_patient_dcd), "\n")
cat("DCD Graft:", sprintf("%.4f", pval_graft_dcd), "\n")
cat("DBD Patient:", sprintf("%.4f", pval_patient_dbd), "\n")
cat("DBD Graft:", sprintf("%.4f", pval_graft_dbd), "\n")

for (lbl in c("DBD graft", "DBD patient", "DCD graft", "DCD patient")) {
  fit_obj <- switch(lbl,
    "DBD graft"   = fit_graft_dbd,
    "DBD patient" = fit_patient_dbd,
    "DCD graft"   = fit_graft_dcd,
    "DCD patient" = fit_patient_dcd
  )
  s <- summary(fit_obj, times = 2)
  for (i in seq_along(s$surv)) {
    cat(sprintf("%s—%s: %.1f%% (%.1f–%.1f%%)\n",
        lbl, labs_short[i], s$surv[i]*100, s$lower[i]*100, s$upper[i]*100))
  }
}

# ---- SUPPLEMENTAL TABLE S4: Survival Outcomes ----
timing_cohort_surv <- timing_cohort %>%
  mutate(PTIME_yrs = PTIME / 365.25, GTIME_yrs = GTIME / 365.25) %>%
  filter(!is.na(donation_type), !is.na(PSTATUS), !is.na(GSTATUS)) %>%
  droplevels()

results_list <- list()
for (dtype in c("DBD", "DCD")) {
  sub_df <- timing_cohort_surv %>% filter(donation_type == dtype)
  for (outcome in c("Patient", "Graft")) {
    time_var  <- ifelse(outcome == "Patient", "PTIME_yrs", "GTIME_yrs")
    event_var <- ifelse(outcome == "Patient", "PSTATUS", "GSTATUS")
    fit <- survfit(as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ seq_flag")), data = sub_df)
    for (yr in c(1, 2)) {
      s <- summary(fit, times = yr)
      for (i in seq_along(s$strata)) {
        grp <- gsub("seq_flag=", "", names(s$strata)[i])
        results_list[[length(results_list) + 1]] <- tibble(
          `Donor Type` = dtype, Outcome = outcome,
          `Time (years)` = yr, Group = grp,
          `Survival (%)` = sprintf("%.1f (%.1f-%.1f)", s$surv[i]*100, s$lower[i]*100, s$upper[i]*100),
          `N at risk` = s$n.risk[i]
        )
      }
    }
  }
}

surv_table <- bind_rows(results_list)
cat("\n=== Supplemental Table S4: Survival Outcomes ===\n")
print(knitr::kable(surv_table, format = "simple"))
