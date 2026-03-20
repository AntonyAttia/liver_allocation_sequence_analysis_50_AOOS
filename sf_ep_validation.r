# sf_ep_validation.r — EP Validation Summary (887 vs OPTN 6-month monitoring report)
# Confirms bypass code 887 identifies policy-compliant EP in PTR data.

EP_BYPASS_CODE <- 887L
POST_START <- as.Date("2021-03-25")
POST_END   <- as.Date("2021-09-25")
PRE_START  <- as.Date("2020-09-25")
PRE_END    <- as.Date("2021-03-24")

# Match-run level info
match_run_info <- ptr_reduced %>%
  filter(!is.na(MATCH_DATE)) %>%
  mutate(match_date = as.Date(MATCH_DATE)) %>%
  distinct(DONOR_ID, MATCH_ID_CODE, match_date)

# Post-policy match runs
post_mr <- match_run_info %>% filter(match_date >= POST_START, match_date <= POST_END)
pre_mr  <- match_run_info %>% filter(match_date >= PRE_START, match_date <= PRE_END)

# Classify post-policy by 887 and 863
post_bypass <- ptr_reduced %>%
  mutate(match_date = as.Date(MATCH_DATE)) %>%
  filter(match_date >= POST_START, match_date <= POST_END, OFFER_ACCEPT == "B")

post_has_887 <- post_bypass %>%
  filter(PRIME_OPO_REFUSAL_ID == 887) %>%
  distinct(DONOR_ID, MATCH_ID_CODE) %>% mutate(has_887 = TRUE)

post_has_863 <- post_bypass %>%
  filter(PRIME_OPO_REFUSAL_ID == 863) %>%
  distinct(DONOR_ID, MATCH_ID_CODE) %>% mutate(has_863 = TRUE)

post_classified <- post_mr %>%
  left_join(post_has_887, by = c("DONOR_ID", "MATCH_ID_CODE")) %>%
  left_join(post_has_863, by = c("DONOR_ID", "MATCH_ID_CODE")) %>%
  mutate(
    has_887 = coalesce(has_887, FALSE),
    has_863 = coalesce(has_863, FALSE),
    process = case_when(
      has_887 & !has_863  ~ "Expedited Pathway (887 only)",
      !has_887 & has_863  ~ "863 Bypass only",
      has_887 & has_863   ~ "Both (887 + 863)",
      TRUE                ~ "Standard - Not Expedited"
    )
  )

table10 <- post_classified %>%
  count(process, name = "yours") %>%
  mutate(report = case_when(
    process == "Expedited Pathway (887 only)" ~ 219L,
    process == "863 Bypass only"              ~ 344L,
    process == "Both (887 + 863)"             ~ 25L,
    process == "Standard - Not Expedited"     ~ 6614L
  ))

# Supplemental table: EP validation (2 rows only)
ep_validation_table <- table10 %>%
  filter(process %in% c("Both (887 + 863)", "Expedited Pathway (887 only)")) %>%
  select(
    Category = process,
    `Our analysis (n)` = yours,
    `OPTN report (n)` = report
  ) %>%
  arrange(match(Category, c("Both (887 + 863)", "Expedited Pathway (887 only)")))

cat("Supplemental Table. EP validation against OPTN 6-month monitoring report.\n\n")
print(ep_validation_table, n = Inf)

# Save to supplement folder
readr::write_csv(ep_validation_table, here("figures", "supplement", "Table_EP_Validation.csv"))
if (requireNamespace("gt", quietly = TRUE)) {
  gt_tbl <- ep_validation_table %>%
    gt::gt() %>%
    gt::tab_header(
      title = "EP validation against OPTN 6-month monitoring report",
      subtitle = "Counts are number of match runs in the 6-month post-policy period (March–September 2021)."
    ) %>%
    gt::tab_options(table.font.size = 14)
  html_ep <- tempfile(fileext = ".html")
  gt::gtsave(gt_tbl, filename = html_ep)
  if (requireNamespace("webshot2", quietly = TRUE)) {
    tryCatch({
      webshot2::webshot(
        url = paste0("file://", normalizePath(html_ep)),
        file = here("figures", "supplement", "Table_EP_Validation.png"),
        vwidth = 800, vheight = 400, zoom = 2
      )
    }, error = function(e) gt::gtsave(gt_tbl, here("figures", "supplement", "Table_EP_Validation.png")))
  } else {
    gt::gtsave(gt_tbl, here("figures", "supplement", "Table_EP_Validation.png"))
  }
}
cat("Saved: figures/supplement/Table_EP_Validation.csv (and .png)\n")
