# fig4_choropleth.r: Figure 4. Choropleth maps
# Panel A: AOOS rate by OPO, Panel B: EP rate by OPO

# Fetch OPO service area polygons
base_url <- "https://gisportal.hrsa.gov/server/rest/services/Organs/OrganProcurementAndTransplantation_FS/MapServer/3/query"
url <- paste0(base_url, "?where=1=1&outFields=OPO_PROVIDER_NUM,OPO_PROVIDER_NM,OPO_STATE_ABBR&outSR=4326&f=geojson")
dsa <- st_read(url, quiet = TRUE)

# OPO-level summary with 4-group counts
opo_summary <- timing_cohort %>%
  mutate(
    OPO_CTR_CODE_current = case_when(
      OPO_CTR_CODE == "00527" ~ "11036",
      OPO_CTR_CODE == "05239" ~ "14353",
      OPO_CTR_CODE == "00093" ~ "19282",
      TRUE ~ OPO_CTR_CODE
    )
  ) %>%
  mutate(seq_flag_safe = case_when(
    seq_flag == "<50 sequence"           ~ "lt50_seq",
    seq_flag == "\u226550 sequence"      ~ "gte50_seq",
    seq_flag == "AOOS (out-of-sequence)" ~ "aoos",
    seq_flag == "EP (Policy 9.10)"       ~ "ep"
  )) %>%
  group_by(OPO_CTR_CODE = OPO_CTR_CODE_current, seq_flag_safe) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = seq_flag_safe, values_from = n, values_fill = 0) %>%
  mutate(
    total = lt50_seq + gte50_seq + aoos + ep,
    aoos_n   = aoos,
    ep_n     = ep,
    aoos_rate = if_else(total > 0, 100 * aoos_n / total, 0),
    ep_rate   = if_else(total > 0, 100 * ep_n / total, 0)
  )

# Join to geography
geo_data <- opo_summary %>%
  left_join(opo_direct_mapping, by = "OPO_CTR_CODE") %>%
  left_join(dsa, by = c("OPO_PROVIDER_NUM" = "OPO_PROVIDER_NUM")) %>%
  filter(!is.na(OPO_PROVIDER_NUM)) %>%
  st_as_sf()

geo_data_filtered <- geo_data %>% filter(!OPO_STATE_ABBR %in% c("AK", "HI", "PR"))

# Panel A: AOOS rate
p_aoos_map <- ggplot(geo_data_filtered) +
  geom_sf(aes(fill = aoos_rate), color = "black", linewidth = 0.5) +
  coord_sf(xlim = c(-130, -65), ylim = c(18, 50), expand = FALSE) +
  scale_fill_gradient2(
    low = "#24A148", mid = "#FFFFFF", high = "#D02670",
    midpoint = median(geo_data$aoos_rate, na.rm = TRUE),
    name = "AOOS Rate (%)"
  ) +
  labs(title = bquote(bold("A") ~ ". AOOS rate by OPO among transplanted recipients (excluding EP)")) +
  theme_void() +
  theme(plot.title = element_text(size = 16, hjust = 0, margin = margin(b = 10)),
        legend.position = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        plot.margin = margin(20, 20, 20, 20))

# Panel B: EP rate
p_ep_map <- ggplot(geo_data_filtered) +
  geom_sf(aes(fill = ep_rate), color = "black", linewidth = 0.5) +
  coord_sf(xlim = c(-130, -65), ylim = c(18, 50), expand = FALSE) +
  scale_fill_gradient2(
    low = "#24A148", mid = "#FFFFFF", high = "#8A3FFC",
    midpoint = median(geo_data$ep_rate, na.rm = TRUE),
    name = "EP Rate (%)"
  ) +
  labs(title = bquote(bold("B") ~ ". EP rate by OPO among transplanted recipients (Policy 9.10)")) +
  theme_void() +
  theme(plot.title = element_text(size = 16, hjust = 0, margin = margin(b = 10)),
        legend.position = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        plot.margin = margin(20, 20, 20, 20))

fig4 <- p_aoos_map / p_ep_map +
  plot_layout(nrow = 2) +
  plot_annotation(caption = "Data: HRSA OPO Service Areas. Alaska, Hawaii, and Puerto Rico excluded.",
                  theme = theme(plot.caption = element_text(size = 10, hjust = 0.5, margin = margin(t = 10))))

print(fig4)
ggsave(fig_path("figures", "main", "Figure4_Choropleth.png"), fig4, width = 14, height = 20, dpi = FIG_DPI, compression = "lzw")

# Per-OPO stats (4 groups)
dsri_by_opo <- timing_cohort %>%
  mutate(OPO_CTR_CODE_current = case_when(
    OPO_CTR_CODE == "00527" ~ "11036", OPO_CTR_CODE == "05239" ~ "14353",
    OPO_CTR_CODE == "00093" ~ "19282", TRUE ~ OPO_CTR_CODE
  )) %>%
  filter(!is.na(DSRI)) %>%
  group_by(OPO_CTR_CODE = OPO_CTR_CODE_current, seq_flag) %>%
  summarise(n = n(), median_dsri = round(median(DSRI, na.rm = TRUE), 2), .groups = "drop")

opo_summary %>%
  select(OPO_CTR_CODE, total, aoos_n, ep_n, aoos_rate, ep_rate) %>%
  arrange(desc(aoos_rate)) %>%
  print(n = Inf)

grand_totals <- opo_summary %>%
  summarise(
    total_transplants = sum(total, na.rm = TRUE),
    total_aoos = sum(aoos_n, na.rm = TRUE),
    total_ep = sum(ep_n, na.rm = TRUE),
    .groups = "drop"
  )
print(grand_totals)

total_opo <- nrow(opo_summary)
n_opo_zero_aoos <- sum(opo_summary$aoos_n == 0, na.rm = TRUE)
n_opo_zero_ep   <- sum(opo_summary$ep_n == 0, na.rm = TRUE)
n_opo_gt10_aoos <- sum(opo_summary$aoos_n > 10, na.rm = TRUE)
n_opo_gt10_ep   <- sum(opo_summary$ep_n > 10, na.rm = TRUE)
cat("Number of OPOs with 0 AOOS (no out-of-sequence): ", n_opo_zero_aoos, "/", total_opo, " OPOs\n", sep = "")
cat("Number of OPOs with 0 EP (EP volume done by remaining OPOs): ", n_opo_zero_ep, "/", total_opo, " OPOs\n", sep = "")
cat("Number of OPOs with >10 AOOS cases: ", n_opo_gt10_aoos, "/", total_opo, " OPOs\n", sep = "")
cat("Number of OPOs with >10 EP cases: ", n_opo_gt10_ep, "/", total_opo, " OPOs\n", sep = "")

opo_high <- opo_summary %>%
  mutate(high_aoos = aoos_n > 100, high_ep = ep_n > 10)
tab <- with(opo_high, table(high_ep, high_aoos))
print(addmargins(tab))
n_high_ep        <- sum(opo_high$high_ep, na.rm = TRUE)
n_high_aoos      <- sum(opo_high$high_aoos, na.rm = TRUE)
n_both           <- sum(opo_high$high_ep & opo_high$high_aoos, na.rm = TRUE)
n_high_ep_only   <- sum(opo_high$high_ep & !opo_high$high_aoos, na.rm = TRUE)
n_high_aoos_only <- sum(!opo_high$high_ep & opo_high$high_aoos, na.rm = TRUE)
cat("\n")
cat("Among OPOs with >10 EP cases (n = ", n_high_ep, "), ", n_both, " also had >100 AOOS cases.\n", sep = "")
cat("Among OPOs with >100 AOOS cases (n = ", n_high_aoos, "), ", n_both, " also had >10 EP cases.\n", sep = "")
cat("OPOs high in EP only (no high AOOS): ", n_high_ep_only, ". OPOs high in AOOS only (no high EP): ", n_high_aoos_only, ".\n", sep = "")

cat("AOOS Rate: median", round(median(opo_summary$aoos_rate, na.rm = TRUE), 1),
    "%, range", round(min(opo_summary$aoos_rate, na.rm = TRUE), 1), "-",
    round(max(opo_summary$aoos_rate, na.rm = TRUE), 1), "%\n")
cat("EP Rate: median", round(median(opo_summary$ep_rate, na.rm = TRUE), 1),
    "%, range", round(min(opo_summary$ep_rate, na.rm = TRUE), 1), "-",
    round(max(opo_summary$ep_rate, na.rm = TRUE), 1), "%\n")

print(dsri_by_opo, n = 30)
