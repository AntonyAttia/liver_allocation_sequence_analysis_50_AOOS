# sf1_cohort_flow.r — SF1: Cohort Flow Diagram (CONSORT-style)
# Run from Rmd (uses existing counts/discarded_complete) or source standalone (loads from saved files).

# Always attach packages this script needs (safe if already loaded)
for (pkg in c("here", "glue", "scales", "DiagrammeR", "DiagrammeRsvg", "rsvg")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Please install.packages('", pkg, "')")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
if (!exists("counts")) {
  counts <- readRDS(here("data_processed", "02_cohort_counts_v2.rds"))
  if (!exists("discarded_complete")) load(here("data_processed", "03_cohort_with_variables_v2.rda"))
}

comma <- scales::comma

# Exclusion cascade from cohort building (02_build_cohort_v2.R)
n_init    <- comma(counts$n_initial)
n_join    <- comma(counts$n_after_join)
n_mo      <- comma(counts$n_after_multiorg)
n_liv     <- comma(counts$n_after_living)
n_spl     <- comma(counts$n_after_split)
n_ded     <- comma(counts$n_after_dedup)
n_dir     <- comma(counts$n_after_directed)
n_fin     <- comma(counts$n_final)

e_join <- comma(counts$n_initial - counts$n_after_join)
e_mo   <- comma(counts$n_after_join - counts$n_after_multiorg)
e_liv  <- comma(counts$n_after_multiorg - counts$n_after_living)
e_spl  <- comma(counts$n_after_living - counts$n_after_split)
e_ded  <- comma(counts$n_after_split - counts$n_after_dedup)
e_dir  <- comma(counts$n_after_dedup - counts$n_after_directed)
e_mis  <- comma(counts$n_after_directed - counts$n_final)

# Timing cohort (00_global_config.r applies timing plausibility filters on final_cohort)
# This is what Table 1 and most figures use, so the flow diagram should end here
if (exists("timing_cohort") && is.data.frame(timing_cohort) && nrow(timing_cohort) > 0) {
  n_timing    <- nrow(timing_cohort)
  e_timing    <- comma(counts$n_final - n_timing)
  n_timing_f  <- comma(n_timing)
  tc_breakdown <- table(timing_cohort$seq_flag)
  n_lt50  <- comma(tc_breakdown["<50 sequence"])
  n_gte50 <- comma(tc_breakdown["\u226550 sequence"])
  n_aoos  <- comma(tc_breakdown["AOOS (out-of-sequence)"])
  n_ep    <- comma(tc_breakdown["EP (Policy 9.10)"])
  has_timing_step <- TRUE
} else {
  fc_breakdown <- counts$final_cohort_breakdown
  e_timing    <- "0"
  n_timing_f  <- n_fin
  n_lt50  <- comma(fc_breakdown["<50 sequence"])
  n_gte50 <- comma(fc_breakdown["\u226550 sequence"])
  n_aoos  <- comma(fc_breakdown["AOOS (out-of-sequence)"])
  n_ep    <- comma(fc_breakdown["EP (Policy 9.10)"])
  has_timing_step <- FALSE
}

# Label with actual ≥ character so Graphviz displays it correctly
gte50_label <- "\u226550 sequence"

# Discarded (not transplanted) section — use when available from 03_cohort_with_variables_v2.rda
if (exists("discarded_complete") && is.data.frame(discarded_complete) && nrow(discarded_complete) > 0) {
  n_disc_aoos       <- comma(sum(discarded_complete$allocation_type == "AOOS (out-of-sequence)", na.rm = TRUE))
  n_disc_sequential <- comma(sum(discarded_complete$allocation_type == "Sequential", na.rm = TRUE))
  n_disc_ep         <- comma(sum(discarded_complete$allocation_type == "EP (Policy 9.10)", na.rm = TRUE))
  n_disc_total      <- comma(nrow(discarded_complete))
  n_disc_dbd        <- comma(sum(discarded_complete$donor_type == "DBD", na.rm = TRUE))
  n_disc_dcd        <- comma(sum(discarded_complete$donor_type == "DCD", na.rm = TRUE))
  has_discarded     <- TRUE
} else {
  n_disc_aoos <- n_disc_sequential <- n_disc_ep <- n_disc_total <- "0"
  n_disc_dbd <- n_disc_dcd <- "0"
  has_discarded <- FALSE
}

# Draft5/6 style: breakdown text inside final cohort node
breakdown_text <- paste0(
  "In-Sequence (<50): ", n_lt50, "  |  ", gte50_label, ": ", n_gte50, "\\n",
  "AOOS: ", n_aoos, "  |  EP: ", n_ep
)
never_accepted_text <- paste0(
  "AOOS: ", n_disc_aoos, "\\n",
  "In-Sequence: ", n_disc_sequential, "\\n",
  "EP: ", n_disc_ep, "\\n",
  "Total: ", n_disc_total
)

consort_dot <- glue('
digraph consort {{
  graph [rankdir=TB, splines=ortho, nodesep=0.35, ranksep=0.45]
  node [shape=box, style="rounded,filled", fillcolor="white", color="#666666", fontname="Arial", fontsize=10]
  edge [color="#666666", arrowsize=0.6]

  # Transplant flow (left) — no connection to discarded
  A [label="All liver transplants in study period\\nN = {n_init}", fillcolor="#F8F8F8"]
  B [label="With accepted PTR record\\nN = {n_join}"]
  C [label="After multi-organ exclusion\\nN = {n_mo}"]
  D [label="After living-donor exclusion\\nN = {n_liv}"]
  E [label="After split-liver exclusion\\nN = {n_spl}"]
  F [label="After PTR duplicate resolution\\nN = {n_ded}"]
  G [label="After directed-donation exclusion\\nN = {n_dir}"]
  G2 [label="After missing data exclusion\\nN = {n_fin}"]
  H [label="Final analysis cohort\\nN = {n_timing_f}\\n\\n{breakdown_text}", fillcolor="#EFEFEF"]

  A -> B -> C -> D -> E -> F -> G -> G2 -> H

  # Exclusions (note shape, to the right; not part of flow)
  node [shape=note, fillcolor="white", color="#AAAAAA"]
  X1 [label="Excluded: No PTR match\\nN = {e_join}", fillcolor="#FFE6E6", color="#CC0000"]
  X2 [label="Excluded: Multi-organ\\nN = {e_mo}"]
  X3 [label="Excluded: Living donor\\nN = {e_liv}"]
  X4 [label="Excluded: Split liver\\nN = {e_spl}"]
  X5 [label="Excluded: PTR duplicates\\nN = {e_ded}"]
  X6 [label="Excluded: Directed\\nN = {e_dir}"]
  X7 [label="Excluded: Missing\\ndonor data\\nN = {e_mis}"]
  X8 [label="Excluded: Timing\\nplausibility\\nN = {e_timing}"]

  A -> X1 [arrowhead=none, color="#CC0000", style=bold]
  B -> X2 [arrowhead=none, color="#BBBBBB"]
  C -> X3 [arrowhead=none, color="#BBBBBB"]
  D -> X4 [arrowhead=none, color="#BBBBBB"]
  E -> X5 [arrowhead=none, color="#BBBBBB"]
  F -> X6 [arrowhead=none, color="#BBBBBB"]
  G -> X7 [arrowhead=none, color="#BBBBBB"]
  G2 -> X8 [arrowhead=none, color="#BBBBBB"]

  # Standalone discarded flow (separate root; not from transplant box)
  node [shape=box, style="rounded,filled", fillcolor="#FFF4E6", color="#666666"]
  J [label="Offers not resulting in transplant\\nN = {n_disc_total}"]
  I2 [label="With donor type data\\nN = {n_disc_total}\\nDBD: {n_disc_dbd}  DCD: {n_disc_dcd}"]
  I [label="{never_accepted_text}"]

  J -> I2 -> I

  {{rank=same; H; I}}
}}
')

svg_raw <- DiagrammeR::grViz(consort_dot) %>% DiagrammeRsvg::export_svg()
rsvg::rsvg_png(charToRaw(svg_raw), here("figures", "supplement", "FigureSF1_Cohort_Flow.png"), width = 2400)
cat("SF1 saved.\n")

# Also render in Rmd
DiagrammeR::grViz(consort_dot)
