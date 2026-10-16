# AOOS analysis pipeline

R scripts used to prepare the analysis cohort and generate the manuscript
figures and tables. OPTN data are not included.

## 1. Build the processed data

Put these raw OPTN SAS files in `data_raw/`:

- `ptr_li.sas7bdat`
- `liver_data.sas7bdat`
- `deceased_donor_data.sas7bdat`
- `labpanels.sas7bdat`

Then run these scripts in order:

1. `01_read_data.r`
2. `02_build_cohort.r`
3. `03_calculate_variables.r`

They write the analysis datasets into `data_processed/` (this folder is not
in the repository; you create it by running the scripts above).

Figure 4 also needs an OPO mapping file (not included).

## 2. Generate figures and tables

Once `data_processed/` exists, open `generate_figures_tables.Rmd` in RStudio
and knit it. TIFF outputs go to `figures/main/` and `figures/supplement/`.

## Main

- Table 1: `table1.r`
- Figures 1-4: `fig1_bypass_trends.r`, `fig2_refusal_timing.r`, `fig3_response_timing.r`, `fig4_choropleth.r`

## Supplement

- Tables S1-S3: `fig1_bypass_trends.r` (S1), `table1.r` (S2-S3)
- Figures S1-S8: `sf_cohort_flow.r`, Figma schematic (S2), `sf_aoos_recovery.r`, `sf_donor_placement_distributions.r`, `sf_survival.r`, `sf_center_aggressiveness.r`, `sf_opo_center_stratified.r`, `sf_utilisation.r`
