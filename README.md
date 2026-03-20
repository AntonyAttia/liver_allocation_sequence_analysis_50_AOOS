## Liver allocation sequence analysis (AOOS) - scripts_v2

This folder contains the v2 data pipeline and the scripts used to generate the main manuscript figures and supplement (see `aoos_draft7.Rmd`).

### Pipeline (data prep)
- `00_global_config.r`: Global theme/colour scheme, anonymisation helpers, timing-cohort subset, and shared plotting utilities (sourced by the Rmd).
- `01_read_data_v2.r`: Loads raw OPTN STAR datasets (PTR, liver transplants, deceased donors, lab panels), reduces to required columns, and saves `01_raw_data_v2.rda` (study period starts at EP implementation, `2021-03-25`).
- `02_build_cohort_v2.r`: Applies exclusion criteria and builds the final cohort with 4 allocation groups: `<50`, `>=50`, `AOOS` (not EP), and `EP` (policy 9.10; code `887`). Saves `02_cohort_built_v2.rda`.
- `03_calculate_variables_v2.r`: Calculates derived variables (e.g., DRI/DSRI, MELD variants, timing/refusal measures, discarded-donor dataset, and machine-perfusion definitions) and saves `03_cohort_with_variables_v2.rda`.

### Manuscript build
- `aoos_draft7.Rmd`: Knits the paper + supplement; sources `table1.r`, all `fig*.r` scripts, and all `sf*.r` scripts and writes outputs to `figures/main` and `figures/supplement`.

### Main figures and tables
- `table1.r`: Table 1 baseline characteristics (and supporting supplemental tables S2/S3 slices).
- `fig1_bypass_trends.r`: Figure 1 bypass code usage over time.
- `fig2_refusal_timing.r`: Figure 2 refusal patterns and placement timing.
- `fig3_response_timing.r`: Figure 3 response timing relative to cross-clamp.
- `fig4_choropleth.r`: Figure 4 geographic distribution (AOOS rate and EP rate).
- `fig5_opo_center.r`: Figure 5 OPO-center pair analysis (AOOS only).
- `fig5_stratified.r`: Figure 5 stratified analysis (AOOS).
- `fig6_center_aggr.r`: Figure 6 center aggressiveness (AOOS).

### Supplemental figures
- `sf1_cohort_flow.r`: Supplemental cohort flow.
- `sf2_survival.r`: Supplemental survival analysis.
- `sf_dri_dsri.r`: Supplemental DRI/DSRI analyses.
- `sf_ep_validation.r`: Supplemental EP (policy 9.10) validation.
- `sf_machine_perfusion.r`: Supplemental machine perfusion analyses.
- `sf_utilisation.r`: Supplemental utilisation analyses.

### Citation
- `CITATION.cff`: Citation file for this analysis.
