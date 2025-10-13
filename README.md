# liver_allocation_sequence_analysis_50_AOOS

01_read_data.r script imports three OPTN STAR datasets—liver waitlist file, deceased donor file, and potential transplant recipient file, sets a lower date bound of 1 January 2018, and filters each dataset to retain only key variables needed for analysis.

02_build_cohort.r script constructs the analysis cohort by identifying all deceased-donor, single-organ liver transplants from 1 January 2018, linking each case to its accepted offer, and excluding multi-organ, living-donor, split-liver, duplicate-donor, and directed-donation cases. It then identifies allocation out-of-sequence (AOOS) offers, classifies transplants into allocation pathway groups (<50, ≥50, AOOS), merges full recipient and donor details, and removes records missing MELD/PELD or essential donor information. Counts are generated for visualising the cohort flow.

