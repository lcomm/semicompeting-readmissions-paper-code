# semicompeting-readmissions-paper-code
This repository contains R code to replicate all analyses for a manuscript on semicompeting risks and an application to hospital readmissions.

# Repository structure

- `data`: directory containing synthetic data (`scrData-fake.rds`) and (permissions allowing) the real data set analyzed in the paper (`scrData.rds`)
- `inst`: necessary configuration files
- `intermediates`: directory for storing saved data objects and other less critical computational outputs
- `output`: plots, tables, and other key outputs directly referenced in the manuscript
- `scripts`: numbered R scripts which recreate analysis outputs for the final manuscript version

# Available toggles

Important code toggles are in upper case near the top of the scripts. These are:
- `QUICK_RUN`: . When `QUICK_RUN = TRUE`, filenames will have a QR suffix appended to them to distinguish between 
- `USE_FAKE_DATA`: whether to use synthetic data (data/scrData-fake.rds) meant to mimic the data set documented in the paper; if not, the real data should be copied to `data/scrData.rds`
- `RUN_DISCREPANCY`: whether to calculate discrepancy metrics, as in early drafts of the paper

Simulations depend on the `batchtools` R package, with a minimal batchtools configuration file for local runs available in `inst/`. Because Bayesian methods can be time- and compute-intensive, it is strongly recommended to configure a different batchtools computation engine unless `QUICK_RUN = TRUE`.

