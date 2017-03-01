[![DOI](https://zenodo.org/badge/81780225.svg)](https://zenodo.org/badge/latestdoi/81780225)

# Increased accuracy of starch granule type quantification using mixture distributions

This repository accompanies the paper "Increased accuracy of starch granule type quantification using mixture distributions".

The repository contains:

**Data**

- Mastersizer data is in `data/mastersizer_input.csv`. 
- Experimental design related variables for field, milling, starch extraction and mastersizer phases is in `data/design.csv`.
- Full output data for only those distributions with three peaks in `data/fitted_mixtures_noextrapeaks.csv`
- Full output data for all distributions in `data/fitted_mixtures.csv`
- Granule type proportions from both methods, combined with design in `data/granule_proportions.csv`

**Scripts**

- The script used to generate output data from mastersizer input (`scripts/fitting_mixtures.R`)
- The script used for fitting the statistical model in the supplementary (`scripts/supp_analysis.R`)
- The script used to generate an example figure (`scripts/distributionExample.R`)
