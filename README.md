# Method to quantify population of starch granule types

This repository accompanies the paper "Method to quantify population of starch granule types".

The repository contains:

**Data**

- Mastersizer data is in `data/mastersizer_input.csv`. 
- Experimental design including field, milling, starch extraction and mastersizer phases is in `data/design.csv`.
- Full output data for only those distributions with three peaks in `data/fitted_mixtures_noextrapeaks.csv`
- Full output data for all distributions in `data/fitted_mixtures.csv`
- Granule type proportions from both both methods, combined with design in `data/granule_proportions.csv`

**Scripts**

- The script used to generate output data from mastersizer input (`scripts/fitting_mixtures.R`)
- The script used for all statistical analysis (`scripts/supp_analysis.R`)
- The script used to generate an example figure (`distributionExample.R`)