library(dplyr)
library(readr)
library(tidyr)
library(granular)
library(purrr)

#Read in mastersizer data with mastersizer experimental design.
ms <- read_csv("data/mastersizer.csv")
ms <- cbind(msID = 1:nrow(ms), ms)

#Define the component peaks, and name them. 
#This is dependent on manual inspection of the data
peaks <- exp(c(C = -0.303, B = 1.5, A = 3.059))

#Make the dataframe 'tidy' by gathering it
ms_g <- ms %>% 
  group_by(msID, StrID) %>% 
  gather(size, prop, 8:ncol(.)) %>%
  mutate_at("size", as.numeric)
  
#Use granular::mix_grp_tbl to fit the mixture distributions
#The tidy dataframe must be grouped by sample
#Uncomment the filter statement to just run a subset
#Running all values takes a long time
ms_mix <- ms_g %>% 
  group_by(StrID, msID) %>% 
  #filter(msID < 5) %>% 
  mix_grp_tbl(prop, size, peaks)

#The returned tbl has a list column for output
#Use unnest to expand it, then put into long form to reshape output stats
ms_mix_g <- ms_mix %>%
  unnest %>% 
  gather(stat_type, stat, -StrID, -msID, -peak) %>% 
  mutate(stat_peak = paste(stat_type, peak, sep = "_")) %>% 
  select(-peak, -stat_type)

#Use spread to put each stat in its own variable
ms_wide <- ms_mix_g %>% 
  spread(stat_peak, stat)

#join with design for writing output
ms_out <- ms_wide %>% 
  left_join(ms[,1:7])

#write output
write_csv(ms_out, "data/fitted_mixtures.csv")

#The fits can be explored with granular::ggfit_grp_tbl
#The data is easier to handle if it is all in one tidy tbl
#the distribution proportion and size data can be more 
#easily handled as list columns
ms_mix_full <- ms_g %>% 
  summarise_at(vars(size, prop), list) %>% 
  inner_join(ms_mix)
  
#The tbl needs to be grouped by sample
ms_plots <- ms_mix_full %>% 
  group_by(msID, StrID) %>% 
  ggfit_grp_tbl(mix_out, prop, size)

#Plots can be inspected
ms_plots$ggfit_plot[ms_plots$msID == 3]