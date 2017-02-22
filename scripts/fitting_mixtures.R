library(dplyr)
library(readr)
library(tidyr)
library(granular)
library(purrr)
library(ggplot2)

#Read in mastersizer data with mastersizer experimental design.
ms <- read_csv("data/mastersizer_input.csv")
ms <- cbind(msID = 1:nrow(ms), ms)


#Make the dataframe 'tidy' by gathering it
ms_g <- ms %>% 
  group_by(msID, StrID) %>% 
  gather(size, prop, 8:ncol(.)) %>%
  mutate_at("size", as.numeric)

ms_g %>% 
  filter(size < 1100) %>% 
ggplot(aes(log(size), prop)) + 
  geom_line(aes(group = msID))

#There are samples with a fourth peak
#These are likely artifacts from the clumping
#For our analysis we want to remove them
#Zoom in on the region of the fourth peak (~200um)
ms_g %>% 
  filter(size < 500, size > 120) %>% 
  ggplot(aes(log(size), prop)) + 
  geom_line(aes(group = msID))

#Looks like the last true curve comes down close to exp(5.3)
#get the size closest to that
limit <- unique(ms_g$size)[which.min(abs(unique(ms_g$size) - exp(5.3)))]

#Find which samples aren't at 0 at that value
ms_ep <- ms_g %>% 
  group_by(msID) %>% 
  mutate(extra_peak = ifelse(prop[size == limit] == 0, FALSE, TRUE))

ms_ep %>% 
  filter(size < 1100,
         size > 40) %>% 
  ggplot(aes(log(size), prop)) + 
  geom_line(aes(group = msID, colour = extra_peak)) +
  facet_wrap(~extra_peak)

length(unique(ms_ep$msID[ms_ep$extra_peak == 0])) #There are 835 samples. Down from 864

#Use threshold method
#C <= 1
#1 < B <= 10
#10 < A
ms_thresh <- ms_ep %>% 
  filter(!extra_peak) %>% 
  summarise(threshold_A = sum(prop[size > 11.481536]) / 100,
            threshold_B = sum(prop[size > 0.954993 & size <= 11.481536]) / 100,
            threshold_C = sum(prop[size <= 0.954993]) / 100)

#Define the component peaks, and name them. 
#This is dependent on manual inspection of the data
peaks <- exp(c(C = -0.303, B = 1.843, A = 3.059))

#Define starting values for pi - this is based on experimentation
#with different values. pi is the proportion each underlying distribution
#makes to the total
pi <- c(C = 0.02344, B = 0.39337, A = 0.58319)

#Define starting values for sigma - this is based on experimentation
#with different values. sigma is the standard deviation
#of each underlying distribution
sigma <- c(C = 0.2675, B = 0.9121, A = 0.4002)

#Use granular::mix_grp_tbl to fit the mixture distributions
#The tidy dataframe must be grouped by sample
#Uncomment the filter statement to just run a subset
#Running all values takes a long time
#These values are filtered based on absence of the fourth peak
ms_mix <- ms_ep %>% 
  filter(extra_peak == FALSE) %>% 
  group_by(StrID, msID) %>% 
  #filter(msID < 5) %>% 
  mix_grp_tbl(prop, size, peaks, pi, sigma, parallel = TRUE)

#The returned tbl has a list column for output
#Use unnest to expand it, then put into long form to reshape output stats
ms_mix_g <- ms_mix %>%
  select(mix_out) %>% 
  unnest %>% 
  gather(stat_type, stat, -StrID, -msID, -peak) %>% 
  mutate(stat_peak = paste(stat_type, peak, sep = "_")) %>% 
  select(-peak, -stat_type)

#Use spread to put each stat in its own variable
ms_wide <- ms_mix_g %>% 
  spread(stat_peak, stat)

#join proportions from mixture and threshold method
#with design for writing output

design <- read_csv("data/design.csv")
ms_out <- ms_wide %>% 
  select(msID, pi_A, pi_B, pi_C) %>% 
  left_join(ms_thresh) %>% 
  left_join(design)

#write output
#uncomment the following line to overwrite the output data
#write_csv(ms_out, "data/granule_proportions.csv") 

#write full parameters from mixture
# ms_wide %>% 
#   left_join(ms[,1:7]) %>% 
#   write_csv("data/fitted_mixtures.csv")
