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

ms_ep %>% 
  filter(extra_peak == 0) %>% 
  count(msID) #There are 835 samples. Down from 864

#Define the component peaks, and name them. 
#This is dependent on manual inspection of the data
peaks <- exp(c(C = -0.303, B = 1.5, A = 3.059))

#Use granular::mix_grp_tbl to fit the mixture distributions
#The tidy dataframe must be grouped by sample
#Uncomment the filter statement to just run a subset
#Running all values takes a long time
#These values are filtered based on absence of the fourth peak
ms_mix <- ms_ep %>% 
  filter(extra_peak == FALSE) %>% 
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

ggplot(ms_out, aes(mu_C)) + geom_density()
#There is a peak around 1 which doesn't look likely.
#Can reinspect those distributions

large_c <- ms_out %>% 
  filter(mu_C > 0.5)

ms_g %>% 
  mutate(large_c = ifelse(msID %in% large_c$msID, TRUE, FALSE)) %>% 
  ggplot(aes(log(size), prop)) + 
  geom_line(aes(group = msID, colour = large_c))

tryagain <- ms_ep %>% 
  filter(extra_peak == FALSE,
         msID %in% large_c$msID) %>% 
  group_by(StrID, msID) %>% 
  mix_grp_tbl(prop, size, peaks, emnum = 20)

ms_mix_g <- tryagain %>%
  unnest %>% 
  gather(stat_type, stat, -StrID, -msID, -peak) %>% 
  mutate(stat_peak = paste(stat_type, peak, sep = "_")) %>% 
  select(-peak, -stat_type)

ms_wide <- ms_mix_g %>% 
  spread(stat_peak, stat)

large_c2 <- ms_wide %>% 
  filter(mu_C > 0.5)

tryagain2 <- ms_ep %>% 
  filter(extra_peak == FALSE,
         msID %in% large_c2$msID) %>% 
  group_by(StrID, msID) %>% 
  mix_grp_tbl(prop, size, peaks, emnum = 500)

ms_mix_g2 <- tryagain2 %>%
  unnest %>% 
  gather(stat_type, stat, -StrID, -msID, -peak) %>% 
  mutate(stat_peak = paste(stat_type, peak, sep = "_")) %>% 
  select(-peak, -stat_type)

ms_wide2 <- ms_mix_g2 %>% 
  spread(stat_peak, stat)

plot(ms_wide2$mu_C)

ms_g %>% 
  filter(msID %in% ms_wide2$msID) %>% 
  ggplot(aes(log(size), prop)) + 
  geom_line(aes(group = msID)) +
  geom_vline(xintercept = c(log(peaks), 1.3))

peaks2 <- peaks
peaks2[2] <- exp(1.3)

tryagain2 <- ms_ep %>% 
  filter(extra_peak == FALSE,
         msID %in% large_c2$msID) %>% 
  group_by(StrID, msID) %>% 
  mix_grp_tbl(prop, size, peaks2, emnum = 500)

ms_mix_g2 <- tryagain2 %>%
  unnest %>% 
  gather(stat_type, stat, -StrID, -msID, -peak) %>% 
  mutate(stat_peak = paste(stat_type, peak, sep = "_")) %>% 
  select(-peak, -stat_type)

ms_wide2 <- ms_mix_g2 %>% 
  spread(stat_peak, stat)

plot(ms_wide2$mu_C)

#There are several that still won't work.
#Take a single example
dist <- ms_g$prop[ms_g$msID == 421]
ps <- unique(ms_g$size)
index_start <- min(which(dist!=0)) # to remove trailing and leading 0 entries.
index_end <- max(which(dist!=0))
dat <- data.frame(ps = log(ps[index_start:index_end]),
                  dist = dist[index_start:index_end])
heights_d <- granular:::get_heights(dat[["dist"]], dat[["ps"]], log(peaks))
comp_weights <- heights_d/(sum(heights_d))

comp_sds <- rep(diff(range(dat$ps))/3, times = 3)

initial_values <- mixdist::mixparam(log(peaks), comp_sds, comp_weights)

temp <- mix(dat, initial_values) # what! that worked
#try it in grp_tbl

ms_g %>% 
  filter(msID == 142) %>% 
  mix_grp_tbl(prop, size, peaks, emnum = 1) -> temp2

sub <- ms_g %>% 
  filter(msID == 142)

#try mix_dist
mix_dist_142 <- mix_dist(sub$prop, sub$size, peaks)
#running with emnum == 1 works :(

tryagain2 <- ms_ep %>% 
  filter(extra_peak == FALSE,
         msID %in% large_c2$msID) %>% 
  group_by(StrID, msID) %>% 
  mix_grp_tbl(prop, size, peaks2, emnum = 1)

ms_mix_g2 <- tryagain2 %>%
  unnest %>% 
  gather(stat_type, stat, -StrID, -msID, -peak) %>% 
  mutate(stat_peak = paste(stat_type, peak, sep = "_")) %>% 
  select(-peak, -stat_type)

ms_wide2 <- ms_mix_g2 %>% 
  spread(stat_peak, stat)

plot(ms_wide2$mu_C)

#titrate emnum
emnum_list <- c(1:10, seq(from = 12, to = 30, by = 2),
                seq(from = 35, to = 70, by = 5), c(80, 90, 100))

test <- map(emnum_list[1], ~ ms_ep %>% 
              filter(extra_peak == FALSE,
                     msID %in% large_c2$msID) %>% 
              group_by(StrID, msID) %>% 
              mix_grp_tbl(prop, size, peaks2, emnum = .x))

#test2 <- setNames(test, emnum_list) %>% 
  bind_rows(test, .id = "emnum") %>% 
  unnest %>% 
  mutate_at("emnum", as.numeric)

ggplot(test2, aes(emnum, mu)) + 
  geom_line(aes(group = msID)) +
  facet_wrap(~peak)

ggplot(test2, aes(emnum, mu)) + 
  geom_line(aes(group = peak, colour = peak)) +
  facet_wrap(~msID)
#repeat with the easy ones
test_easy <- map(emnum_list, ~  {
  print(.x)
  ms_ep %>% 
    filter(extra_peak == FALSE,
           !msID %in% large_c2$msID) %>% 
    filter(msID < 11) %>% 
    group_by(StrID, msID) %>% 
    mix_grp_tbl(prop, size, peaks2, emnum = .x)
})

test_easy2 <- setNames(test_easy, emnum_list) %>% 
  bind_rows(test_easy, .id = "emnum") %>% 
  unnest %>% 
  mutate_at("emnum", as.numeric)

ggplot(test2, aes(emnum, mu)) + 
  geom_line(aes(group = msID)) +
  facet_wrap(~peak)

ggplot(test_easy2, aes(emnum, mu)) + 
  geom_line(aes(group = peak, colour = peak)) +
  facet_wrap(~msID)

#What about if we change the dist slightly?
hard <- ms_ep %>% 
  filter(extra_peak == FALSE,
         msID %in% large_c2$msID)

hard <- hard %>% 
  mutate(prop_mod = prop * rnorm(n(), 1, 0.01),
         prop_mod.1 = prop * rnorm(n(), 1, 0.1),
         prop_mod10 = prop * 10)

hard_mix <- map(seq(from = 1, to = 10, by = 2), ~ {
  print(.x)
  hard %>% 
    group_by(StrID, msID) %>% 
    mix_grp_tbl(prop_mod, size, peaks2, emnum = .x)
})

hard_df <- setNames(hard_mix, seq(from = 1, to = 10, by = 2)) %>% 
  bind_rows(hard_mix, .id = "emnum") %>% 
  unnest %>% 
  mutate_at("emnum", as.numeric)

ggplot(hard_df, aes(emnum, mu)) + 
  geom_line(aes(group = msID)) +
  facet_wrap(~peak)

ggplot(hard_df, aes(emnum, mu)) + 
  geom_line(aes(group = peak, colour = peak)) +
  facet_wrap(~msID)

hard_mix2 <- map(seq(from = 1, to = 10, by = 2), ~ {
  print(.x)
  hard %>% 
    group_by(StrID, msID) %>% 
    mix_grp_tbl(prop_mod.1, size, peaks2, emnum = .x)
})

hard_df2 <- setNames(hard_mix2, seq(from = 1, to = 10, by = 2)) %>% 
  bind_rows(hard_mix2, .id = "emnum") %>% 
  unnest %>% 
  mutate_at("emnum", as.numeric)

ggplot(hard_df2, aes(emnum, mu)) + 
  geom_line(aes(group = msID)) +
  facet_wrap(~peak)

ggplot(hard_df2, aes(emnum, mu)) + 
  geom_line(aes(group = peak, colour = peak)) +
  facet_wrap(~msID)

hard_mix10 <- map(setNames(seq(from = 1, to = 10, by = 3),
                         seq(from = 1, to = 10, by = 3)),
                ~ {
                  print(.x)
                  hard %>% 
                    group_by(StrID, msID) %>% 
                    mix_grp_tbl(prop_mod10, size, peaks2, emnum = .x)
                })

hard_df10 <- hard_mix10 %>% 
  bind_rows(.id = "emnum") %>% 
  unnest %>% 
  mutate_at("emnum", as.numeric)

ggplot(hard_df10, aes(emnum, mu)) + 
  geom_line(aes(group = msID)) +
  facet_wrap(~peak)

ggplot(hard_df10, aes(emnum, mu)) + 
  geom_line(aes(group = peak, colour = peak)) +
  facet_wrap(~msID)

#Let's compare a single hard one with all the easy ones, see what's different
#I'll take 493
#First the easy ones
easy <- ms_ep %>% 
  filter(extra_peak == FALSE,
         !msID %in% large_c2$msID)

comp493 <- hard %>% 
  filter(msID == 493) %>% 
  ungroup %>% 
  select(prop493 = prop, size) %>% 
  full_join(easy) %>% 
  filter(prop493 > 0 & prop > 0)

diff493 <- comp493 %>% 
  group_by(msID) %>% 
  mutate(prop_diff = abs(prop - prop493)) %>% 
  summarise(sum_diff = sum(prop_diff),
            max_diff = max(prop_diff))

diff493 %>% arrange(sum_diff)
diff493 %>% arrange(desc(sum_diff))
ms_g %>% 
  filter(msID %in% c(hard$msID, 541, 309),
         size < 100) %>% 
  ggplot(aes(log(size), prop)) + 
  geom_line(aes(group = msID, colour = factor(msID)))
#writing temp output
#need ms_wide2, ms_wide that aren't in ms_wide2 and ms_out that aren't in either

ms_out2 <- ms_wide2 %>% 
  bind_rows(ms_wide %>% 
              filter(!msID %in% ms_wide2$msID)) %>% 
  inner_join(ms[,1:7])
ms_out_all <- ms_out %>% 
  filter(!msID %in% ms_out2$msID) %>% 
  bind_rows(ms_out2) %>% 
  arrange(msID)
write_csv(ms_out_all, "data/fitted_mixtures_noextrapeak.csv") 
#write output
#uncomment this line to overwrite the output data
write_csv(ms_out, "data/fitted_mixtures_noextrapeak.csv") 

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
ms_plots$ggfit_plot[ms_plots$msID == 27]
