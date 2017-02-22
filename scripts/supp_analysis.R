library(asreml)
library(dplyr)
dat <- read.csv("data/granule_proportions.csv") %>% 
  mutate_at(vars(-dplyr::contains("_")), as.factor)

dat$zero <- rep(0, nrow(dat))
dat$logit_piC <- log(dat$pi_C / (1 - dat$pi_C))
dat$logit_piB <- log(dat$pi_B / (1 - dat$pi_B))
dat$logit_piA <- log(dat$pi_A / (1 - dat$pi_A))
dat$C_logit <- log(dat$threshold_C / (1 - dat$threshold_C))
dat$B_logit <- log(dat$threshold_B / (1 - dat$threshold_B))
dat$A_logit <- log(dat$threshold_A / (1 - dat$threshold_A))

asrout_m1_pi3 <- asreml(A_logit ~ Type, random=~id +
                          Leftbay:zero + and(Leftbay, 0.5) + and(Rightbay, 0.5) +
                          Leftcol:zero + and(Leftcol, 0.5) + and(Rightcol, 0.5) +
                          Leftrow:zero + and(Leftrow, 0.5) + and(Rightrow, 0.5) +
                          Leftplot:zero + and(Leftplot, 0.5) + and(Rightplot, 0.5) +
                          MillRep + MillDay + MillDay:MillOrder + 
                          ExtractDay + ExtractDay:ExtractOrd + 
                          MS.Day + MS.Day:MS.Ord, data=dat, na.method.X='include')
asrout_m1_pi3 <- update(asrout_m1_pi3)

asrout_m2_pi3 <- asreml(logit_piA~ Type, random=~id +
                          Leftbay:zero + and(Leftbay, 0.5) + and(Rightbay, 0.5) +
                          Leftcol:zero + and(Leftcol, 0.5) + and(Rightcol, 0.5) +
                          Leftrow:zero + and(Leftrow, 0.5) + and(Rightrow, 0.5) +
                          Leftplot:zero + and(Leftplot, 0.5) + and(Rightplot, 0.5) +
                          MillRep + MillDay + MillDay:MillOrder + 
                          ExtractDay + ExtractDay:ExtractOrd + 
                          MS.Day + MS.Day:MS.Ord, data=dat, na.method.X='include')
asrout_m2_pi3 <- update(asrout_m2_pi3)

asrout_m1_pi2 <- asreml(B_logit ~ Type, random=~id +
                          Leftbay:zero + and(Leftbay, 0.5) + and(Rightbay, 0.5) +
                          Leftcol:zero + and(Leftcol, 0.5) + and(Rightcol, 0.5) +
                          Leftrow:zero + and(Leftrow, 0.5) + and(Rightrow, 0.5) +
                          Leftplot:zero + and(Leftplot, 0.5) + and(Rightplot, 0.5) +
                          MillRep + MillDay + MillDay:MillOrder + 
                          ExtractDay + ExtractDay:ExtractOrd + 
                          MS.Day + MS.Day:MS.Ord, data=dat, na.method.X='include')
asrout_m1_pi2 <- update(asrout_m1_pi2) 

asrout_m2_pi2 <- asreml(logit_piB ~ Type, random=~id +
                          Leftbay:zero + and(Leftbay, 0.5) + and(Rightbay, 0.5) +
                          Leftcol:zero + and(Leftcol, 0.5) + and(Rightcol, 0.5) +
                          Leftrow:zero + and(Leftrow, 0.5) + and(Rightrow, 0.5) +
                          Leftplot:zero + and(Leftplot, 0.5) + and(Rightplot, 0.5) +
                          MillRep + MillDay + MillDay:MillOrder + 
                          ExtractDay + ExtractDay:ExtractOrd + 
                          MS.Day + MS.Day:MS.Ord, data=dat, na.method.X='include')
asrout_m2_pi2 <- update(asrout_m2_pi2)

asrout_m1_pi1 <- asreml(C_logit ~ Type, random=~id +
                          Leftbay:zero + and(Leftbay, 0.5) + and(Rightbay, 0.5) +
                          Leftcol:zero + and(Leftcol, 0.5) + and(Rightcol, 0.5) +
                          Leftrow:zero + and(Leftrow, 0.5) + and(Rightrow, 0.5) +
                          Leftplot:zero + and(Leftplot, 0.5) + and(Rightplot, 0.5) +
                          MillRep + MillDay + MillDay:MillOrder + 
                          ExtractDay + ExtractDay:ExtractOrd + 
                          MS.Day + MS.Day:MS.Ord, data=dat, na.method.X='include')
asrout_m1_pi1 <- update(asrout_m1_pi1)

asrout_m2_pi1 <- asreml(logit_piC ~ Type, random=~id +
                          Leftbay:zero + and(Leftbay, 0.5) + and(Rightbay, 0.5) +
                          Leftcol:zero + and(Leftcol, 0.5) + and(Rightcol, 0.5) +
                          Leftrow:zero + and(Leftrow, 0.5) + and(Rightrow, 0.5) +
                          Leftplot:zero + and(Leftplot, 0.5) + and(Rightplot, 0.5) +
                          MillRep + MillDay + MillDay:MillOrder + 
                          ExtractDay + ExtractDay:ExtractOrd + 
                          MS.Day + MS.Day:MS.Ord, data=dat, na.method.X='include')
asrout_m2_pi1 <- update(asrout_m2_pi1)

acc <- function(asrout) {
  blup <- asrout$coefficients$random
  pev <- asrout$vcoeff$random 
  vc <- summary(asrout, nice=T)$nice
  gg <- as.character(levels(dat$id))
  diagblup.df <- data.frame(blup=blup[grep('id', names(blup))], 
                            pev=pev[grep('id', names(blup))],
                            sigg=rep(vc[['id']], each=length(gg)))
  diagblup.df$acc <- sqrt(1 - diagblup.df$pev / diagblup.df$sigg)
  return(diagblup.df)
}

acc_m1_pi1 <- acc(asrout_m1_pi1)
acc_m1_pi2 <- acc(asrout_m1_pi2)
acc_m1_pi3 <- acc(asrout_m1_pi3)
acc_m2_pi1 <- acc(asrout_m2_pi1)
acc_m2_pi2 <- acc(asrout_m2_pi2)
acc_m2_pi3 <- acc(asrout_m2_pi3)
