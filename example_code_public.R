## Code for example ##

library(tidyverse)
library(survival)
library(ggthemr)

tau <- 4*365.25

# read in data ----
maindat <- read.csv("/Users/jkedwar/BoxSync/Research/HIVmeths/outcome_misclassification/survival/data/main_internal_small.csv") %>% 
  mutate(tstar = pmin(tstar, tau),  #censor at tau
         jstar = ifelse(tstar == tau, 0, jstar), 
         t = pmin(t, tau), 
         j = ifelse(t == tau, 0, j))
true <- survfit(Surv(t, j) ~ 1, data = maindat)
summary(true, times = tau, extend = T)

# read in val data ----
valdat <- maindat %>%  
  filter(r == 1)


# perform analysis ----
# estimate a and fp from val data 
est_a <- sum(valdat$tp)/sum(valdat$j)
fp <- with(valdat, ifelse((j == 0  & jstar == 1) | (tstar < t &  j == 1 & jstar ==1), 1, 0))
pt <- with(valdat, pmin(t, tstar))
est_fp_rate <- sum(fp)/sum(pt)

# estimate corrected  risk
taus <- sort(unique(maindat$tstar * maindat$jstar))
obsrisk <- survfit(Surv(tstar, jstar) ~ 1, data = data)
obs <- data.frame(time = summary(obsrisk)$time, 
                    Ftstar = 1 - summary(obsrisk)$surv)
cor_fxn <- obs %>% 
  mutate(Ftstar = risk,
         b = (1 - exp(-fp_rate * time)),
         a = est_a, 
         risk = cummax((Ftstar - b)/(a - b)))


