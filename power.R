library(epiR) # for epi.sscohortt

options(scipen = 999) # just for figure formatting

# areas

# Thames Valley (Hampshire)
# Durham (Humberside)
# West Midlands (Greater Manchester)

# inputs
# ......

n_intervention <- 1000L
n_comparison <- 1000L

# follow-up duration - minimum 12 months, maximum 24 months (for participants recruited at the start)
# assume matched on cohort entry

follow_up_min <- 365L
follow_up_max <- 730L

# alpha (for power)

alpha <- 0.05

# exposure groups: (1a) PDD area: diverted; (1b) PDD: not diverted; (2) comparison
# primary outcome: admissions related to acute mental health and drug or alcohol problems - HES (PHE)

# outcome frequency (in events per person-year)
# using data from a study of c.100,000 people with records of illicit opioid use in primary care
# https://wellcomeopenresearch.org/articles/5-282
# (calculated specifically for this power estimate)
# rate of admission due to mental health / drugs & alcohol = 0.28 per person-year
# this compares to 0.04 in the general population
# also note that 58% had such an admission in mean 9 years' follow-up

outcome_frequency <- 0.28

# simulate data
# .............

daily_frequency <- outcome_frequency / 365.25

sim <- function (effect = 0.7) {
  follow_up_intervention <- sample(x = follow_up_min:follow_up_max, size = n_intervention, replace = T)
  events_intervention <- rbinom(n = n_intervention, size = follow_up_intervention, prob = daily_frequency * effect)
  events_comparison <- rbinom(n = n_intervention, size = follow_up_intervention, prob = daily_frequency)
  data.frame(exposure = c(rep(0, n_comparison), rep(1, n_intervention)),
             follow_up = rep(follow_up_intervention, times = 2),
             outcome = c(events_comparison, events_intervention))}

# example simulation

set.seed(304)
example_sim <- sim(effect = 0.5)
sapply(split(example_sim, f = example_sim$exposure), function(x) with(x, c(
  events = sum(outcome), 
  person_years = sum(follow_up) / 365, 
  rate = sum(outcome) / (sum(follow_up)/365)))) # rate is approx. half in unexposed group

# power calculation function
# ..........................

# uses a GLM to calculate the rate ratio

power <- function (..., B = 1000L) { # B is the number of simulation. ... provided for effect estimate
  d <- lapply(seq_len(B), function (x) sim (...))
  m <- lapply(d, function (x) glm(outcome ~ exposure, data = x, family = 'poisson'))
  t(sapply(m, function (x) summary(x)$coef[2, c(1, 4)]))
}

# example power calculation for IRR = 0.9 (i.e. rate of outcome in diversion area is 0.9x that in comparison area)

set.seed(14)
power7 <- power(effect = 0.9, B = 1000L)
hist(power7[,1][power7[,1] < 0], xlab = 'log IRR', ylab = 'Number of simulations', main = 'Effect of diversion in 1000 simulations', xlim = c(-0.5, 0.2), breaks = seq(-0.5, 0, 0.02), col = '#D4E157', border = '#004D40')
hist(power7[,1][power7[,1] >= 0], add = T, breaks = seq(0, 0.2, 0.02), col = '#4DB6AC', border = '#004D40')
abline(v = 0)
mean(power7[,2] < alpha) # power = 31%

# compare to epi.sscohortt (unsure of validity of this formula)

epi.sscohortt(irexp0 = outcome_frequency, irexp1 = outcome_frequency * 0.9, FT = 1.5, n = 2000, power = NA)

# multiple power calculations for different effect sizes

logIRR <- seq(-0.5, 0, 0.01)
effects <- exp(logIRR)
ns <- seq(500, 2500, 100)

set.seed(103)

t0 <- proc.time()

power_calcs <- lapply(ns, function(y) {
  sapply(effects, function (x) { # prints off value of effect size currently being simulated
    n_intervention <<- y
    n_comparison <<- y
    print(paste0(y, '; ', x))
    d <- power(effect = x, B = 1000L)
    c(target = x, mean_IRR = exp(mean(d[,1])), power = mean(d[,2] < alpha))
  })
})

t1 <- proc.time()
t1 - t0

power_calcs <- data.frame(t(power_calcs))
with(power_calcs, plot(log(target), power, type = 'b', pch = 19, ylim = c(0, 1), xlim = range(logIRR), xlab = 'Effect size (log IRR)', ylab = 'Power'))

# using epiR::epi.sscohortt

epi.sscohortt(irexp0 = outcome_frequency, irexp1 = outcome_frequency * 0.8, FT = 1.5, power = 0.8)
