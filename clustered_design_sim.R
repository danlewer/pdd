setwd("C:/Users/rmhidle/OneDrive - University College London/pdd/design")

# -------------------------------
# general functions and libraries
# -------------------------------

library(lme4) # for REML/ML mixed effects regression
library(meta) # for meta analysis of matched analyses

# convert risk to odds

risk2odds <- function (risk) risk / (1-risk)
odds2risk <- function (odds) odds / (1+odds)
or2rr <- function (OR, p0) OR / ((1 - p0) + (p0 * OR))
rr2or <- function(RR, p0) RR * (1 - p0) / (1 - RR * p0)

# lots of letters

ll <- do.call(paste0, expand.grid(LETTERS, LETTERS))

# ----------------------------
# force-level reoffending rate
# ----------------------------

fs <- read.csv("https://raw.githubusercontent.com/danlewer/pdd/main/force_stats.csv")
hist(fs$drug_reoffending)
mean(fs$drug_reoffending)
sd(fs$drug_reoffending)

m <- mean(fs$drug_reoffending) # 0.21
s <- sd(fs$drug_reoffending) # 0.03

x <- seq(0.14, 0.28, 0.001)
y <- dnorm(x, mean = m, sd = s)

png('reoffending_hist.png', height = 5, width = 6, units = 'in', res = 300)
hist(fs$drug_reoffending, breaks = 7:14/50, xlab = 'Reoffending rate', main = NA)
lines(x, y, col = 'red')
abline(v = c(m, m - s, m + s), col = 'red')
dev.off()

# expected odds ratio with risk difference of 3ppts
# rr is 0.18/0.21 = 0.857
rr2or(0.18/0.21, p0 = 0.21) # OR = 0.825784

# -------------
# police forces
# -------------

pf <- data.frame(
  force = ll[1:16],
  pdd = rep(0:1, each = 8),
  n = rep(1000, 16)
)

# --------------------
# simulation functions
# --------------------

# simulate data

sim <- function (
    
  outcome_baseline = 0.21,
  outcome_force_sd = 0.03,
  pdd_risk_difference = 0.03,
  pdd_force_sd = 0.007) {
  
  force_baseline_risk <- rnorm(nrow(pf), mean = outcome_baseline, sd = outcome_force_sd)
  force_baseline_odds <- risk2odds(force_baseline_risk)
  force_pdd_risk_diff <- rnorm(nrow(pf), mean = pdd_risk_difference, sd = pdd_force_sd)
  force_pdd_rr <- (force_baseline_risk - force_pdd_risk_diff) / force_baseline_risk
  force_pdd_or <- rr2or(force_pdd_rr, force_baseline_risk)
  
  d <- data.frame(
    force = rep(pf$force, pf$n),
    pdd = rep(pf$pdd, pf$n),
    force_baseline_risk = rep(force_baseline_risk, pf$n),
    baseline_odds = rep(force_baseline_odds, pf$n),
    force_or = rep(force_pdd_or, pf$n))
  
  d$force_or[d$pdd == 0] <- 1
  d$target_odds <- d$baseline_odds * d$force_or 
  d$target_risk <- odds2risk(d$target_odds)
  d$outcome <- rbinom(nrow(d), 1, d$target_risk)
  
  return (d[, c('force', 'pdd', 'force_baseline_risk', 'outcome')])}

d <- sim() # example simulation
tapply(d$outcome, d$pdd, mean)

# regression-based adjustment (one-step)
# note nAGQ = 0 and confint.merMod method = 'Wald' approximations

f <- function(B = 100, mf = 'glm', form = 'outcome ~ pdd', actual = log(rr2or(0.18/0.21, p0 = 0.21)), ...) {
  Bs <- sapply(1:B, function (x) {
    if (x %% 10 == 0) print (x)
    d <- sim(...)
    m <- if (mf == 'glmer') {
      get(mf)(form, data = d, family = 'binomial', nAGQ = 0)
    } else {
      get(mf)(form, data = d, family = 'binomial')
    }
    c(summary(m)$coef[2,1], suppressMessages(confint(m, parm = 'pdd', method = 'Wald')), summary(m)$coef[2, 4]) # method is ignored for glm
  })
  c(mean_coef = mean(Bs[1,]),
    OR = exp(mean(Bs[1,])),
    coverage = mean(actual > Bs[2,] & actual < Bs[3,]),
    power = mean(Bs[4,] < 0.05))
}

# matched analysis (two-step)

fm <- function (B = 100, pdd_forces = pf$force[pf$pdd == 1], mf = 'glm', form = 'outcome ~ pdd + force_baseline_risk', actual = log(rr2or(0.18/0.21, p0 = 0.21)), ...) {
  Bs <- sapply(1:B, function (z) {
    if (z %% 10 == 0) print (z)
    d <- sim(...)
    e <- sapply(pdd_forces, function (x) {
      cases <- d[d$force == x,]
      pool <- d[d$pdd == 0,] # this is where exact matching would happen
      controls <- pool[sample(1:nrow(pool), nrow(cases) * 2, replace = T),]
      d2 <- rbind(cases, controls)
      m <- if (mf == 'glmer') {
        get(mf)(form, data = d, family = 'binomial', nAGQ = 0)
      } else {
        get(mf)(form, data = d, family = 'binomial')
      }
      summary(m)$coef[2, c(1, 2, 4)]
    })
    ma <- metagen(TE = e[1,], seTE = e[2,])
    c(FEpoint = ma$TE.common, FElower = ma$lower.common, FEupper = ma$upper.common, FEp = ma$pval.common,
      REpoint = ma$TE.random, RElower = ma$lower.random, REupper = ma$upper.random, REp = ma$pval.random)
  })
  c(mean_coef = mean(Bs[1,]),
    OR = exp(mean(Bs[1,])),
    coverage = mean(actual > Bs[2,] & actual < Bs[3,]),
    power = mean(Bs[4,] < 0.05))
}

# -------
# results
# -------

set.seed(378)

# 1: no cluster effects (coverage should be 0.95; power is good)
# --------------------------------------------------------------

r1 <- f(B = 1000, pdd_force_sd = 0, outcome_force_sd = 0)

# 2: cluster effects but no accounting for cluster effects (coverage is poor)
# ---------------------------------------------------------------------------

r2 <- f(B = 1000, pdd_force_sd = 0, pdd_risk_difference = 0.01, actual = log(rr2or(0.20/0.21, 0.21))) 
r3 <- f(B = 1000, pdd_force_sd = 0, pdd_risk_difference = 0.03) 
r4 <- f(B = 1000, pdd_force_sd = 0, pdd_risk_difference = 0.06, actual = log(rr2or(0.15/0.21, 0.21))) # very low coverage, very high power

# 3: random intercepts for police forces - coverage is good but power is low (RE washes out PDD)
# ----------------------------------------------------------------------------------------------

r5 <- f(B = 1000, mf = 'glmer', form = 'outcome ~ pdd + (1|force)', pdd_force_sd = 0, pdd_risk_difference = 0.01, actual = log(rr2or(0.20/0.21, 0.21)))
r6 <- f(B = 1000, mf = 'glmer', form = 'outcome ~ pdd + (1|force)', pdd_force_sd = 0, pdd_risk_difference = 0.03)
r7 <- f(B = 1000, mf = 'glmer', form = 'outcome ~ pdd + (1|force)', pdd_force_sd = 0, pdd_risk_difference = 0.06, actual = log(rr2or(0.15/0.21, 0.21)))

# versions with baseline_risk fail to estimate random intercepts - but work ok just with the force_baseline_risk (as this does not seek to estimate)

r5a <- f(B = 1000, mf = 'glmer', form = 'outcome ~ pdd + force_baseline_risk + (1|force)', pdd_force_sd = 0, pdd_risk_difference = 0.01, actual = log(rr2or(0.20/0.21, 0.21)))
r6a <- f(B = 1000, mf = 'glmer', form = 'outcome ~ pdd + force_baseline_risk + (1|force)', pdd_force_sd = 0, pdd_risk_difference = 0.03)
r7a <- f(B = 1000, mf = 'glmer', form = 'outcome ~ pdd + force_baseline_risk + (1|force)', pdd_force_sd = 0, pdd_risk_difference = 0.06, actual = log(rr2or(0.15/0.21, 0.21)))

# 5: matched analysis (1:2)
# -------------------------

# not adjusting for force baseline risk (biased)

r8 <- fm(B = 1000, form = 'outcome ~ pdd', pdd_force_sd = 0, pdd_risk_difference = 0.01, actual = log(rr2or(0.20/0.21, p0 = 0.21)))
r9 <- fm(B = 1000, form = 'outcome ~ pdd', pdd_force_sd = 0, pdd_risk_difference = 0.03, actual = log(rr2or(0.18/0.21, p0 = 0.21)))
r10 <- fm(B = 1000, form = 'outcome ~ pdd', pdd_force_sd = 0, pdd_risk_difference = 0.06, actual = log(rr2or(0.15/0.21, p0 = 0.21)))

# adjusting for force baseline risk

r11 <- fm(B = 1000, pdd_force_sd = 0, pdd_risk_difference = 0.01, actual = log(rr2or(0.20/0.21, p0 = 0.21)))
r12 <- fm(B = 1000, pdd_force_sd = 0, pdd_risk_difference = 0.03, actual = log(rr2or(0.18/0.21, p0 = 0.21)))
r13 <- fm(B = 1000, pdd_force_sd = 0, pdd_risk_difference = 0.06, actual = log(rr2or(0.15/0.21, p0 = 0.21)))

# including random intercept (this takes a long time as has to do RE model for each intervention force)

r14 <- fm(B = 1000, pdd_force_sd = 0, mf = 'glmer', form = 'outcome ~ pdd + (1|force)', pdd_risk_difference = 0.01, actual = log(rr2or(0.20/0.21, p0 = 0.21)))
r15 <- fm(B = 1000, pdd_force_sd = 0, mf = 'glmer', form = 'outcome ~ pdd + (1|force)', pdd_risk_difference = 0.03, actual = log(rr2or(0.18/0.21, p0 = 0.21)))
r16 <- fm(B = 1000, pdd_force_sd = 0, mf = 'glmer', form = 'outcome ~ pdd + (1|force)', pdd_risk_difference = 0.06, actual = log(rr2or(0.15/0.21, p0 = 0.21)))

# summary

res <- rbind(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16)
write.csv(res, 'design_sims.csv')

# notes

# ONS publishes force level reoffending for drug offenses - 21%, it does vary by police forces - with sd of 0.03 across forces
# in primary analysis of diversion policies, the intervention is allocated at the force level, and you cannot ignore police forces
# if you ignore police forces and pretend the intervention is allocated to individuals, the result will be highly biased (see low coverage)
# you cannot adjust for baseline reoffending because the interventions precede the study baseline
# random effects will wash out the intervention (low power in random intercepts model)
# you can adjust for stable factors that may determine reoffending, eg. deprivation, ?drug related offending, ?drug use
