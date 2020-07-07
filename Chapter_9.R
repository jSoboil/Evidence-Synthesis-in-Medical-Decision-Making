# ========================================================================================
# Mixed treatment and indirect treatment comparisons ----------------------
# ========================================================================================
# Why go beyond 'direct' head-to-head trials? Several reasons can be advanced for taking
# a wider view of what the legitimate evidence base should be.

# First, it may be that there are no A vs B trials, but instead an 'indirect' estimate can
# be formed from the results of A vs C and B vs C trials:

# d[ind[AB]] = d[dir[AC]] - d[dir[BC]]

# A second reason might be that, even if direct AB evidence exists, it may be sparse; the 
# volume of indirect evidence can be much greater. This is, in fact, a very common 
# situation. 

# The Cochrane Collaboration has suggested that direct and indirect estimates should be 
# kept dinstinct. However, once the legitamacy of indirect evidence has been conceeded 
# when direct data are lacking or sparse, it is difficult to justify excluding it all. And
# even if the direct evidence is sparse, why not pool the direct and indirect using for 
# example an inverse-weighted average? It could even be the case that the quality of 
# indirect evidence is better than the available direct evidence!

# This pooling of direct and indirect evidence is referred to as a mixed treatment 
# comparison. Using P (precision) to denote the reciprocal of the variance of the 
# estimates:

# d[pooled[AB]] = ((P[dir[AB]] * d[dir[AB]]) + (P[dir[AB]] * d[dir[AB]])) / 
#                                           (P[dir[AB]] + P[Indir[AB]])

# Such pooling can be seen to be no more than an extension of the widely accepted rubric 
# that all the available evidence should be marshalled, to increase precision and to avoid
# selection biases. But there is a still more cogent reason for turning to mixed treatment 
# comparison synthesis.

# Clinicians and national decision-making bodies frequently have to make a choice between 
# several different alternative treatments: in this situation there seems to be no 
# alternative other than to combine the data on all the pairwise comparisons within a 
# single unified analysis, producing an internally consistent set of relative treatment 
# effects.

# ========================================================================================
# A fixed treatment effects model for MTC ---------------------------------
# ========================================================================================
# We will begin with trials generating binomial data. The models are, in fact, recognisable
# as logistic regression models; however, it is useful to describe the model making use of 
# the terminology of the Confidence Profile Method.

# Suppose we have four treatments A, B, C and D, where A is the ‘reference’ or standard 
# treatment. We define three treatment effect parameters representing the Log Odds Ratios 
# (LORs) of B, C and D relative to reference treatment A: d[AB], d[AC], d[AD]. These are 
# the basic parameters, and they will need to be given priors. 
# There are  4 * (4 - 1) / 2 = 6 potential constrasts as if there are K treatments, there 
# are K(K - 1)/2 possible pairwise comparisons: so with six treatments there are 15 
# contrasts of potential interest to a decision maker). In the data set, there is at least
# one trial direclty informing every contrast. The three remaining contrasts d[BC], d[BD],
# d[CD], are represented by functional parameters, and are therefore defined in terms of 
# the basic parameters:

# d[BC] = d[AC] - d[AB]
# d[BD] = d[AD] - d[AB]
# d[CD] = d[AD] - d[AC]

# These consistency equations capture the idea behind MTC: simply put, if (b - a) = 2 and
# (c - a) = 3, then (c - b) must be = 1. If we abandon the consistency relations we would 
# revert to a model in which there were six LORs which are quite unrelated to each other.

# First, some notation: we will adopt the convention that d[XY] is the effect of Y relative
# to X, and we will always express these relative effects with Y alphabetically following
# X, because d[YX] = -d[XY]. The full Fixed Effect model will be as follows, for treatment 
# k in trial j:

#    r[jk] ~ Binomial(p[jk], n[jk])
#    logit(p[jk]) = µ[j] + d[XY] * I(k = Y)
#    d[BC] = d[AC] - d[AB]
#    d[BD] = d[AD] - d[AB]
#    d[CD] = d[AD] - d[AC]

#    µ[j], d[AB], d[AC], d[AD] ~ Normal(0, 100^2)

# Note that in a trial comparing treatments X and Y, the model statement sets the log-odds 
# of an outcome in trial j on treatment k equal to trial ‘baseline’ µ[j] when k = X, and 
# µ[j] + d[XY] when k = Y. The trial baselines are all given unrelated, vague priors. The 
# (basic) parameters for treatment effects relative to the reference treatment A are also 
# given vague priors, while the remaining (functional) parameters are defined in terms of
# the basic parameters.

# Example 9.1 Smoking cessation MCT ---------------------------------------
# An example of a set of 24 trials comparing four interventions for smoking cessation: no 
# contact (A), self-help (B), individual counselling (C) and group counselling (D). There 
# is direct evidence on every one of the possible six pairwise contrasts, although the 
# majority of the data compare no contact (A) with individual counselling (C).

# Coding for pairwise meta-analysis introduced in Chapter 4 had separate likelihood 
# statements for treatment and control arms. This would become clumsy with three or more 
# treatments, so we adopt the approach taken in Chapter 6, and number the treatments 1, 2, 
# 3,... and so on. d[AB],d[AC],d[AD] will now be labelled d[2], d[3], d[4]. d[1] 
# corresponds to d[AA], the effect of treatment A relative to itself, and this is set to 
# zero.
library(R2jags)
library(rjags)
library(bayesplot)
library(tidyverse)


# Model data: 
data <- tribble(
~s, ~t,  ~r,    ~n,  ~b,
1,    1,    9,     140,   1,
1,    3,    23,    140,   1,
1,    4,    10,    138,   1,
2,    2,    11,    78,    2,
2,    3,    12,    85,    2,
2,    4,    29,    170,   2,
3,    1,    75,    731,   1,
3,    3,    363,   714,   1,
4,    1,    2,     106,   1,
4,    3,    9,     205,   1,
5,    1,    58,    549,   1,
5,    3,    237,   1561,  1,
6,    1,    0,     33,    1,
6,    3,    9,     48,    1,
7,    1,    3,     100,   1,
7,    3,    31,    98,    1,
8,    1,    1,     31,    1,
8,    3,    26,    95,    1,
9,    1,    6,     39,    1,
9,    3,    17,    77,    1,
10,   1,    79,    702,   1,
10,   2,    77,    694,   1,
11,   1,    18,    671,   1,
11,   2,    21,    535,   1,
12,   1,    64,    642,   1,
12,   3,    107,   761,   1, 
13,   1,    5,     62,    1,
13,   3,    8,     90,    1,
14,   1,    20,    234,   1,
14,   3,    34,    237,   1, 
15,   1,    0,     20,    1,
15,   4,    9,     20,    1,
16,   1,    8,     116,   1,
16,   2,    19,    149,   1, 
17,   1,    95,    1107,  1,
17,   3,    143,   1031,  1,
18,   1,    15,    187,   1,
18,   3,    36,    504,   1,
19,   1,    78,    584,   1, 
19,   3,    73,    675,   1, 
20,   1,    69,    1177,  1, 
20,   3,    54,    888,   1,  
21,   2,    20,    49,    2,
21,   3,    16,    43,    2,
22,   2,    7,     66,    2,
22,   4,    32,    127,   2,
23,   3,    12,    76,    3,
23,   4,    20,    74,    3,
24,   3,    9,     55,    3,
24,   4,    3,     26,    3
)
data <- as.matrix(data)
data

# The following codes the functional parameters in an efficient and compact way using 
# indexing. Thus, avoiding the need to actually having to specify the functional parameters
# explicitly at all!
model_String <- "model{
  # Likelihood:
   for (i in 1:50) {
    r[i] ~ dbin(p[i], n[i])
      
  # Sampling model:
   logit(p[i]) <- mu[s[i]] + d[t[i]] - d[b[i]]
   }
    
 # Priors on baselines:
 for (j in 1:24) {
  mu[j] ~ dnorm(0, .0001)
 }
 # Priors on treatment
 # effects:
 for (k in 2:4) {
  d[k] ~ dnorm(0, .0001)
 }
 
  # Set d[AA] to 0:
 d[1] <- 0
 
}
"
writeLines(text = model_String, con = "FixedMCT.txt")

data_JAGS <- list(s = data[, 1], t = data[, 2], r = data[, 3], n = data[, 4], 
                  b = data[, 5])
data_JAGS

#initial 1
inits <- list(
 list(
  d = c(NA, 0, 0, 0), mu = c(rep(0, 24))),
 list(
  d = c(NA, .1, -1, -.2), mu = c(1, -1, -2, 0,
                                 0, -2, 1, 0, 
                                 2, 2, 1, -1, 
                                 -2, 0, 0, -2, 
                                 1, 0, 2, 2,
                                 -2, -.5, -3, .5)
 ))

params <- c("d")

jags_Mod <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "FixedMCT.txt", inits = inits,
                 n.chains = 2, n.iter = 20000, n.burnin = 10000)
jags_Mod

# Visual Inspection of posterior:
posterior <- as.array(jags_Mod$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("d[2]", "d[3]", "d[4]"), 
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[,, 1:5], window = c(100, 150), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("d[2]", "d[3]", "d[4]"))

color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("d[2]", "d[3]", "d[4]"), 
           off_diag_args = list(size = 1.5))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("d[2]", "d[3]", "d[4]"), lags = 50)

# Absolute treatment effects ----------------------------------------------
# Cost-effectiveness analysis requires information on the absolute treatment difference on
# the probability scale. We can generate this from the LORs, but only if we have 
# information on the ‘baseline’ probability of the outcome.  In the context of Example 9.1, 
# the question is: what is the probability of smoking cessation in the ‘no treatment’ 
# group? This information could be based on one or more cohort studies that are considered
# to be representative of the target population, or on the more contemporary trials, or a 
# combination of both. The important point is that the issue of a suitable baseline for 
# the treatment A (no contact) strategy should be kept as a separate issue from the 
# calculation of relative treatment effects, which should be based on RCT evidence using 
# methods that respect the randomisation. For the sake of simplicity we assume that a 
# separate analysis has been conducted, which delivers a posterior distribution for the 
# log odds of smoking cessation under treatment A, Normal(-2.6, .38^2), which corresponds 
# to a median estimate of 7.5% with a 95% credible interval (3.4 to 14). We can then 
# construct absolute effects for the other treatments as follows, adding the relative 
# treatment effect to this baseline on the log-odds scale, then converting back to the 
# probability scale.

# A ~ Normal(-2.6, .38^2)
# logit(T[k]) = A + d[k]

model_String <- "model{
  # Likelihood:
   for (i in 1:50) {
    r[i] ~ dbin(p[i], n[i])
      
  # Sampling model:
   logit(p[i]) <- mu[s[i]] + d[t[i]] - d[b[i]]
   }
   
  # Absolute treatment
  # effect sampling 
  # model:
  for (k in 1:4) {
  logit(T[k]) <- A + d[k]
  }
 
 # Set d[1] to 0:
 d[1] <- 0
 
 # Vague Priors on baseline:
 for (j in 1:24) {
  mu[j] ~ dnorm(0, .0001)
 }
 # Vague Priors on treatment
 # effects:
 for (k in 2:4) {
  d[k] ~ dnorm(0, .0001)
 }
 # Prior Absolute treatment
 # effects:
 A ~ dnorm(-2.6, precA)
 
 # Transformation of
 # Absolute treatment
 # var to prec:
 precA <- pow(.38, -2)
 
}
"
writeLines(text = model_String, con = "AbsoluteMCT.txt")

params <- c("d", "T")

jags_Mod <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "AbsoluteMCT.txt", inits = inits,
                 n.chains = 2, n.iter = 20000, n.burnin = 10000)
jags_Mod
# ...really bad deviance for these mdoels. I assume this will be addressed later in the 
# chapter?

# Visual Inspection of posterior:
posterior <- as.array(jags_Mod$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("T[1]", "T[2]", "T[3]", "T[4]"),
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[,, 1:5], window = c(100, 150), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("T[1]", "T[2]", "T[3]", "T[4]"))

color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("T[2]", "d[2]"),
           off_diag_args = list(size = 1.5))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("d[2]", "d[3]", "d[4]"), lags = 50)

mcmc_hist(x = posterior, c("T[2]", "T[3]"))

# Relative treatment efficacy and ranking ---------------------------------
# The code presented so far delivers posterior LORs for the efficacy of treatments B, C, D
# relative to A, but it does not address the question of how to make inferences about 
# treatment efficacy. We can clearly generate the LORs for any treatment relative to any 
# other, and also monitor the Odds Ratios:
model_String <- "model{
  # Likelihood:
   for (i in 1:50) {
    r[i] ~ dbin(p[i], n[i])
      
  # Sampling model:
   logit(p[i]) <- mu[s[i]] + d[t[i]] - d[b[i]]
   }
   
  # Absolute treatment
  # effect sampling 
  # model:
  for (k in 1:4) {
  logit(T[k]) <- A + d[k]
  
  }
 # Set d[1] to 0:
 d[1] <- 0
 
 # Vague Priors on baseline:
 for (j in 1:24) {
  mu[j] ~ dnorm(0, .0001)
 }
 # Vague Priors on treatment
 # effects:
 for (k in 2:4) {
  d[k] ~ dnorm(0, .0001)
 }
 # Prior Absolute treatment
 # effects:
 A ~ dnorm(-2.6, precA)
 
 # Transformation of
 # Absolute treatment
 # var to prec:
 precA <- pow(.38, -2)
 
 for (c in 1:3) { 
 # All pair-wise comparison log odds ratios:
  for (k in (c + 1):4) { 
  # and single study odds ratios:
  logOR[c, k] <- d[k] - d[c] 
  OR[c, k] <- exp(logOR[c, k])
  }
 }
 
}
"
writeLines(text = model_String, con = "relativeTreat_MCT.txt")

params <- c("d", "T", "OR")

jags_Mod <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "relativeTreat_MCT.txt", inits = inits,
                 n.chains = 2, n.iter = 20000, n.burnin = 10000)
jags_Mod

# Visual Inspection of posterior:
posterior <- as.array(jags_Mod$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("OR[1,2]", "OR[1,3]", "OR[2,3]", "OR[1,4]",
                               "OR[2,4]", "OR[3,4]"),
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[,, 1:5], window = c(100, 150), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("OR[1,2]", "OR[1,3]", "OR[2,3]", "OR[1,4]",
                               "OR[2,4]", "OR[3,4]")
                  )

color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("OR[1,2]", "OR[1,3]"),
           off_diag_args = list(size = 1.5))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("OR[1,2]", "OR[1,3]", "OR[2,3]", "OR[1,4]",
                               "OR[2,4]", "OR[3,4]"),
         lags = 50)

color_scheme_set("mix-teal-pink")
mcmc_hist(x = posterior, c("OR[1,2]", "OR[1,3]", "OR[2,3]", "OR[1,4]",
                               "OR[2,4]", "OR[3,4]"))

# But we are then confronted with a typical multiple comparisons problem. The approach we 
# suggest is to rank the treatments and examine the posterior distributions of the ranks, 
# and also to calculate the probability that each treatment is the best treatment. In the 
# following code the rank(v,s) function returns the number of elements of the vector v
# whose value is less than or equal to the sth element. Then best[k] takes the value 1 
# when treatment k has the highest cessation rate and 0 otherwise. We also generate the 
# LORs and Odds Ratios for all six contrasts.
model_String <- "model{
  # Likelihood:
   for (i in 1:50) {
    r[i] ~ dbin(p[i], n[i])
      
  # Sampling model:
   logit(p[i]) <- mu[s[i]] + d[t[i]] - d[b[i]]
   }
   
  # Absolute treatment
  # effect sampling 
  # model:
  for (k in 1:4) {
  logit(T[k]) <- A + d[k]
  
  }
 # Set d[1] to 0:
 d[1] <- 0
 
 # Vague Priors on baseline:
 for (j in 1:24) {
  mu[j] ~ dnorm(0, .0001)
 }
 # Vague Priors on treatment
 # effects:
 for (k in 2:4) {
  d[k] ~ dnorm(0, .0001)
 }
 # Prior Absolute treatment
 # effects:
 A ~ dnorm(-2.6, precA)
 
 # Transformation of
 # Absolute treatment
 # var to prec:
 precA <- pow(.38, -2)
 
 # Log-Odds calculations for
 # each comparison
 for (c in 1:3) { 
 # All pair-wise comparison 
 # log odds ratios:
  for (k in (c + 1):4) { 
  
  # and single study comparison 
  # odds ratios:
  OR[c, k] <- d[k] - d[c] 
  log(LOR[c, k]) <- OR[c,k]
  }
 }
 # Rank treatment effect
 # (where 1 = best):
 rk <- 5 - rank(T[])
 # & record the best treatment:
 best <- equals(rk, 1)
}
"
writeLines(text = model_String, con = "rankCessationMCT.txt")

params <- c("d", "T", "LOR", "best", "rk", "best")

jags_Mod <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "rankCessationMCT.txt", inits = inits,
                 n.chains = 2, n.iter = 40000, n.burnin = 20000)
jags_Mod

# Visual Inspection of posterior:
posterior <- as.array(jags_Mod$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("mix-teal-pink")
mcmc_hist(x = posterior, c("rk[1]", "rk[2]", "rk[3]", "rk[4]"),
          binwidth = .5)

# The information in the rankings is neatly summarised in the vector "best". The posterior 
# means probabilities are : best[1] = 0, best[2] = 0, best[3] = 0.33, best[4] = 0.67. 
# These results appear to rule out no contact and self-help. Group counselling has the 
# highest probability of being the ‘best’.

# ========================================================================================
# Random Effects MTC Models -----------------------------------------------
# ========================================================================================
# Both the model and the code introduced below are the natural extension of the Random 
# Effects model for pairwise comparisons, discussed in Chapter 4, and the Fixed Effect MTC
# model of the previous section.

# In the Random Effects model each trial j on treatment contrast XY estimates a distinct 
# LOR, ∂[jXY], which is drawn from a common distribution ∂[jXY] ~ N(d[XY], sigma^2). We 
# will make the simplyfying assumption that the between-trial variance for all six 
# contrasts are equal, such that sigma^2[XY] = sigma^2. Adding a vague Uniform prior for
# sigma the full model becomes:

# r[jk] ~ dbin(p[jk], n[jk])
# logit(p[jk]) = µ[j] + ∂[jXY] * I(k = Y)
# ∂[jXY] ~ N(d[XY], sigma^2)
# 
# d[BC] = d[AC] - d[AB]
# d[BD] = d[AD] - d[AB]
# d[CD] = d[AD] - d[AC]
# 
# µ[j], d[AB], d[AC], d[AD] ~ N(0, 100^2)
# sigma ~ dunif(0, 2)

# Example 9.1 revisited ---------------------------------------------------
# There are a number of ways to code a Random Effects model. One approach is to modify the
# code for the Fixed Effect version, as follows:
model_String <- "model{
  # Likelihood:
   for (i in 1:50) {
    r[i] ~ dbin(p[i], n[i])
      
  # Sampling model:
   logit(p[i]) <- mu[s[i]] + delta[i] * (1 - equals(t[i], b[i]))
   # Random effects distribution:
    delta[i] ~ dnorm(md[i], prec)
    
   # Mean of random 
   # effect distribution:
   md[i] <- d[t[i]] - d[b[i]]

   }
   
  # Absolute treatment
  # effect sampling 
  # model:
  for (k in 1:4) {
  logit(T[k]) <- A + d[k]
  
  }
 # Set d[1] to 0:
 d[1] <- 0
 
 # Vague Priors on baseline:
 for (j in 1:24) {
  mu[j] ~ dnorm(0, .0001)
 }
 # Vague Priors on treatment
 # effects:
 for (k in 2:4) {
  d[k] ~ dnorm(0, .0001)
 }
 # Prior Absolute treatment
 # effects:
 A ~ dnorm(-2.6, precA)
 # Transformation of
 # Absolute treatment
 # var to prec:
 precA <- pow(.38, -2)
 
 # Vague prior on RE sd:
 sd ~ dunif(0, 2)
 # Var function of RE 
 # sd:
 tau.sq <- sd * sd
 # RE precision:
 prec <- 1 / tau.sq
 
 # Log-Odds calculations for
 # each comparison
 for (c in 1:3) { 
 # All pair-wise comparison 
 # log odds ratios:
  for (k in (c + 1):4) { 
  
  # and single study comparison 
  # odds ratios:
  OR[c, k] <- d[k] - d[c] 
  log(LOR[c, k]) <- OR[c,k]
  }
 }
 # Rank treatment effect
 # (where 1 = best):
 rk <- 5 - rank(T[])
 # & record the best treatment:
 best <- equals(rk, 1)
 
}
"
writeLines(text = model_String, con = "rEffectsMCT.txt")

params <- c("d", "T", "LOR", "best", "sd")

jags_Mod <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "rEffectsMCT.txt", inits = inits,
                 n.chains = 2, n.iter = 40000, n.burnin = 20000)
jags_Mod

# Visual Inspection of posterior:
posterior <- as.array(jags_Mod$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("mix-teal-pink")
mcmc_hist(x = posterior, c("best[2]", "best[3]", "best[4]"),
          binwidth = .6)

# ========================================================================================
# Model choice and consistency of MTC evidence ----------------------------
# ========================================================================================
# An important finding from the random effect analysis concerns the sigma parameter. Not 
# only is its mean value of the same order as the mean treatment effects, but the lower 
# credible limit, 0.54, is so high as to effectively rule out the hypothesis that sigma is
# close to zero. This points us firmly in the direction of the Random Effects model. This 
# can be put on a slightly more formal basis by comparing the Fixed and Random Effects 
# models using some of the model critique methods from Chapter 4.

# Random effects model:
model_String <- "model{
  # Likelihood:
   for (i in 1:50) {
    r[i] ~ dbin(p[i], n[i])
      
  # Sampling model:
   logit(p[i]) <- mu[s[i]] + delta[i] * (1 - equals(t[i], b[i]))
   # Random effects distribution:
    delta[i] ~ dnorm(md[i], prec)
    
   # Mean of random 
   # effect distribution:
   md[i] <- d[t[i]] - d[b[i]]

   }
   
  # Absolute treatment
  # effect sampling 
  # model:
  for (k in 1:4) {
  logit(T[k]) <- A + d[k]
  
  }
 # Set d[1] to 0:
 d[1] <- 0
 
 # Vague Priors on baseline:
 for (j in 1:24) {
  mu[j] ~ dnorm(0, .0001)
 }
 
 # Vague Priors on treatment
 # effects:
 for (k in 2:4) {
  d[k] ~ dnorm(0, .0001)
 }
 
 # Prior Absolute treatment
 # effects:
 A ~ dnorm(-2.6, precA)
 # Transformation of
 # Absolute treatment
 # var to prec:
 precA <- pow(.38, -2)
 
 # Vague prior on RE sd:
 sd ~ dunif(0, 2)
 # Var function of RE 
 # sd:
 tau.sq <- sd * sd
 # RE precision:
 prec <- 1 / tau.sq
 
 # Log-Odds calculations for
 # each comparison
 for (c in 1:3) { 
 # All pair-wise comparison 
 # log odds ratios:
  for (k in (c + 1):4) { 
  
  # and single study comparison 
  # odds ratios:
  OR[c, k] <- d[k] - d[c] 
  log(LOR[c, k]) <- OR[c,k]
  }
 }
 # Rank treatment effect
 # (where 1 = best):
 rk <- 5 - rank(T[])
 # & record the best treatment:
 best <- equals(rk, 1)
 
    # Model deviance 
    # calculations:
     for (i in 1:50) {
 
    # Predicted model 
    # deviance:
    rhat[i] <- p[i] * n[i]
    # Deviance of each 
    # data point:
    dev[i] <- 2 * (r[i] * (log(r[i]) - log(rhat[i])) + 
    (n[i] - r[i]) * (log(n[i] - r[i]) - log(n[i] - rhat[i])))
     }
     
     # Residual deviance:
     resdev <- sum(dev[])
}
"

writeLines(text = model_String, con = "RandEff_dev_MCT.txt")

params <- c("resdev")

jags_Mod <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "RandEff_dev_MCT.txt", inits = inits,
                 n.chains = 2, n.iter = 40000, n.burnin = 20000)
jags_Mod

# ... versus fixed effects model:
model_String <- "model{
  # Likelihood:
   for (i in 1:50) {
    r[i] ~ dbin(p[i], n[i])
      
  # Sampling model:
   logit(p[i]) <- mu[s[i]] + d[t[i]] - d[b[i]]
   }
   
  # Absolute treatment
  # effect sampling 
  # model:
  for (k in 1:4) {
  logit(T[k]) <- A + d[k]
  
  }
 # Set d[1] to 0:
 d[1] <- 0
 
 # Vague Priors on baseline:
 for (j in 1:24) {
  mu[j] ~ dnorm(0, .0001)
 }
 
 # Vague Priors on treatment
 # effects:
 for (k in 2:4) {
  d[k] ~ dnorm(0, .0001)
 }
 
 # Prior Absolute treatment
 # effects:
 A ~ dnorm(-2.6, precA)
 # Transformation of
 # Absolute treatment
 # var to prec:
 precA <- pow(.38, -2)
 
 # Log-Odds calculations for
 # each comparison
 for (c in 1:3) { 
 # All pair-wise comparison 
 # log odds ratios:
  for (k in (c + 1):4) { 
  
  # and single study comparison 
  # odds ratios:
  OR[c, k] <- d[k] - d[c] 
  log(LOR[c, k]) <- OR[c,k]
  }
 }
 
 # Rank treatment effect
 # (where 1 = best):
 rk <- 5 - rank(T[])
 # & record the best treatment:
 best <- equals(rk, 1)
 
    # Model deviance 
    # calculations:
    for (i in 1:50) {
 
    # Predicted model 
    # deviance:
    rhat[i] <- p[i] * n[i]
    # Deviance of each 
    # data point:
    dev[i] <- 2 * (r[i] * (log(r[i]) - log(rhat[i])) + 
    (n[i] - r[i]) * (log(n[i] - r[i]) - log(n[i] - rhat[i])))
    }
    
    # Residual deviance:
    resdev <- sum(dev[])

}
"
writeLines(text = model_String, con = "FixEff_dev_MCT.txt")

params <- c("resdev")

jags_Mod <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "FixEff_dev_MCT.txt", inits = inits,
                 n.chains = 2, n.iter = 40000, n.burnin = 20000)
jags_Mod

# Deviance calculated in for the RE model is close to the number of observations (50), and
# we would be justified in concluding that the Random Effects model provides an adequate 
# fit to the data. Of course, a Random Effects model is extremely tolerant. The variance 
# term will happily stretch to fit trials whose values are far from the mean without 
# producing any sign that the model fit is poor. As we saw with the magnesium 
# meta-analysis, a globally poor fit can only be obtained if one or two very large trials 
# are distinctly far from the mean of the others. Therefore, the comparison of Fixed and 
# Random Effects models tells us mostly about the level of between-trial heterogeneity 
# within the different comparison types. It may not tell us much about whether the key 
# consistency assumptions are being met.

# End file ----------------------------------------------------------------