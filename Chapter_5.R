# Chapter concerned with meta-regression analysis, and the use of covaraities in evidence
# synthesis
library(R2jags)
library(rjags)
library(bayesplot)
library(ggplot2)
library(magrittr)
library(parallel)
library(tidyverse)

options(mc.cores = detectCores())

# Remember: BUGS/JAGS requires five key elements: specification of the likelihood, a 
# description of the statistical model linking parameters to the likelihood, 
# specification of prior distributions for the parameters, the data, and initial values 
# for the sampler.

# "Meta-regression relates the size of the treatment effect to numerical characteristic(s) 
# of the trials included in a meta-analysis. The relationship modelled is like a standard 
# (e.g. least-squares) regression, only the precision of each study’s outcome estimate is 
# taken into account (just like in standard meta-analysis models). Either fixed or random 
# effect meta-analysis models can be extended to include covariates, however only the 
# random effects version is pursued here (although the fixed effect version is a 
# straightforward simplification of the model presented). Generally, ran- dom effects 
# meta-regression models are preferred since any residual heterogeneity, not explained by 
# the included covariates, is allowed for via the random effect term."

# ================================================================================
# Generic Random-Effect Meta-Regression Model -----------------------------
# ================================================================================

# Section 5.2.2 -----------------------------------------------------------
# Random effects meta-regression model for Odds Ratio (OR) outcomes using a Binomial 
# likelihood, i.e. a logistic regression, as it is a special case of the binomial model.
# Data based on trials of the BCG vaccine to prevent TB. The model is given by:

# r[Ai] ~ Binomial(p[Ai], n[Ai])         r[Bi] ~ Binomial(p[Bi], n[Bi])
# logit(p[Ai]) = µ[i]                    logit(p[Bi]) = µ[i] + ∂[i] + ßx[i]
# ∆[i] ~ Normal(d, tau^2)                i = 1, ..., k

# ß * x[i] (the regression coefficient multiplied by the covariate value for the i'th 
# study) has been added to the linear predictor for the effect in the treatment group:
Model_string <- "model {
   # Binomial Likelihood 
   # specification of 
   # i'th study:
   for (i in 1:Nstud) {
    rA[i] ~ dbin(pA[i], nA[i])
    rB[i] ~ dbin(pB[i], nB[i])
  
    # Logistic regression
    # model, for i'th study,
    # linking to Binomial 
    # statistic p:
     logit(pA[i]) <- mu[i]
     logit(pB[i]) <- mu[i] + delta[i] + beta * lat[i]
  
   # Statistical model, 
   # linking the
   # parameters to likelihood:
    delta[i] ~ dnorm(d, prec)
  
  }
 # Prior on µ[i], the prior 
 # log(odds) of an event for 
 # i'th study:
  for (i in 1:Nstud) {
   mu[i] ~ dnorm(0.0, 1.0e-5)
  
  }
 # Prior on intercept (d), 
 # the pooled treatment 
 # effect where 
 # covariate lat = 0:
  d ~ dnorm(0.0, 1.0e-6)
 # Data transformation for d:
  OR <- exp(d)
  
 # Prior on tau, sd for
 # between studies:
  tau ~ dunif(0, 10)
 # tau.sq, between-study 
 # variance:
  tau.sq <- tau * tau
 # Precision as a 
 # fucntion of tau:
  prec <- 1 / (tau.sq)
  
 # Prior on Beta, the
 # coefficient:
  beta ~ dnorm(0.0, 1.0e-6)
}
"
writeLines(text = Model_string, con = "Model.5.2.2.txt")

# Treatment vaccinated = B
# Treatment non-vaccinated = A
Data_jags <- list(
 rB = c(4, 6, 3, 62, 33, 180, 8, 505, 29, 17, 186, 5, 27),
 rA = c(11, 29, 11, 248, 47, 372, 10, 499, 45, 65, 141, 3, 29),
 nB = c(123, 306, 231, 13598, 5069, 1541, 2545, 88391, 7499, 1716, 50634, 2498, 16913),
 nA = c(139, 303, 220, 12867, 5808, 1451, 629, 88391, 7277, 1665, 27338, 2342, 17854),
 lat = c(44, 55, 42, 52, 13, 44, 19, 13, -27, 42, 18, 33, 33),
 Nstud = 13
)

inits <- list(d = 0, tau = 1, delta = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
 mu = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), beta = 0
)

jags_Model <- jags(data = Data_jags, model.file = "Model.5.2.2.txt", 
                   parameters.to.save = c("beta", "d", "tau", "tau.sq", "OR"),
                  n.iter = 20000, n.burnin = 1000)
jags_Model

posterior <- as.array(jags_Model$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("d", "beta", "tau"), 
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[,, 1:5], window = c(100, 150), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("d", "beta"))

color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("d", "beta"),
           off_diag_args = list(size = 1.5))
# ... seems there to be a high negative autocorrelation between the beta and d variables.
color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("d", "beta"), lags = 50)

# When such plots are obtained, the possibility of lack of convergence should always be 
# considered. In this instance, as will be explained, the issue here is a poor mixing 
# MCMC chain rather than lack of convergence. In such instances, the MCMC sampler will 
# give valid parameter estimates (as indeed are those presented above), if the model is 
# run for enough iterations, but since successive iterations, even at a lag of 50, are 
# still highly correlated, the sampler is very inefficient.

# Such inefficiencies could be much more limiting in more complex models so an 
# explanation of what causes the problem together with a solution are presented below.

 # A scatter plot of pairs of estimates for the model parameters d and ß as sampled by 
# the model are presented. The shape formed by this scatter of points indicates that 
# there is a strong negative association between associated estimates sampled for these 
# parameters (i.e. sampling a high value of ß results in a low value of d and vice 
# versa). It is this correlation between intercept and slope parameter of the regression
# which is responsible for the poor mixing of the MCMC chains.

# Fortunately, there is a relatively simple solution to this problem, which is achieved 
# by centring the covariates. This means taking the mean covariate value away from each 
# covariate value, which has the effect of moving  the origin of the regression to the 
# mean covariate value in the centre of the data. This a must practice for Bayesian
# analysis in regression models! See McElreath 'Statistical Rethinking'.

# Centring the covariates can be done outside BUGS/JAGS before analysis, or actually 
# within the analysis code by replacing the line of code:

# logit(pB[i]) <- mu[i] + delta[i] + beta * lat[i]

# with...

# logit(pB[i]) <- mu[i] + delta[i] + beta * (lat[i] - mean(lat[]))

Model_string <- "model {
   # Binomial Likelihood for the 
   # i'th study:
   for (i in 1:Nstud) {
    rA[i] ~ dbin(pA[i], nA[i])
    rB[i] ~ dbin(pB[i], nB[i])
  
    # Logistic regression for
    # the i'th study:
     logit(pA[i]) <- mu[i]
     logit(pB[i]) <- mu[i] + delta[i] + beta * (lat[i] - mean(lat[]))
  
   # Statistical model, 
   # linking the
   # parameters to likelihood:
    delta[i] ~ dnorm(d, prec)
  
  }
 # Prior on µ, the estimate 
 # log(odds) of an event for 
 # i'th study:
  for (i in 1:Nstud) {
   mu[i] ~ dnorm(0.0, 1.0e-5)
  
  }
 # Prior on intercept (d), 
 # the mean pooled treatment 
 # effect:
  d ~ dnorm(0.0, 1.0e-6)
 # Data transformation for d:
  OR <- exp(d)
  
 # Prior on tau, sd for
 # between studies:
  tau ~ dunif(0, 10)
 # tau.sq, between-study 
 # variance:
  tau.sq <- tau * tau
 # Precision as a 
 # fucntion of tau:
  prec <- 1 / (tau.sq)
  
 # Prior on Beta, the
 # coefficient:
  beta ~ dnorm(0.0, 1.0e-6)
}
"
writeLines(text = Model_string, con = "Model.5.2.2.txt")

Data_jags <- list(
 rB = c(4, 6, 3, 62, 33, 180, 8, 505, 29, 17, 186, 5, 27),
 rA = c(11, 29, 11, 248, 47, 372, 10, 499, 45, 65, 141, 3, 29),
 nB = c(123, 306, 231, 13598, 5069, 1541, 2545, 88391, 7499, 1716, 50634, 2498, 16913),
 nA = c(139, 303, 220, 12867, 5808, 1451, 629, 88391, 7277, 1665, 27338, 2342, 17854),
 lat = c(44, 55, 42, 52, 13, 44, 19, 13, -27, 42, 18, 33, 33),
 Nstud = 13
)

inits <- list(d = 0, tau = 1, delta = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
 mu = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), beta = 0
)

jags_Model <- jags(data = Data_jags, model.file = "Model.5.2.2.txt", 
                   parameters.to.save = c("beta", "d", "tau", "tau.sq", "OR", "rA",
                                          "rB", "nA", "nB"))
jags_Model

posterior <- as.array(jags_Model$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("d", "beta", "tau"), 
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("OR", "d", "beta"), 
           window = c(100, 150), size = 1, ) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("d", "beta"))

color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("d", "beta"),
           off_diag_args = list(size = 1.5))
# ... the above indicates the autocorrelation has been *greatly* reduced!
color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("d", "beta"), lags = 10)

# Note: once the model is fitted with this new line of code and parameter estimates 
# obtained, it is important to interpret the results correctly and, if necessary, 
# transform the parameters to remove the centring (i.e. the model should produce exactly
# the same results but the latter is a different parameterisation of the first and the 
# parameters need transforming to give the same interpretation as the uncentred model).

# As fitted, the interpretation of the parameter d is the treatment effect at the mean 
# covariate value, where before it was the effect when the covariate value was 0. To 
# obtain the latter, beta x mean(lat[]) is subtracted from the intercept (d) to remove 
# centring and this can be done within BUGS/JAGS through the creation of a new node, 
# i.e:

#  d.uncent <- d - beta*mean(lat[ ])

Model_string <- "model {
   # Binomial Likelihood for the 
   # i'th study:
   for (i in 1:Nstud) {
    rA[i] ~ dbin(pA[i], nA[i])
    rB[i] ~ dbin(pB[i], nB[i])
  
   # Logistic regression for
   # the i'th study:
    logit(pA[i]) <- mu[i]
    logit(pB[i]) <- mu[i] + delta[i] + beta * (lat[i] - mean(lat[]))
  
   # Statistical model, 
   # linking the
   # parameters to likelihood:
    delta[i] ~ dnorm(d, prec)
  
  }
 # Prior on µ, the estimate 
 # log(odds) of an event for 
 # i'th study:
  for (i in 1:Nstud) {
   mu[i] ~ dnorm(0.0, 1.0e-5)
  
  }
 # Prior on intercept (d), 
 # the mean pooled treatment 
 # effect:
  d ~ dnorm(0.0, 1.0e-6)
 # Uncentering d:
  d.uncent <- d - beta * mean(lat[])
 # Data transformation for d:
  OR <- exp(d)
  
 # Prior on tau, sd for
 # between studies:
  tau ~ dunif(0, 10)
 # tau.sq, between-study 
 # variance:
  tau.sq <- tau * tau
 # Precision as a 
 # fucntion of tau:
  prec <- 1 / (tau.sq)
  
 # Prior on Beta, the
 # coefficient:
  beta ~ dnorm(0.0, 1.0e-6)
}
"
writeLines(text = Model_string, con = "Model.5.2.2.txt")

Data_jags <- list(
 rB = c(4, 6, 3, 62, 33, 180, 8, 505, 29, 17, 186, 5, 27),
 rA = c(11, 29, 11, 248, 47, 372, 10, 499, 45, 65, 141, 3, 29),
 nB = c(123, 306, 231, 13598, 5069, 1541, 2545, 88391, 7499, 1716, 50634, 2498, 16913),
 nA = c(139, 303, 220, 12867, 5808, 1451, 629, 88391, 7277, 1665, 27338, 2342, 17854),
 lat = c(44, 55, 42, 52, 13, 44, 19, 13, -27, 42, 18, 33, 33),
 Nstud = 13
)

inits <- list(d = 0, tau = 1, delta = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
 mu = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), beta = 0
)

jags_Model <- jags(data = Data_jags, model.file = "Model.5.2.2.txt", 
                   parameters.to.save = c("beta", "d", "d.uncent", "tau", "tau.sq", 
                                          "OR", "prec"))
jags_Model
#  It can be seen that a strong negative association is estimated (ß is estimared to be
# 0.02 (95% CrI -0.03 to -0.01)), which implies that the Log Odds Ratio (LOR) decreases 
# (linearly) the further away from the equator you are – implying the treatment increases
# in effectiveness. Hence, this means the meta-analysis has not produced a single pooled 
# effectiveness estimate, but implies effectiveness varies with location on the globe.

# The estimate for the between study heterogeneity, tau.sq, is much smaller than was 
# estimated before the inclusion of the covariate, which confirms that much of the 
# between study heterogeneity – but not all – has been explained by the inclusion of the 
# covariate.

posterior <- as.array(jags_Model$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("d", "d.uncent", "beta", "tau"), 
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("OR", "d", "d.uncent", "beta"), 
           window = c(100, 150), size = 1, ) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("d", "beta"))

color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("d", "beta"),
           off_diag_args = list(size = 1.5))
# ... the above indicates the autocorrelation has been *greatly* reduced!
color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("d", "beta"), lags = 10)

# ====================================================================================
# Limitations of meta-regression ------------------------------------------
# ====================================================================================
# However, this example is not typical, and in many instances exploration of 
# heterogeneity via covariates will be much less conclusive. The typical numbers of 
# studies in a meta-analysis severely limits the power such analyses have in many 
# contexts.  Where analyses are based on aggregate patient level covariates (e.g.
# % white, average age, etc.), typically, such analyses will have much less power than 
# if the analysis was carried out at the patient (rather than the study) level through 
# the use of Individual Patient Data (IPD). It is recommended where feasible, and 
# exploring the effect of patient level covariates is a context where IPD will be 
# particularly valuable. 

# A further issue in fitting regression models is that often there is concern that the 
# severity of illness interacts with the effectiveness of the intervention. For this 
# reason the effect of ‘baseline risk’ on outcome is commonly explored through 
# meta-regression. However, a complication with this analysis is that the covariate of 
# interest is actually also part of the outcome definition. A specific model which allows
# for this is considered in the next.

# Baseline Risk -----------------------------------------------------------
# Heterogeneity in baseline risk among trials is likely to reflect differences in patient
# characteristics (e.g. age, medical history, co-morbidities, etc.) and thus baseline 
# risk may be considered a proxy, or surrogate, for such underlying factors; data for 
# which may not be available to the analyst.

# However, if we formally define baseline risk to be the risk of the outcome event for a 
# patient under the control condition, since this indicates the average risk of patient 
# in that trial if they were not treated, then the covariate information would be derived
# from the outcome data from the control group in each trial. Hence, the same data are 
# used to inform both outcome and covariate values (i.e. both y and x in a regression). 
# This leads to structural dependence within the regression equation. Additionally, both 
# the covariate and outcome are estimated from the trials with finite sample sizes and 
# therefore are estimates rather than true values (i.e. they are measured with error). 
# The issues of structural dependency and measurement error (when combined) present 
# problems associated with regression to the mean, also known as regression dilution 
# bias. If a standard meta-regression model of the form described above is fitted for
# baseline risk, which ignores regression dilution bias, then the association between 
# covariate and outcome can be overestimated.


# Example 5.2 -------------------------------------------------------------
# BCG TB endoscopic sclerotherapy in the prevention of bleeding for patients with 
# cirrhosis and oesagogastric varices.  The outcome of interest is the number of patients
# who go on to develop a bleed.

data_Subjects <- tribble(
~Control, ~Treatment, 
36, 35,
53, 56, 
18, 16, 
22, 23, 
46, 49, 
60, 53, 
60, 53, 
69, 71, 
41, 41, 
20, 21, 
41, 42, 
35, 33,
138, 143, 
51, 55, 
72, 73, 
16, 13, 
28, 21, 
19, 18, 
24, 22)

nA <- data_Subjects[, 1]
nA <- nA$Control
nA

nB <- data_Subjects[, 2]
nB <- nB$Treatment
nB

data_Bleeds <- tribble(
 ~Control, ~Treatment,
 22, 3, 
 30, 5, 
 6, 5, 
 9, 3, 
 31, 11, 
 9, 19, 
 26, 17, 
 29, 10, 
 14, 12, 
 3, 0, 
 13, 9, 
 14, 13, 
 23, 31, 
 19, 20, 
 13, 13, 
 12, 3, 
 5, 3, 
 0, 4, 
 2, 6
) 
data_Bleeds

rA <- data_Bleeds[, 1]
rA <- rA$Control
length(rA)
length(nA)

rB <- data_Bleeds[, 2]
rB <- rB$Treatment
length(rB)
length(nB)

Model_string <- "model {
 for (i in 1:19) {
  # Binomial likelihood
  # sampling model:
  rA[i] ~ dbin(pA[i], nA[i])
  rB[i] ~ dbin(pB[i], nB[i])
  
  # Logistic regression
  # model:
  logit(pA[i]) <- mu[i]
  logit(pB[i]) <- mu[i] + delta[i] + beta * (mu[i] - mean(mu[]))
  
  # Statistical sampling
  # model:
  delta[i] ~ dnorm(d, prec)
 }
 # Prior on µ:
 for (i in 1:19) {
  mu[i] ~ dnorm(0.0, 1.0e-5)
 }
 # Prior on intercept d, 
 # pooled effect estimate
 d ~ dnorm(0.0, 1.0e-6)
 
 # Prior on tau, 
 # between study 
 # sd:
 tau ~ dunif(0, 10)
 # Variance:
 tau.sq <- tau * tau
 # Precision as a
 # function of var:
 prec <- 1 / (tau.sq)
 
 # Prior on ß:
 beta ~ dnorm(0.0, 1.0e-6)
}
"
writeLines(text = Model_string, con = "Base_Risk.txt")

Data_jags <- list(
 rB = rB,
 rA = rA,
 nB = nB,
 nA = nA
)
Data_jags

inits <- list(
 d = 0, tau = 1, delta = c(rep(0, 19)), mu = c(rep(0, 19)), beta = 0
)
inits

jags_Model <- jags(data = Data_jags, model.file = "Base_Risk.txt", 
                   parameters.to.save = c("beta", "d", "tau.sq", "prec")
                   )
jags_Model

# covariate.

posterior <- as.array(jags_Model$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("d", "beta", "tau.sq"), 
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("d", "beta"), 
           window = c(100, 150), size = 1, ) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("d", "beta"))

color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("d", "beta"),
           off_diag_args = list(size = 1.5))
# ... the above indicates the autocorrelation has been *greatly* reduced!
color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("d", "beta"), lags = 10)

# Note: that no separate column of data is specified for this covariate, i.e. it is the 
# same µ[i] that has always been included in the modelling and is treated as a random 
# variable (since it is defined as the log odds of pA[i]), thus its value will vary at 
# each iteration of the MCMC simulation. By doing so, the structural dependence between
# the intercept and slope of this model are correctly accounted for (Note, in other 
# contexts, such a model, which allows for uncertainty in covariate values, is referred 
# to as a measurement error model.)

# Hence, as predicted the ‘na ̈ıve’ analysis overestimates the association, although the
# bias is modest in this example. In instances where there are small numbers of small 
# trials, the bias is potentially much greater.

# End file ----------------------------------------------------------------