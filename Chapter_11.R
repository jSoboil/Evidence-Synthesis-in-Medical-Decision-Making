library(R2jags)
library(rjags)
library(bayesplot)
library(ggplot2)
library(magrittr)
library(tidyverse)

# ==========================================================================================
# Generalised Evidence Synthesis ------------------------------------------
# ==========================================================================================
# This chapter focuses on integrating non-RCT evidence. However, it should be noted that,
# despite its pros, generalised evidence synthesis is still very much in an experimental stage.

# Example 11.1 ------------------------------------------------------------
# This example focuses on how to account for potential biases in the evidence chosen. This
# first approach considers a meta-analysis of the RCT evidence in which an informative prior
# distribution on the treatment effect is derived from the observational data. Ibrahim and 
# Chen have presented a framework in which a power transform of the data likelihood is 
# considered:

#  P(theta | Data) = L(theta | RCTs) * [L(theta | Obs)] ^ alpha * P(theta),

# where P(theta | RCTs) is the posterior distribution for a set of model parameters, 
# L(theta | RCTs) is the likelihood function for the RCT evidence, and L(theta | Obs) is the
# prior based on the observational evidence.

# Note how the likelihood function for the observational data is raised to the power alpha.
# If alpha takes a value between 0 and 1 it reduces the likelihood and thus discounts the 
# observational data and essentially removes it from the analysis. A value of 1 for alpha 
# implies the observational evidence is accepted at 'face value' and is not down-weighted at 
# all. Note that prior distributions, P(theta), for all parameters q are still required. The
# challenge in using this approach is deciding on the appropriate value for alpha to use. 
# Given this difficulty, one approach is to evaluate the meta-analysis for a range of values 
# of alpha, and explore how the estimate changes with alpha in a sensitivity analysis. This 
# approach gives the potential to identify threshold values for alpha for which decisions 
# change and decide how plausible such values are. An example of this approach, where a 
# meta-analysis of the observational evidence is conducted prior to the final meta-analysis
# is shown on page 231.

# ==========================================================================================
# Hierarchical Models for Evidence from different study designs -----------
# ==========================================================================================
# This section considers an extension of the random effects model presented in Chapter 4. Here
# we model an extra level of variation to allow for variability in effect sizes between 
# different sources of evidence (in addition to allowing for variability between study 
# estimates within each study type). The model can be used when there are three or more 
# different study types to include in the synthesis.

# As an illustrative example, consider (some of) the evidence for the effect of electronic 
# fetal heart rate monitoring (EFM) on perinatal mortality. EFM has not been shown to reduce
# perinatal mortality in the 9 RCTs available at the time of the initial analysis (around 
# 2000), although this may be due to low power since perinatal mortality is rare; thus the 
# potential opportunity to consider a wider evidence base. In the evidence considered here, 
# there are the aforementioned RCTs, comparative cohort studies and before and after studies.

# The analysis was conducted on the Risk Difference (RD) scale, although detailed 
# consideration of this data set elsewhere may now suggest that the Relative Risk may have 
# been a better choice due to issues with excessive heterogeneity on the RD scale. 

# If we assume that the studies of each particular type of are exchangeable and that the 
# pooled estimates obtained from each study type are themselves exchangeable, then this can
# be expressed formally using the following model:

# Y[jl] ~ Normal(∂[jl], V[jl])    j = 1, ..., J[l] studies and l = 1, ... L (study type)
# ∂[jl] ~ Normal(ø[l], psi^2[l])    ø[l] ~ Normal[d, tau^2]
# psi^2[l] ~ [-, -]  d ~ [-, -]  tau^2 ~ [-, -]

# where Y[jl] and V[jl] are the effect size and variance in the jth study of type l. In the
# context of the EFM example, the effect sizes are RDs, and l = 1, 2, 3, which relate to RCTs,
# comparative cohort studies and before and after studies, respectively. d is the overall 
# pooled effect over all sources of evidence, and tau^2 is the between study type variance.
# ø[l] is the pooled effect within study type l, and psi^2 is the between study variance 
# within study type l. ∂[jl] is the true underlying effect estimated by the observed Y[jl] 
# with variance V[jl]. When fitted, psi^2, d, and tau^2 require prior distributions with the 
# notation indicating these are to be specified by the user. 

# The code used to fit this model can be split into four sections. The first three sections 
# consider the RCTs, cohort studies and before and after studies respectively and the code is 
# identical in structure for all three. Data are entered as RDs and associated standard errors
# are combined assuming studies of the same design that are exchangeable using a Random 
# effects model. For each type of evidence, a different pooled mean (theta[1], theta[2, 
# theta[3]]) and between study (within type) standard deviation (sd.theta[1], sd.theta[2], 
# sd.theta[3]) is estimated. In the fourth section of the code, these three pooled estimates
# are assumed to be themselves exchangeable, with overall mean and standard deviation.

model_String <- "
model {
 
 # Randomised Control Trials:
 for (i in 1:R) {
 # Likelihood
  rct.rd[i] ~ dnorm(rct.psi[i], rct.prec[i])
  # Random effects model
  rct.psi[i] <- theta[1] + (rct.z[i] * sd.theta[1])
  # Prior
  rct.z[i] ~ dnorm(0, 0.001)
  # Transformation
  rct.prec[i] <- 1 / (rct.serd[i] * rct.serd[i])
 }
 
 # Comparative cohort studies:
 for (i in 1:C) {
# Likelihood
  coh.rd[i] ~ dnorm(coh.psi[i], coh.prec[i])
  # Random effect model
  coh.psi[i] <- theta[2] + (coh.z[i] * sd.theta[2])
  # Prior
  coh.z[i] ~ dnorm(0, 0.001)
  # Transformation
  coh.prec[i] <- 1 / (coh.serd[i] * coh.serd[i])
 }
 
 # Before and after studies:
 for (i in 1:B) {
 # Likelihood
  ba.rd[i] ~ dnorm(ba.psi[i], ba.prec[i])
  # Random effects model
  ba.psi[i] <- theta[3] + (ba.z[i] * sd.theta[3])
  # Prior
  ba.z[i] ~ dnorm(0, 0.001)
  # Transformation
  ba.prec[i] <- 1 / (ba.serd[i] * ba.serd[i])
 }
 
 # Combining all 3 sources of information:
 for (i in 1:T) {
 # Regression effect model
  theta[i] <- mean + (u[i] * sd.mean)
  u[i] ~ dnorm(0, 1 / u.prec[i]^2)
  # Prior
  u.prec[i] ~ dunif(0, 1)
  # Prior
  sd.theta[i] ~ dnorm(0, 0.001)T(0, )
  # Transformation
  var.theta[i] <- sd.theta[i] * sd.theta[i]
  # Transformation
  prec.theta[i] <- 1 / (sd.theta[i] * sd.theta[i])
  
 }
 # Hyperpriors:
 mean ~ dnorm(0, mean.prec)
 sd.mean ~ dnorm(0, 1 / sd.prec^2)T(0, )
 var.mean <- sd.mean * sd.mean
 prec.mean <- 1 / (sd.mean * sd.mean)
 sd.prec ~ dunif(0, 1)
 mean.prec ~ dunif(0, 1)
 
 mean.Pr <- exp(mean)
 }
"
writeLines(text = model_String, con = "GES_chap11.txt")


# Data:
data_JAGS <- list(R = 9, C = 7, B = 10, T = 3,
                  rct.rd = c(-10.51552, -2.028398, 4.115085, 6.479482, 0.0078609, 0, 
                             2.247191, -5.817028, -3.984064),
                  rct.serd = c(4.762193, 2.871006, 7.142432, 5.032322, 0.8079891, 8.058098,
                               3.1075, 44.53912, 5.587013),
                  ba.rd = c(-4.036327, 2.304048, -.6941801, -3.186446, -7.431126, -1.458522,
                            -4.036984, -1.613824, -1.461775, -.1177738),
                  ba.serd = c(2.242277, 3.579612, 0.6056279, 0.9381518, 2.121014, 0.5100973,
                              1.072718, 0.6358061, 0.507642, 0.1981163),
                  coh.rd = c(-1.41, -2.19, -4.34, -2.84, -2.53, -.23, -.46),
                  coh.serd = c(1.433, 4.71, 1.914, 1.052, 3.081, 0.232, 0.123)
                  )
data_JAGS

# Initial chain values:
inits <- list(
 list(
 mean = -.001, sd.mean = 0.001, sd.theta = c(0.001, 0.001, 0.001)
 ),
 list(
  mean = -.001, sd.mean = 0.001, sd.theta = c(0.001, 0.001, 0.001)
  ),
  list(
 mean = -.001, sd.mean = 0.001, sd.theta = c(0.001, 0.001, 0.001)
 ),
 list(
  mean = -.001, sd.mean = 0.001, sd.theta = c(0.001, 0.001, 0.001)
  )
 )
inits

# Parameters to monitor:
params <- c("mean.Pr", "mean", "sd.mean", "theta[1]", "theta[2]", "theta[3]")

n.iter <- 10000
n.burnin <- 5000
n.thin <- floor((n.iter-n.burnin)/500)

# Model run:
jags_Model <- jags(data = data_JAGS, model.file = "GES_chap11.txt",
                   parameters.to.save = params, n.chains = 4, n.iter = n.iter, 
                   n.burnin = n.burnin, n.thin = n.thin, jags.seed = 22)
jags_Model

autoJAGS_Update <- autojags(jags_Model, Rhat = 1.001, parallel = TRUE, n.cores = 4,
                            jags.seed = 22)
autoJAGS_Update

# Posterior visual checks:
posterior <- as.array(autoJAGS_Update$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("mean", "theta[1]", "theta[2]", "theta[3]"), 
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[,, 1:4], window = c(250, 500), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("mean", "theta[1]", "theta[2]", "theta[3]"))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("mean", "theta[1]", "theta[2]", "theta[3]")
         , lags = 50)

# End file ----------------------------------------------------------------