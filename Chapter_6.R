# ====================================================================================
# Libraries ---------------------------------------------------------------
# ====================================================================================
# This chapter is concerned with model critique and evidence consistency in random effects
# meta-analysis
library(R2jags)
library(rjags)
library(bayesplot)
library(tidyverse)
library(parallel)

options(mc.cores = detectCores())
set.seed(20200407)
sessionInfo()

# ====================================================================================
# The Random Effects model revisited --------------------------------------
# ====================================================================================
# This is generally a more efficient way to parse data into JAGS. We require two new 
# variables, si and ti, which indicate which study number and which treatment arm is 
# represented by row i. So, row 1 is treatment arm 2 (magnesium) of study 1, row 2 is
# treatment arm 2 (placebo) of study 1, row 3 is treatment 2 (magnesium) of study 2, and 
# so on. This format allows for easy extension to trials with three or more arms and 
# analysis of more than two treatment options.

# Data:
s <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 
       12, 12, 13, 13, 14, 14, 15, 15, 16, 16)

n <- c(40, 36, 135, 135, 200, 200, 48, 46, 150, 148, 59, 56, 25, 23, 22, 21, 76, 75,
       27, 27, 89, 80, 23, 33, 130, 122, 1159, 1157, 107, 108, 29011, 29039)

r <- c(1, 2, 9, 23, 2, 7, 1, 1, 10, 8, 1, 9, 1, 3, 0, 1, 6, 11, 1, 7, 2, 12, 5, 13, 
       4, 8, 90, 118, 4, 17, 2216, 2103)

t <- rep(c(2, 1), 16)

# Bind data:
jags_Data <- list(s = s, t = t, r = r, n = n)
as.data.frame(jags_Data)
# The Random Effects model can be written as:

# r[i] ~ Binomial(p[i], n[i]), for i = 1, ..., 31

# logit(p[i]) = µ[s[i]] + ∆ * I[t[i = 2]], where I[t[i = 2]] 
                                           # := if 1 t[i] = 2; if 0 t[i] = 1

# Because we have a generic likelihood statement for each arm, the logistic regression 
# equation needs to pick out the right baseline and treatment parameters. This is done 
# using nested indexing. µ[s[i]] picks out µ1 when s[i] = 1 (i.e. rows i = 1 and 2), 
# µ2 when s[i] = 2 (i.e. rows i = 3 and 4) and so on. Similarly for ∆[s[i]], however we 
# only add on a treatment effect for the magnesium arm (t[i] = 2), which is achieved by
# multiplying ∆[s[i]] by zero when row i represents the placebo arm (t[i] = 1). The 
# indicator function I[t[i = 2]] returns a 1 when the data represent the active arm 
# (t[i] = 2), and a 0 otherwise.

# The model for study level random effects ∆[j] is exactly the same as before:

# ∆[j] ~ Normal(d, tau.sq), j = 1, ..., 16

# and we give Normal priors to baseline and mean treatment effect parameters, and a 
# Half-Normal prior to between study standard deviation:

# µ[j] ~ Normal(0, 1000^2); d ~ Normal(0, 1000^2); tau ~ Half-Normal(0, 1000^2)

# The WinBUGS code makes use of the equals(t[i],2) function that returns a 1 if t[i] = 2 
# and a 0 otherwise (in place of the indicator function I[t[i = 2]]):
Model_string <- "model {
    for (i in 1:32) {
     # Likelihood:
     r[i] ~ dbin(p[i], n[i])
      
      # Logistic regression model:
       logit(p[i]) <- mu[s[i]] + delta[s[i]] * equals(t[i], 2)
    }
  for (j in 1:16) {
   # Hierarchical random 
   # effects model:
    delta[j] ~ dnorm(d, prec) 
   
 # Prior for baseline
 # random effect, µ:
  mu[j] ~ dnorm(0, 1.0e-6)
 }
 # Prior for Pop. 
 # treatment effect:
  d ~ dnorm(0, 1.0e-6)
 # Population OR:
  OR <- exp(d)
 # Prior for btw. studies
 # sd:
  tau ~ dnorm(0, 1.0e-6)T(0, )
 # Variance:
  tau.sq <- (tau * tau)
 # Precision:
  prec <- 1 / (tau.sq)
 
 # Replicate LOR for 
 # prediction:
  delta.new ~ dnorm(d, prec)
 # Mean random
 # effect:
  delta[19] <- exp(d)
 # Predictive dist.
 # for plotting:
  delta[20] <- delta.new
}
"
writeLines(text = Model_string, con = "Model.6.txt")

inits <- list(
 list(d = 0, tau = 0.25, mu = c(rep(1, 16)), delta = c(rep(1, 16), rep(NA, 4)), 
      delta.new = 0), 
 list(d = 0, tau = 1, mu = c(rep(0, 16)), delta = c(rep(0, 16), rep(NA, 4)), 
      delta.new = 0), 
 list(d = 0, tau = 0.5, mu = c(rep(-1, 16)), delta = c(rep(-1, 16), rep(NA, 4)), 
      delta.new = 0)
 )

jags_Model <- jags(data = jags_Data, model.file = "Model.6.txt", 
                   parameters.to.save = c("d", "OR", "tau", "tau.sq", "prec"),
                   n.iter = 10000, n.burnin = 5000, n.chains = 3, 
                   inits = inits)
jags_Model

posterior <- as.array(jags_Model$BUGSoutput$sims.array)
dimnames(x = posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("d", "tau"), 
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("d", "tau"), 
           size = 1, ) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("d", "tau"))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("d", "prec", "tau.sq"), lags = 10)

# ====================================================================================
# Assessing model fit -----------------------------------------------------
# ====================================================================================
# Instead of point estimates of these RSS or residual deviance, we will obtain posterior
# distributions for them. These can be summarised in various ways, but typically the 
# posterior mean of these distributions is reported.

# Deviance ----------------------------------------------------------------
# The deviance statistic measures the fit of the predictions made by a particular model 
# to the observed data using the likelihood function. The likelihood function measures 
# how ‘likely’ observed data are given a particular model and so is a natural measure to 
# focus on to assess model fit - in fact Frequentist statistics revolves around 
# maximising the likelihood function. Deviance, D[model], is defined as:

# D[model] = -2 * Loglik[model]

# For a given model and observed data, the larger the likelihood then the closer the 
# model fit. Similar for log-likelihood. Multiplying by -2 reverses this, so the smaller 
# the deviance, D[model], then the closer the model fit. D[model] measures how far the 
# model predictions deviate from the observed data. The deviance is simply a function of 
# model parameters that can be written down and calculated for each iteration of a Markov
# Chain Monte Carlo (MCMC) simulation. Posterior summaries for D[model] can then be 
# obtained as for other parameters.

# Residual deviance -------------------------------------------------------
# The smaller the deviance statistic, D[model], then the better the model fit. But how 
# small is small? The disadvantage of using the raw deviance statistic, D[model], is that
# there is no clear answer to this question.

# Instead we define the residual deviance, D[res], which helps us gauge how good the 
# model fit is, by providing a reference point. The residual deviance is equal to the 
# deviance for a given model, D[model], minus the deviance for a saturated model, D[sat]:

# D[res] = D[model] - D[sat]

# Note: a saturated model is one where all of the predictions from the model are equal to
# the observed data values.

# There is no pre-set node for the residual deviance, and it is not included on the DIC 
# tool. Instead we need to write out the relevant formula for a given likelihood, D[res], 
# as new nodes in the BUGS/JAGS code. See residual deviance formulae for likelihood 
# functions.

# The BUGS/JAGS code to calculate this consists of two lines of code, the first to 
# calculate the individual contributions to the residual deviance (within a loop over i),
# and the second to sum over these:
Model_string <- "model {
    for (i in 1:32) {
     # Likelihood:
     r[i] ~ dbin(p[i], n[i])
     
     # Predicted model 
     # deviance:
     rhat[i] <- p[i] * n[i]
     # Deviance of each 
     # data point:
     dev[i] <- 2 * (r[i] * (log(r[i]) - log(rhat[i]))
     + (n[i] - r[i]) * (log(n[i] - r[i]) - log(n[i] - rhat[i])))
     
      # Logistic model:
       logit(p[i]) <- mu[s[i]] + delta[s[i]] * equals(t[i], 2)
    }
    # Residual deviance:
    resdev <- sum(dev[])
    
  for (j in 1:16) {
   # Hierarchical random 
   # effects model:
    delta[j] ~ dnorm(d, prec) 
   
 # Prior for mean
 # random effect, µ:
  mu[j] ~ dnorm(0, 1.0e-6)
 }
 # Prior for Pop. 
 # treatment effect:
  d ~ dnorm(0, 1.0e-6)
 # Population OR:
  OR <- exp(d)
 # Prior for btw. studies
 # sd:
  tau ~ dnorm(0, 1.0e-6)T(0, )
 # Variance:
  tau.sq <- (tau * tau)
 # Precision:
  prec <- 1 / (tau.sq)
 
 # Replicate LOR for 
 # prediction:
  delta.new ~ dnorm(d, prec)
 # Mean random
 # effect:
  delta[19] <- exp(d)
 # Predictive dist.
 # for plotting:
  delta[20] <- delta.new

}
"
writeLines(text = Model_string, con = "Model.6.txt")

inits <- list(
 list(d = 0, tau = 0.25, mu = c(rep(1, 16)), delta = c(rep(1, 16), rep(NA, 4)), 
      delta.new = 0), 
 list(d = 0, tau = 1, mu = c(rep(0, 16)), delta = c(rep(0, 16), rep(NA, 4)), 
      delta.new = 0), 
 list(d = 0, tau = 0.5, mu = c(rep(-1, 16)), delta = c(rep(-1, 16), rep(NA, 4)), 
      delta.new = 0)
 )

jags_Model <- jags(data = jags_Data, model.file = "Model.6.txt", 
                   parameters.to.save = c("d", "OR", "tau", "resdev"),
                   n.iter = 10000, n.burnin = 5000, n.chains = 3, 
                   inits = inits)
jags_Model
# There are 32 unconstrained data points in this example (16 studies # 2 arms), so we 
# would expect the posterior mean of the residual deviance to be close to 32 if the 
# model predictions fit the data well. There is no evidence of lack of fit for the 
# Random Effects model; resdev ≈ 29.5

posterior <- as.array(jags_Model$BUGSoutput$sims.array)
dimnames(x = posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("d", "resdev"), 
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("d", "tau"), 
           size = 1, ) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("d", "resdev"))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("d", "resdev"), lags = 10)

# ====================================================================================
# Model comparison --------------------------------------------------------
# ====================================================================================
# A model selection process systematically compares the fit of a set of models. For 
# example, we can simplify our Random Effects model by reducing it to a Fixed Effect 
# model (model contraction), or we could complicate it further by adding in covariates 
# in a meta-regression (model expansion). Both the deviance and residual deviance 
# statistics can be used to compare the fit of different models. In fact, because the 
# two measures differ only by a constant term, Dsat, which does not depend on the model 
# fitted, then when we look at the difference in either measure between two models, we 
# get exactly the same result, i.e.:

# D[res1] - D[res2] = (D[model1] - D[sat]) - (D[model2] - D[sat]) = D[model1] - D[model2]

# Model fit is not the only consideration to make when selecting a model. The more
# parameters we include (i.e. the more complicated the model gets), obviously the better 
# the model fit will be. However, problem of bias...

# Instead, we would like to select a model that remains as simple as possible for good
# predictive power, whilst still fitting well to the observed data – the parsimony 
# principle.

# The DIC extends the AIC to handle non-nested hierarchical models by defining the 
# effective number of parameters, pD:

# pD = mean(D[model]) - D(estimated theta)

# The effective number of parameters is the posterior mean deviance for a given model 
# minus the deviance calculated at some plug-in value for the parameters. However, when 
# the relationship between the model predictions and the basic model parameters is 
# highly nonlinear (as is often the case in the multi-parameter evidence syntheses), 
# then it is more appropriate to use the posterior mean of the model predictions as the 
# plug-in, and calculate pD externally in R.

# Deviance Information Criteria -------------------------------------------
# The DIC is equal to the posterior mean deviance plus the effective number of 
# parameters, pD. It has been suggested that differences in DIC over 5 are important, 
# whereas if there are only small differences (less than 3) in DIC there is probably 
# little to choose between two models – although one should check robustness of 
# conclusions to choice of model.

# Model consitency --------------------------------------------------------
# In this section we describe cross-validation which can be viewed as a gold- standard 
# approach to assessing consistency, then go on to look at mixed predictive p-values as
# an alternative method that approximates cross-validation.

# Cross-validation extends naturally to the Bayesian context where the predictions for 
# the omitted data are based on the posterior predictive distributions.

# The key step is how to make predictions for the omitted data. We shall restrict 
# attention to the situation where a single data point is omitted. For a Fixed Effect
# model we would expect the true treatment effect in an omitted study to be the fixed 
# effect ∆[new] = d, whereas in a Random Effects model we would expect the true treatment
# effect to be drawn from the random effects distribution of effects, i.e. the predictive
# distribution:

#   ∆[new] ~ Normal(d, tau.sq)

# The predictive distribution represents two different sources of uncertainty: 
# uncertainty in the mean value, d, and uncertainty as to where in the random effects 
# distribution the new study will lie. This predicts the treatment effect we would 
# expect to see in an infinitely sized new study. However, the omitted study has a (known) 
# finite sample size on each arm, and so there will be an additional element of 
# uncertainty in the observed treatment effect that is due to sampling error. Let 
# the observed data for the omitted study be r[c] out of n[c] on the control arm and r[t]
# out of n[t] on the treatment arm. The predicted number of responders, r[t[new]], is
# drawn from a binomial distribution with the predicted probability of response on the 
# treatment arm, p[t[new]] and known denominator, n[t]:

# r[t[new]] ~ Binomial(p[t[new]], n[t])

# This describes the sampling error that we would expect in a study with n[t] on the 
# treatment arm. The predicted probability of response, p[t[new]], follows the logistic
# regression relationship:

# logit(p[t[new]]) = logit(p[c[new]]) + ∆[new]

# The predicted log-odds on the treatment arm is equal to the predicted log-odds on the 
# control arm plus the predicted treatment effect. We can use the observed data to 
# characterise the baseline probability of response and its uncertainty. We draw the 
# baseline probability of response on the control arm, p[c[new]], from a Beta 
# distribution:

# p[c[new]] ~ Beta(r[c], n[c] - r[c])

# The Beta distribution describes the uncertainty in the estimate of a proportion for a 
# given number of responders, r[c], and nonresponders, (n[c] - r[c]).

# The predicted number of responders, r[t[new]], can, from the above, then be compared 
# with the observed number, r[t]. If the observed data, r[t], is supported by the 
# posterior distribution for r[t[new]], then that data point is consistent with the model
# predictions based on the remaining data alone. We can compare r[t] with r[t[new]] by 
# forming a Bayesian p-value, which measures the probabilitiy that r[t[new]] exceeds the
# observed value r[t], Pr(r[t[new]] > r[t]). This is achieved using the step(e) function,
# which records a 1 if its argument e ≥ 0 and a zero otherwise. So, the BUGS/JAGS code

# p.crossval <- step(r.new - r)

# creates a node which is a string of 0’s (on iterations where r[t[new]] < r[t]) and 1's
# (on iterations where r[t[new]] > r[t]). The posterior mean of the node p.crossval gives
# the posterior probability that the model predictions exceed the observed value. Because
# the number of responders is discrete (i.e. can only take whole number values), we could
# get the situation where the model predictions and observed data are equal, i.e. 
# (r[t[new]] == r[t]). In such a case, we can make a continuity correction, where we only
# record 0.5 rather than 0 or 1 on iterations where r[t[new]] == r[t]. The BUGS/JAGS
# code is:

# p.crossval <- step(r.new - r) - .5 * equals(r.new, r)

# which uses another BUGS/JAGS function equals(a,b), which records a 1 if its arguments 
# are equal (a = b) and a zero otherwise.

# Example 6.1 revisited ---------------------------------------------------
# We are interested in whether the ISIS-4 trial results are consistent with the results 
# from the other 15 trials, so we perform cross-validation on ISIS-4.

# ISIS-4 can be omitted from the random effects meta-analysis model in BUGS/JAGS by 
# simply changing the loop over data from 1:32 to 1:30 and the loop over study from 1:16 
# to 1:15. We do not need to delete the data from the data array. We will also need to 
# adjust the initial values, as we only need 15 deltas and mus rather than 16.

#  Note that the observed data from the omitted trial ISIS-4 is r[31] and n[31] from the 
# control arm and r[32] and n[32] from the magnesium arm:
Model_string <- "model {
    for (i in 1:30) {
     # Likelihood:
     r[i] ~ dbin(p[i], n[i])
      
      # Logistic regression model:
       logit(p[i]) <- mu[s[i]] + delta[s[i]] * equals(t[i], 2)
    }
  for (j in 1:15) {
   # Hierarchical random 
   # effects model:
    delta[j] ~ dnorm(d, prec) 
   
 # Prior for baseline
 # random effect, µ:
  mu[j] ~ dnorm(0, 1.0e-6)
  }
 
 # Prior for Pop. 
 # treatment effect:
  d ~ dnorm(0, 1.0e-6)
 # Population OR:
  OR <- exp(d)
 # Prior for btw. studies
 # sd:
  tau ~ dnorm(0, 1.0e-6)T(0, )
 # Variance:
  tau.sq <- (tau * tau)
 # Precision:
  prec <- 1 / (tau.sq)
 
 # Replicate LOR for 
 # prediction:
  delta.new ~ dnorm(d, prec)
 # Mean random
 # effect:
  delta[19] <- exp(d)
 
 # Estimated no. deaths
 # in control group:
  a <- r[31]
 # Estimated no. survivors
 # in control group:
  b <- n[31] - r[31]
 # Draw on new control
 # probability:
  p[31] ~ dbeta(a, b)
 
 # Form new treatment
 # probability:
  logit(p[32]) <- logit(p[31]) + delta.new
  
 # Draw new no. of deaths
 # in treatment group:
  r.new ~ dbin(p[32], n[32])
 # Record whether predicted no.
 # of deaths exceeds observed:
 p.crossval <- step(r.new - r[32]) - .5 * equals(r.new, r[32])

}
"
writeLines(text = Model_string, con = "Model.6.2.txt")

inits <- list(
 list(d = 0, tau = 0.25, mu = c(rep(1, 15)), delta = c(rep(1, 15), rep(NA, 4)), 
      delta.new = .5), 
 list(d = 0, tau = 1, mu = c(rep(0, 15)), delta = c(rep(0, 15), rep(NA, 4)), 
      delta.new = 0), 
 list(d = 0, tau = 0.5, mu = c(rep(-1, 15)), delta = c(rep(-1, 15), rep(NA, 4)), 
      delta.new = -.5)
 )

jags_Model <- jags(data = jags_Data, model.file = "Model.6.2.txt", 
                   parameters.to.save = c("p.crossval", "r.new", "p[32]"),
                   n.iter = 10000, n.burnin = 5000, n.chains = 3, 
                   inits = inits)
jags_Model

posterior <- as.array(jags_Model$BUGSoutput$sims.array)
dimnames(x = posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("r.new", "p[32]"), 
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("p[32]", "r.new"), 
           size = 1, ) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("p[32]", "r.new"))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("p[31]", "r.new"), lags = 10)

# The cross-validation p-value is p.crossval = 0.073, which shows that r.new only 
# exceeds r[32] about 7% of the time. There is therefore only weak evidence that the 
# results of ISIS-4 are inconsistent with results from the other 15 studies for a Random
# Effects model. This result is due to the wide credible interval for the predictive 
# distribution (as a result of the high level of heterogeneity), which still just 
# includes the ISIS-4 result even when based on only the remaining 15 studies. Note that 
# evidence of inconsistency does not identify which data point is ‘wrong’. For example 
# here, it is perfectly feasible that the ISIS-4 result is showing the real treatment 
# effect, whereas the other 15 studies are in some way biased.

# ====================================================================================
# Mixed Predictive Checks -------------------------------------------------
# ====================================================================================
# Cross-validation quickly becomes computationally expensive if we want to find 
# cross-validation p-values for each data point in our meta-analysis. For instance in the
# magnesium example we would have to run 16 separate models for each of the 16 trials.

# An alternative is to find a way to calculate approximate cross-validation p-values in a
# single run of the MCMC simulation, so that we only need to run the model once. Marshall
# and Spiegelhalter proposed mixed predictive p-values when working with hierarchical 
# models as approximations to the cross-validation p-values.

# The idea is to fit a random effects meta-analysis model to all of the data points. The 
# resulting model is used to make predictions for each of the individual observed data 
# points, which are then compared with the observed data value to form a p- value. Again,
# the key step is in forming the model predictions. Again, a naive approach would be to 
# use the study specific treatment and baseline parameters estimated from the model to 
# predict the observed data:

# r[i[post]] ~ dbin(p[i], n[i])

# where the probability of response for data point i is given by a logistic regression.
# The posterior probability that r[i[post]] exceeds r[i], Pr(r[i[post]] > r[i]), is known
#  as the posterior predictive p-value. However it gives a very conservative estimate of
# the inconsistency between data point i and the remaining evidence.

# A less conservative approach is to predict the true treatment effect in each study 
# using the predictive distribution for a ‘new’ study:

# ∆[new] ~ dnorm(d, tau.sq)

# as we used in cross-validation. This removes some of the influence of the individual 
# data point on its predicted value – although the observed data will still have some 
# influence on the mean and variance of the random effects distribution. We use the 
# estimated baseline for each study from the model, so that the predicted probability of 
# response is:

# logit(p[i[new]]) = µ[s[i]] + ∆[new] * I[t[i] = 2]; where I[t[i] = 2], 1 = t[2] and
# 0 = t[1]

# This is identical to the logistic regression used in the previous exercise, *except* 
# that the study specific treatment effect is replaced by the treatment effect predicted 
# in a new distribution from the same random effects population.

# The predicted number of responders, r[i[mxd]], is drawn from a Binomial distribution
# with the predicted probability response, p[i[new]], and known denominator, n[i]:

# r[i[mxd]] ~ dbin([p[i[new]], n[i])

# These will still be conservative compared with the ‘gold-standard’ method of cross-
# validation, however can be used as a first step to get an indication of the consistency
# of the evidence, and guide individual data points on which to run cross-validation. 
# Because we obtain p-values for each study, we are in danger of incorrectly interpreting
# the significance of these multiple hypothesis tests. If we calculate several p-values 
# on the same data, we would expect these to follow a Uniform distribution on the 
# interval (0,1). Plotting the ordered p-values for each study against the relevant 
# Uniform order statistics provides a tool to gauge whether any of the p-values are 
# unusually small or large.

# Example 6.1 revisited Magnesium vs placebo for MI ----------------------
# The BUGS/JAGS code to calculate the mixed predictive p-values is:
Model_string <- "model {
    for (i in 1:32) {
     # Likelihood:
      r[i] ~ dbin(p[i], n[i])
      
     # Logistic regression model:
      logit(p[i]) <- mu[s[i]] + delta[s[i]] * equals(t[i], 2)
        
       # Predicted likelihood
       # for deaths:
        r.mxd[i] ~ dbin(p.new[i], n[i])
       
       # Mixed predictve 
       # p-value:
        p.mxd[i] <- step(r.mxd[i] - r[i]) - 0.5 * equals(r.mxd[i], r[i])
        
       # Predicted probability
       # of death:
        logit(p.new[i]) <- mu[s[i]] + delta.new * equals(t[i], 2)

    }
  for (j in 1:16) {
   # Hierarchical random 
   # effects model:
    delta[j] ~ dnorm(d, prec) 
   
 # Prior for baseline
 # random effect, µ:
  mu[j] ~ dnorm(0, 1.0e-6)
  }
 
 # Prior for Pop. 
 # treatment effect:
  d ~ dnorm(0, 1.0e-6)
 # Population OR:
  OR <- exp(d)
 # Prior for btw. studies
 # sd:
  tau ~ dnorm(0, 1.0e-6)T(0, )
 # Variance:
  tau.sq <- (tau * tau)
 # Precision:
  prec <- 1 / (tau.sq)
 
 # Predicted average 
 # treatment effect:
  delta.new ~ dnorm(d, prec)
  
 # Mean random
 # effect:
  delta[19] <- exp(d)

}
"
writeLines(text = Model_string, con = "Model.6.3.txt")

inits <- list(
 list(d = 0, tau = 0.25, mu = c(rep(1, 16)), delta = c(rep(1, 16), rep(NA, 3)), 
      delta.new = .5), 
 list(d = 0, tau = 1, mu = c(rep(0, 16)), delta = c(rep(0, 16), rep(NA, 3)), 
      delta.new = 0), 
 list(d = 0, tau = 0.5, mu = c(rep(-1, 16)), delta = c(rep(-1, 16), rep(NA, 3)), 
      delta.new = -.5)
 )

jags_Model <- jags(data = jags_Data, model.file = "Model.6.3.txt", 
                   parameters.to.save = c("p.new"),
                   n.iter = 10000, n.burnin = 5000, n.chains = 3, 
                   inits = inits)
jags_Model
# ISIS-4 has a mixed predictive p-value of 0.073, which as we expect is bigger (i.e. more
# conservative) than the cross-validation p-value of 0.04.

# End file ----------------------------------------------------------------