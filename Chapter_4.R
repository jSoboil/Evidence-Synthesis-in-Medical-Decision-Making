# ===========================================================================================
# Simple Fixed Effects Evidence Synthesis Example -------------------------
# ===========================================================================================
# Example found in Evidence Synthesis for Medical Decision Making, Example 4.1, pages 78-81.

# Set working directory - best practice for when saving model files etc.
# setwd(getwd())

# Load jags package:
library(rjags)
library(R2jags)
# Load post-inspection utilities:
library(bayesplot)
library(ggplot2)
library(magrittr)

# Declare model logic in a string: 
modelString <- "model {
# Sampling model/likelihood:
 for (i in 1:Nstud) {
  P[i] <- 1/V[i]
  Y[i] ~ dnorm(d, P[i])
 }
 # Prior:
 d ~ dnorm(0, 1.0e-5)
 
 # Converting pooled LOR (d) to standard OR:
 OR <- exp(d)
}"

# Convert and save object string to current directory. Note: make sure directory is set, 
# otherwise R will not be able to find it, as it saves file outside the current environment:
writeLines(text = modelString, con = "ExModel.txt")

# Import mean data variables to be used in simulation:
jags_Data <- list(Y = c(-.3289, -.3845, -.2196, -.2222, -.2255, .1246, -.1110), 
                  V = c(.0389, .0412, .0205, .0648, .0352, .0096, .0015), Nstud = 7)

# Declare starting Markov Chain Monte Carlo values. I.e. where the initial sampler starts in 
# the 'sample space'. Note that starting values can also be given randomised distributions, 
# which is best practice when conducting large, complex models to help avoid autocorrelation:
inits <- list(d = 0)

# Run Jags model, setting the number of simulations to 10000, throwing away 5000 samples from 
# this, due to conditional sampling, i.e. we want to simulate independent sampling. Momnitor
# the variables of interest, which in this case is d, the overall pooled effect size. The OR
# is also monitored:
EviSynth_Mod <- jags(data = jags_Data, model.file = "ExModel.txt", 
     n.iter = 10000, n.burnin = 5000, parameters.to.save = c("OR", "d"))

# ... thus this is the 'average' marginal effect/averaged over 10000 simulations, given the 
# declared normal probability distributions and the likelihood/sampling model. Note that due
# to the very simple model and vague prior, little to none deviance occurs (see sd for 
# declared sd for d in model, sigma equivalent to .00001, thus when marginalised/averaged out
# over 10000 simulations, this has a minor effect).

# Convert Jags object to MCMC object for post-inspection:
mcmc_Output <- as.mcmc(EviSynth_Mod$BUGSoutput)

plot(mcmc_Output)
mcmcplot(mcmc_Output)

# ===========================================================================================
# Simple Random Effects Evidence Synthesis Example ------------------------
# ===========================================================================================
# Since variances cannot go negative, a Normal distribution with mean zero and large variance
# is not a viable option:

#   tau ~ dunif(0, 10)

# ... a value of 10 for the between study standard deviation is very large on the LOR scale 
# and thus this prior distribution covers all plausible values.

model_string <- "model {
  for (i in 1:nStud) {
    # Transforming V into precision variable:
    P[i] <- 1 / V[i]
    # Likelihood:
    Y[i] ~ dnorm(delta[i], P[i])
    
    # Random sampling model for between 
    # study variability:
    delta[i] ~ dnorm(d, prec)
  }
  # Prior for mean effect (d):
  d ~ dnorm(0, 1.0E-5)
  # Converting LOR back to OR:
  OR <- exp(d)
  
  # Prior for between study sd:
  tau ~ dunif(0, 10)
  # Transforming tau into var
  # and precision for delta: 
  tau.sq <- tau * tau
  prec <- 1 / (tau.sq)
}"

writeLines(text = model_string, con = "RandomEffModel.txt")

# Data:
jags_Data <- list(Y = c(-.3289, -.3845, -.2196, -.2222, -.2255, .1246, -.1110), 
                  V = c(.0389, .0412, .0205, .0648, .0352, .0096, .0015), nStud = 7)


# Initial values:
inits <- list(d = 0, tau = 1, delta = c(0, 0, 0, 0, 0, 0, 0))

R_effects_Mod <- jags(data = jags_Data, model.file = "RandomEffModel.txt", 
     n.iter = 30000, n.burnin = 10000, 
     parameters.to.save = c("OR", "d", "tau.sq", "tau"))
R_effects_Mod


# Posterior Inspection via Visualisation -----------------------------------------
posterior <- as.array(R_effects_Mod$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("OR", "d"), 
           facet_args = list(ncol = 1, strip.position = "left"))


color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[,, 1:5], window = c(950, 1000), size = 1) +
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("OR", "d"))

color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("OR", "tau"),
           off_diag_args = list(size = 1.5))
# Note that in this example, there is considerable uncertainty in the estimation of the 
# between study variance parameter, this is typical when the number of studies in the 
# meta-analysis is small. This can be visually seen in the large spread of tau, forming a
# large, right skewed fat tail of uncertainty.

# End file ----------------------------------------------------------------