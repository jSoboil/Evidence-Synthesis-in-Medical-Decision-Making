# ==========================================================================================
# Markov models: Modelling study effects ----------------------------------
# ==========================================================================================
# Unfortunately, when modelling transition parameters, it is rare that we can work directly
# with transition probabilities. Although logistic regression models that are typically used
# with binary data can be extended to multinomial outcomes, the necessity for transition 
# probabilities to sum to 1 leads to model parameters that do not have a natural 
# interpretation. Instead, it is much more convenient to model transition rates, which are
# only constrained to be positive. Log-linear models can then be used to model the rates, 
# and regression parameters can be interpreted as hazard ratios.

# Evidence reported in an event history format directly informs rate parameters through 
# Poisson likelihoods. Evidence reported in aggregate format informs transition 
# probabilities for a given observation cycle length, which in turn can be written down as 
# functions of the underlying rate parameters.

# The transition rates are the basic parameters, and the transition probabilities are the 
# functional parameters. Treating transition rates, rather than transition probabilities, as
# the basic parameters has the added advantage that we can fit models with fewer parameters.

# We use as an example the three-state forwards model. There are three transition rate 
# parameters, gamma[1, 2], gamma[1, 3], and gamma[2, 3]. The transition probability 
# parameters can be written in terms of these three rates, for given observed cycle length, 
# T, using solutions to Kolmorgorov’s forward equations, such that, for example, moving
# from a rate to probability is

# P = 1 - exp(- r * t)

# where r is the rate and t is the time period. In other words the transition probability 
# matrix can be obtained by taking the exponential of t times the transition rate matrix.

# It is useful to re-parameterise the rates from state 1 as follows:

# gamma[1, 2] = lambda[1](1 - rho); gamma[1, 3] = lambda[1](rho)

# where lambda[1] is interpreted as the rate at which individuals leave state 1, and rho as
# the probability that individuals leaving state 1, go directly to state 3 rather than via 
# state 2.

# There are therefore three basic parameters, lambda[1], rho, and gamma[2, 3], for which we
# will from which we can derive all other functional parameters (transition rates and
# probabilities). As in standard meta-analysis, we can fit either Fixed Effect or Random 
# Effects models. However, there are now three parameters on which there may be random 
# study-specific effects. We shall illustrate how to fit a model where the rate that 
# individuals leave state 1, lambda[1, s], depends on study, s, but all other parameters,
# rho, and gamma[2, 3], are assumed to be fixed across studies. Other models for study 
# effects could be investigated, and compared using goodness of fit measures.

# *It is most natural to model rates on the log scale, so we put a Random Effects model on 
# the log rate that individuals leave state 1*:

#      log(lambda[1, s]) ~ Normal(L, tau.sq)

# We give a flat Beta prior for rho, the probability of going to state 3 rather than 
# state 2, having left state 1. We give an Exponential prior to gamma[2, 3] the transition
# rate from state 2 to state 3. Finally, we give a flat Normal prior for the random effects 
# mean, L, and a uniform prior on tau, the between study standard deviation in the log rate 
# of leaving state 1.

# rho ~ Beta(1, 1); 
# gamma[2, 3] ~ Exponential(.001);
# L ~ Normal(0, 100^2)
# tau ~ Uniform(0, 1)

# The code for this model is as follows:
model_String <- "
 model {
 for (s in 1:5) {
  # Random effects model
  # on log(lambda[1]):
  loglam1[s] ~ dnorm(L, prec)
  
  # Derive lambda1 and transition 
  # rates 1 to 2 and 1 to 3:
  log(lambda1[s]) <- loglam1[s]
  gamma12[i] <- lambda1[i] * (1 - rho)
  gamma13[i] <- lambda[i] * (rho)
 }
 # Priors for...
 
 # Transition rate 2 to 3:
 gamma23 ~ dexp(.001)
 # Random effects mean:
 L ~ dnorm(0, 1.0e-5)
 # Between study sd:
 tau ~ dunif(0, 1)
 # Precision as 
 # an inverse function 
 # of var:
 prec <- 1 / (tau * tau)
 # Probability of tp to 3 
 # rather than 2, having 
 # left 1:
 rho ~ dbeta(1, 1)
}
"

# ==========================================================================================
# Synthesis of studies reporting aggregate data ---------------------------
# ==========================================================================================
# I did not finish this section as solving kolgomorov's forward equations in complex models
# can be daunting. Therefore, have moved on to section that includes studies with event
# history.











