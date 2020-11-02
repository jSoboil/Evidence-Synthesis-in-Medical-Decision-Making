library(rjags)
library(R2jags)
library(bayesplot)
library(tidyverse)

# ============================================================================================
# Expected Value of Information for research prioritisation and s --------
# ============================================================================================
# The underlying rationale for performing an evidence synthesis to inform a CEA is that the 
# resulting decision reflects the body of evidence available, producing joint parameter 
# estimates and corresponding uncertainties. Under uncertainty, we average over the joint
# distribution that reflects parameter uncertainty to obtain Expected (average) net benefit.

# Accordingly, the best decision in the face of uncertainty is to recommend the intervention
# with the highest E[NB]. However, if the E[NB] lost from choosing a suboptimal intervention
# is very small, there will be little to gain from reducing uncertainty in the decision. The
# value of carrying out further research therefore depends on both the uncertainty in the 
# decision, and the E[NB] lost from making the 'wrong' decision. This is the 'Expected Value
# of Information' (EVI) calculation, which helps to prioritise and design new research 
# studies.

# Note: methods for prioritising further research are intrinsically linked to evidence 
# synthesis and CEA. Moreover, the process is inherently Bayesian. The posteriors of an 
# analysis become the priors for an EVI analysis and for a new synthesis that incorporates 
# new evidence collected in a new study.

# The following chapter introduces Expected Value of Perfect Information (EVPI), which 
# measures the value of collecting infinite quantities of evidence on all parameters to 
# eliminate decision uncertainty; Expected Value of Partial Perfect Information (EVPPI), 
# which measures the value of collecting infinite quantities of evidence on just a subset of
# parameters, whilst retaining uncertainty in the remaining parameters; Expected Value of 
# Sample Information (EVSI) which measures the value of collecting evidence from a given 
# study design to reduce, but not eliminate, decision uncertainty; and Expected Net Benefit 
# of Sampling (ENBS), which measures the net value of running a particular study minus the 
# costs of such a study, thus providing a basis for determining optimal study design. This
# extends on the model originally introduced in Chapter 8. See original text for description
# of additional parameters.

# The decision problem is to determine the most cost-effective of the two screening 
# strategies; targeted testing where only pregnant women in the high risk groups are tested, 
# or universal testing. T is the cost of screening and M is the net maternal benefit of an
# early diagnosis. The incremental net benefit (INB) function for a population with N 
# pregnancies per year is the product along the decision tree of INBs minus the 
# Incremental costs multiplied by the probability of each branch on the tree.

# The parameters in the Net Benefit function can be classified into two groups, 
# epidemiology parameters and economic parameters. N, M, T are economic; a, b, e, h are 
# epidemiological.

# We set the number of pregnancies per year at N = 105000; the unit cost of an HIV test at 
# T = 3 British sterling; the net benefit of an early maternal diagnosis M, comes from the
# previous model:

# M = 6000012 - 54296 * Y; 
# where Y ~ dgamma(0.56, 3)T(0, 2)



model_String <- "
model {

  # Binomial likelihood:
  for (i in 1:12) {
   r[i] ~ dbin(p[i], n[i])
  
  }
  
 # Declared relationships
 # between basic and
 # functional parameters:
  p[1] <- a
  p[2] <- b
  p[3] <- c
  p[4] <- d
  p[5] <- (b * d + (1 - a - b) * e / (1 - a))
  p[6] <- (a * c + b * d + (1 - a - b) * e)
  p[7] <- (a * c * f) / ((a * c * f) + (b * d * g) + (e * h * (1 - a - b)))
  p[8] <- (b * d * g) / ((b * d * g) + (e * h * (1 - a - b)))
  p[9] <- ((f * c * a) + (g * d * b) + h * e * (1 - a - b)) / 
                                                 ((c * a) + (d * b) + e * (1 - a - b))
  p[10] <- g
  p[11] <- w
  p[12] <- ((b * d) + (w * e) * (1 - a - b)) / ((b * d) + e * (1 - a - b))
 
# Estimated basic
# prior parameters:
 a ~ dbeta(1, 2)
 c ~ dbeta(1, 1)
 d ~ dbeta(1, 1)
 e ~ dbeta(1, 1)
 f ~ dbeta(1, 1)
 g ~ dbeta(1, 1)
 h ~ dbeta(1, 1)
 w ~ dbeta(1, 1)
 z ~ dbeta(1, 1)

# Estiamted functional
# prior parameter:
b <- z * (1 - a)

# Distribution for NB of Maternal Diagnosis:
 M <- 60012 - 54296 * Y

# Maternal diagnosis:
Y ~ dgamma(0.56, 3)T(0, 2)

# Net Benefit:
nb[1] <- 0
nb[2] <- 105000 * (1 - a - b) * (M * e * (1 - h) - 3.0 * (1 - e * h))

}
"
writeLines(text = model_String, con = "Chapter_12.txt")

jags_data <- list(
 r = c(11044, 12, 252, 10, 74, 254, 43, 4, 87, 12, 14, 5),
 n = c(104577, 882, 15428, 473, 136139, 102287, 60, 17, 254, 15, 118, 31)
)


params <- c("nb[1]", "nb[2]")

jags_Mod <- jags(data = jags_data, model.file = "Chapter_12.txt", 
                 parameters.to.save = params, n.iter = 50000, 
                 n.burnin = 10000, n.chains = 2)
jags_Mod

# The E[NB] is positive, indicating that the optimal strategy based on current information 
# is universal testing, thus k = 2.

# ==========================================================================================
# Expected value of perfect information -----------------------------------
# ==========================================================================================
# EVPI measures the ave. gain in NB from knowing perfectly what values the parameters take.
# In other words, EVPI is the gain from eliminating decision uncertainty altogether. One can
# think this is equivalent to running an infinitely sized study. Although unrealistic, EVPI
# is still a useful calculation because: 1) it is straightforward to calculate, and 2) it
# provides an upper bound for the value of running a study of any given design. So, if there
# is little value in running an infinitely sized study, there will certainly be no value in
# running a smaller trial that only measures some of the model parameters.

# let NB(\theta, k) be the Net Benefit function for intervention k, given values for the 
# input parameters, \theta. Then EVPI is given by:

# EVPI = E[max{NB(\theta, k)}] - max{E[NB(\theta, k)]},

# where the first parameter is based on perfect information and the second on current 
# information. For simulation, however, it is better to transform the above into

# EVPI = E[max{NB(\theta, k) - NB(\theta, k)}]

# In the HIV example, we assume a 10-year horizon for the screening strategy, and so EVPI is
# multiplied by 10. However, it is common to discount returns in decision models. Using a
# discount rate of 5% we multiply EVPI by 7.7127 (rather than 10). in this example there are
# only two strategies, so the code is:

model_String <- "
model {

  # Binomial likelihood:
  for (i in 1:12) {
   r[i] ~ dbin(p[i], n[i])
  
  }
  
 # Declared relationships
 # between basic and
 # functional parameters:
  p[1] <- a
  p[2] <- b
  p[3] <- c
  p[4] <- d
  p[5] <- (b * d + (1 - a - b) * e / (1 - a))
  p[6] <- (a * c + b * d + (1 - a - b) * e)
  p[7] <- (a * c * f) / ((a * c * f) + (b * d * g) + (e * h * (1 - a - b)))
  p[8] <- (b * d * g) / ((b * d * g) + (e * h * (1 - a - b)))
  p[9] <- ((f * c * a) + (g * d * b) + h * e * (1 - a - b)) / 
                                                 ((c * a) + (d * b) + e * (1 - a - b))
  p[10] <- g
  p[11] <- w
  p[12] <- ((b * d) + (w * e) * (1 - a - b)) / ((b * d) + e * (1 - a - b))
 
# Estimated basic
# prior parameters:
 a ~ dbeta(1, 2)
 c ~ dbeta(1, 1)
 d ~ dbeta(1, 1)
 e ~ dbeta(1, 1)
 f ~ dbeta(1, 1)
 g ~ dbeta(1, 1)
 h ~ dbeta(1, 1)
 w ~ dbeta(1, 1)
 z ~ dbeta(1, 1)

# Estiamted functional
# prior parameter:
b <- z * (1 - a)

# Distribution for NB of Maternal Diagnosis:
 M <- 60012 - 54296 * Y

# Maternal diagnosis:
Y ~ dgamma(0.56, 3)T(0, 2)

# Net Benefit:
nb[1] <- 0
nb[2] <- 105000 * (1 - a - b) * (M * e * (1 - h) - 3.0 * (1 - e * h))

# EVPI:
nb.sorted <- sort(nb[])
vpi <- 7.7127 * (nb.sorted[N.k] - nb[k.current])

}
"
writeLines(text = model_String, con = "Chapter_12.txt")

jags_data <- list(
 r = c(11044, 12, 252, 10, 74, 254, 43, 4, 87, 12, 14, 5),
 n = c(104577, 882, 15428, 473, 136139, 102287, 60, 17, 254, 15, 118, 31),
 N.k = 2,
 k.current = c(1, 2)
)


params <- c("nb[2]", "vpi[2]", "vpi[1]")

jags_Mod <- jags(data = jags_data, model.file = "Chapter_12.txt", 
                 parameters.to.save = params, n.iter = 50000, 
                 n.burnin = 10000, n.chains = 4)
jags_Mod

# Posterior visual checks:
posterior <- as.array(jags_Mod$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("nb[2]"))

# ==========================================================================================
# Expected value of partial perfect information ---------------------------
# ==========================================================================================





