# ========================================================================================
# Multi-parameter Evidence Synthesis --------------------------------------
# ========================================================================================
library(rjags)
library(R2jags)
library(bayesplot)
library(tidyverse)
library(parallel)

options(mc.cores = detectCores())
set.seed(20202204)

# ========================================================================================
# Prior and posterior simulation in a probabilistic model: Maple Syrup Urine Disease (MSUD)
# ========================================================================================
# Each of these basic parameters requires prior distributions, and as they are probabilities,
# Beta distributions are a natural choice. We choose minimally informative (uniform) 
# distributions, Beta(1,1). Next we form expressions for the functional parameters that can
# be used to assess the impact of screening on the numbers of MSUD-associated cases of MR. 
# Functional parameters can, by definition, be written as functions of basic parameters.

# Pr (mental retardation in MSUD with screening):
# theta[sm] = phi[s] * theta[em] + (1 - psi[s]) * theta[lm]

# Pr (mental retardation in MSUD without screening):
# theta[nm] = phi[n] * theta[em] + (1 - psi[n]) * theta[lm]

# Birth rate of MSUD-associated retardation with screening:
# theta[s] = 100, 000 * r * theta[sm]

# Birth rate of MSUD-associated retardation without screening:
# theta[n] = 100, 000 * r * theta[nm]

# Reduction in MSUD-associated retardation due to screening:
# e[d] = theta[s] - theta[n]

model_String <- "model {
  # Likelihood models:
   # for no. of 
   # MSUD cases:
   r.r ~ dbin(r, n.r)
   # no. detected early
   # (screening):
   r.s ~ dbin(phi.s, n.s)
   # no. deteced early
   # (no screening):
   r.n ~ dbin(phi.n, n.n)
   # no. retarded
   # (early detection):
   r.em ~ dbin(theta.em, n.em)
   # no retarded
   # (late detection):
   r.lm ~ dbin(theta.lm, n.lm)
 
 # Priors for basic parameters:
 r ~ dbeta(1, 1)
 phi.s ~ dbeta(1, 1)
 phi.n ~ dbeta(1, 1)
 theta.em ~ dbeta(1, 1)
 theta.lm ~ dbeta(1, 1)
 
 # Functional parameters:
 theta.sm <- phi.s * theta.em + (1 - phi.s) * theta.lm
 theta.nm <- phi.n * theta.em + (1 - phi.n) * theta.lm
 theta.s <- (theta.sm * r) * 100000
 theta.n <- (theta.nm * r) * 100000
 e.d <- theta.s - theta.n
}
"
writeLines(text = model_String, con = "MSUDmod.txt")


Data_jags <- list(n.r = 724262, r.r = 7, n.s = 276, r.s = 253, n.n = 18,
                  r.n = 8, n.em = 10, r.em = 2, n.lm = 10, r.lm = 10)
Data_jags

params <- c("e.d", "theta.sm", "theta.nm", "theta.s", "theta.n")

jags_Mod <- jags(data = Data_jags, model.file = "MSUDmod.txt", parameters.to.save = params, 
     n.chains = 2, n.iter = 30000, n.burnin = 10000)
jags_Mod

posterior <- as.array(jags_Mod$BUGSoutput$sims.array)

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = params, 
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = params, 
           size = 1, ) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = params)

color_scheme_set("pink")
mcmc_pairs(posterior, pars = params,
           off_diag_args = list(size = 1.5))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = params, lags = 10)

# Note that in this particular case all the data directly informs five basic parameters, 
# and these parameters are sufficient to define the target quantity e.d, the reduction in 
# retardation that would follow screening for MSUD. The simplicity of the problem lies in 
# the fact that there are no data directly informing the functional parameters. Further, 
# because the Beta priors are conjugate with the Binomial likelihood, we can obtain the 
# posterior distributions for the basic parameters in closed form. A Beta(a, b) prior 
# together with binomial data with numerator y and denominator n gives us a 
# Beta(a + y, b + n ! y) posterior. It is this that makes it possible to sample from the 
# Bayesian posterior distributions using standard forward MC simulation. In the next 
# example, however, forward MC simulation would not be sufficient and we must appeal to 
# more powerful methods for Bayesian updating, such as MCMC, which do not rely on conjugacy.

# ========================================================================================
# A model for Prenatal HIV testing ----------------------------------------
# ========================================================================================
# The following analysis is simplified and somewhat stylised, but it illustrates an 
# important characteristic of the data sources that must often be used in studies of 
# screening, namely that data may not be available on the key parameters, but may instead 
# be available on functions of the parameters. Also, and in contrast to the MSUD example, 
# there may be data on more functions of parameters than there are basic parameters, which 
# raises the possibility of conflict between sources of evidence, or – to put it another 
# way – a chance to validate parameter estimates against independent evidence.

# The critical parameters in an economic analysis of screening are those that relate to the 
# prevalence of the undetected condition. Similarly, in a comparison of universal and 
# targeted testing, the key parameters are those relating to the prevalence of undiagnosed 
# infection in the Low Risk group, as the higher risk groups will be covered by both 
# screening strategies. The crucial parameters are, therefore, e, the prevalence of HIV in 
# the Low Risk, and h, the proportion of Low Risk HIV infected who are already diagnosed. 
# This intuition is confirmed by the Net Benefit equation. If M is the net monetary benefit 
# of an early diagnosis, T is the cost of an HIV test, then the Incremental Net 
# Benefit (INB) of universal relative to targeted screening for the women is:

# INB = (1 - a - b) * (Me(1 - h) - T(1 - eh))

# If this example was like MSUD, an investigator would be able to look through the 
# literature and find data sources that informed each of the four parameters in the Net 
# Benefit function. These could then be represented by suitable Beta distributions. Samples
# would then be drawn repeatedly from these distributions and INB calculated on each cycle,
# and finally a distribution of INB would be obtained. The problem, however, is that there 
# are no data that directly inform the key parameters e and h, and it is difficult to 
# imagine that there ever could be.

# An introductory synthesis exercise --------------------------------------
# We set ourselves the task of writing BUGS code to estimate the five parameters a, b, c,
# d and e from the first six items of data. We begin by assigning priors to these five 
# basic parameters. For c, d and e we can straightforwardly assign vague uniform priors 
# Beta(1, 1). For a and b we have to add a constraint that (a + b ≤ 1) because the 
# proportions in each risk group must each be greater or equal to zero, but must always sum 
# to 1. This is a Dirichlet distribution. A convenient way to represent a Dirichlet 
# distribution in BUGS is to express it as a series of binomial distributions. For example, 
# the code

# a ~ dbeta(1, 2)
# z ~ dbeta(1, 1)
# b <- z * (1 - a)

# sets up three variables: a, z = (1 - a), and b = z * (1 - a), where a is a probability, 
# and b is a proportion of 1 - a, so that the required constraint is satisfied whatever 
# values a and b may take. Note also that the priors give a, b and (1 - a - b) equal unit 
# weight, which represents a vague uniform Dirichlet(1, 1, 1).

# The remaining basic parameters can be given uniform distributions:

# c ~ dbeta(1, 1)
# d ~ dbeta(1, 1)
# e ~ dbeta(1, 1)

# The next step is to define the relationships between the basic parameters and the 
# functional parameters that the data estimates directly. Here it is convenient to create a
# vector p[ ], which in effect monitors the fitted values:

# p[1] <- a
# p[2] <- b
# p[3] <- c
# p[4] <- d
# p[5] <- (b * d + (1 - a - b) * e / (1 - a))
# p[6] <- (a * c + b * d + (1 - a - b) * e)


# The final step is to specify the Binomial likelihood, looping over the six items, and then
# to supply the data, in the form of a vector r[ ] for the numerators and a vector n[ ] for
# the denominators.

# So, the declared model logic takes the form:

model_String <- "model {
  # Binomial likelihood:
  for (i in 1:6) {
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
 
# Estimated basic
# prior parameters:
a ~ dbeta(1, 2)
c ~ dbeta(1, 1)
d ~ dbeta(1, 1)
e ~ dbeta(1, 1)
z ~ dbeta(1, 1)
# Estiamted functional
# prior parameters:
b <- z * (1 - a)

}
"
writeLines(text = model_String, con = "SyntEx.txt")

Data_jags <- list(
 r = c(11044, 12, 252, 10, 74, 254),
 n = c(104577, 882, 15428, 473, 136139, 102287)
)
Data_jags

params <- c("a", "b", "c", "d", "e", "z")

jags_Mod <- jags(data = Data_jags, model.file = "SyntEx.txt", 
                 parameters.to.save = params, n.iter = 50000, 
                 n.burnin = 10000, n.chains = 2)
jags_Mod

posterior <- as.array(jags_Mod$BUGSoutput$sims.array)

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = params, 
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = params, 
           size = 1, ) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = params)

color_scheme_set("pink")
mcmc_pairs(posterior, pars = params,
           off_diag_args = list(size = 1.5))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = params, lags = 10)

# ========================================================================================
# Model criticism in multi-parameter models -------------------------------
# ========================================================================================
#  As we shall see, statistical investigation can reveal that inconsistency is present, but 
#  is unlikely by itself to tell us which data input is the source of that inconsistency.
#  A further issue to consider in advance of any analysis is the number of potential 
#  inconsistencies that there can be.

#  A natural approach is to compare the number of parameters (nine) with the number of 
#  functions of parameters on which data is available, and regard the difference (three) as 
#  a measure of the number of independent inconsistencies there could be in the data set. 
#  We might say that there are three degrees of freedom for inconsistency.

# Inconsitency in the HIV data --------------------------------------------
#  The formal way to examine this is by cross-validation. Taking as an example, data point 
#  4, we need to form a predictive distribution for the numerator r[4], based only on data 
#  points 1 – 3 and 5 – 12, that is excluding data point 4, and then compare this predictive
#  distribution with the original data.

#  To do this the priors and the model sections of the BUGS/JAGS code remain unaltered, but
#  the likelihood is rewritten to exclude r[4]. The data list also remains unaltered. The 
#  cross-validation method then requires us to create a ‘replicate’ binomial variable based
#  on the p[4] (d), now estimated without the benefit of data point 4, and n[4], and 
#  compare this to the original r[4].

model_String <- "model {
  # Binomial likelihood:
  for (i in 1:3) {
   r[i] ~ dbin(p[i], n[i])
  }
  for (i in 5:12) {
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
 p[9] <- ((a * c * d) + (b * d * g) + (e * h * (1 - a - b)) /
             ((a * c) + (b * d) + e * (1 - a - b)))
 p[10] <- g
 p[11] <- w
 p[12] <- ((b * d / ((b * d) + e * (1 - a - b)))) + (e * w) * (1 - a - b) / 
             ((b * d) + e * (1 - a - b))
 
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
# prior parameters:
b <- z * (1 - a)

r.rep ~ dbin(p[4], n[4])
p.xval <- step(r.rep - r[4]) - 0.5 * equals(r.rep, r[4])
dev[4] <- 0

}
"
writeLines(text = model_String, con = "ModCrit.txt")

Data_jags <- list(
 r = c(11044, 12, 252, 10, 74, 254, 43, 4, 87, 12, 14, 5),
 n = c(104577, 882, 15428, 473, 136139, 102287, 60, 17, 254, 15, 118, 31)
)
Data_jags

params <- c("p", "p.xval", "r.rep")

jags_Mod <- jags(data = Data_jags, model.file = "ModCrit.txt", 
                 parameters.to.save = params, n.iter = 50000, 
                 n.burnin = 10000, n.chains = 2)
jags_Mod
# The posterior for p.xval and hence the probability of observing a sample r[4] as high as 
# 10 with a sample of n[4], given our model and the remaining data is 0.004.

posterior <- as.array(jags_Mod$BUGSoutput$sims.array)

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("r.rep", "p[4]"), 
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("r.rep", "p[4]"),
           size = 1, ) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("p[4]", "r.rep"))

color_scheme_set("pink")
mcmc_hist(x = posterior, pars = c("r.rep"), binwidth = 1)

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("r.rep"), lags = 10)

# ========================================================================================
# Summary of key points ---------------------------------------------------
# ========================================================================================
# We distinguish between forward Monte Carlo simulation from distributions of parameters, 
# and Monte Carlo sampling from the joint posterior distribution of parameters.

# Multi-parameter evidence synthesis combines evidence not just from different sources, but 
# on different parameters, within a coherent overall model. This may often involve more 
# sources of evidence, or more types of sources, than there are parameters, introducing the 
# possibility of inconsistency between different evidence sources.

# But, equally, the possibility of inconsistency is also an opportunity for independent 
# validation of a model.

# We emphasise the advantages of assessing the potential sources of bias in evidence before 
# attempting synthesis, due to the dangers of post hoc adjustments.

# Bayesian multi-parameter evidence synthesis meets the requirements we should expect from 
# methods under-pinning ‘evidence-based’ policy: incorpo- ration of all available evidence, 
# assessment of inconsistency, and uncertainty propagation.

# The danger with post hoc adjustments intended to resolve inconsistency is that the same 
# inconsistency can be resolved equally well by removing any one of a number of data points, 
# or adding bias parameters in several different parts of the network. Each different 
# approach constitutes an alternative model with different parameter values, and potentially
# with fundamentally different public health impli- cations, and yet there is no statistical
# way of choosing between them. The entire analysis becomes prey to selective 
# interpretation. Therefore, it is clearly preferable to identify potential sources of bias 
# attaching to each data point prior to putting data together in a synthesis.

# Thus! the power and convenience of the method must not be allowed to lure the user away 
# from a careful consideration of the assumptions being made, nor from careful attention to 
# convergence.

# End file ----------------------------------------------------------------