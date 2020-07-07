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
# I will focus on this type of generalised evidence synthesis. I think that combining RCT and
# observational evidence can be overly complex and opaque. Borrowing strength via hierarhical
# modelling makes more sense, to me.

