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

# A fixed treatment effects model for MTC ---------------------------------
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


















