# mlesurv, Bias Correction of the ML Estimates in the Proportional Hazards Model with Frailty 
Alexander Begun
## Description
Fits Parametric Frailty Models with Proportional Hazards Functions by maximum marginal likelihood.
Provides parameter estimates, their standard errors, bias corrected estimates, times-to-failure, empirical, estimated and bias corrected survivals. 
Possible baseline hazards: exponential, Weibull, Gompertz.
Possible frailty distributions: gamma, inverse gaussian.
## Details
Maximum likelihood parameter estimates in the proportional hazards model with random effect suffer under bias. The Cox-Snell and Cordeiro-Klein methodology 
allows us to estimate this bias and correct the estimates.

The data are clustered. It is assumed that all subjects in a cluster share the same unobserved risk of failure (frailty) that is either
the gamma either the inverse gaussian distributed random variable. The exponential, Weibull and the Gompertz baseline hazard functions are implemented. 
Expectations needed for the Cox-Snell methodoligy are calculated using the Monte Carlo integration.
