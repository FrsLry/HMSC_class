# THIS SCRIPT EVALUATES MCMC CONVERGENCE OF HMSC MODELS FOR THE BIRD EXAMPLE (SECTION 11.1) OF THE BOOK
# Ovaskainen, O. and Abrego, N. 2020. Joint Species Distribution Modelling - With Applications in R. Cambridge University Press.

# The preliminaries are as in script S1

# Change the working directory if necessary
# wd = "C:/HMSC_book/R-scripts/Section_11_1_birds"
# setwd(wd)

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")
library(Hmsc)
set.seed(1)

# The previous script saved the finished models saved the results in
# 'model.directory', and you can look at its content with
list.files(model.directory)

# We recommed running the script first with setting thin = 1
# Then rerun it with thin = 10, thin = 100, ... to examine how the convergence statistics improve

nChains = 2
samples = 100
thin = 1 # try with thin = 1, thin = 10, thin = 100, etc.
filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

# We restrict here the study of MCMC convergence to the examination of the potential scale reduction factor
# of the beta parameters and the Omega parameters. For the latter, we take a subsample of 200 randomly selected
# species pairs to avoid excessive computations.

mpost = convertToCodaObject(m, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
tmp = mpost$Omega[[1]]
z = ncol(tmp[[1]])
sel = sample(z, size=200)
# Here we take the subset of species pairs. We loop over the 2 MCMC chains.
for(i in 1:length(tmp)){ 
  tmp[[i]] = tmp[[i]][,sel]
}
psrf.omega = gelman.diag(tmp,multivariate=FALSE)$psrf

par(mfrow=c(1,2))
hist(psrf.beta, xlab = "psrf (beta)")
hist(psrf.omega, xlab = "psrf (Omega)")

# The MCMC convergence diagnostics can be considered satisfactory if for
# most parameters the potential scale reduction factor is close to the ideal value of one.
# For example, if most of the values are smaller than the value of 1.1. If you are
# evaluating this script with thin=1, the largest values you see are likely to be greater than ten.
# This means that MCMC convergence is not achieved, which means that the posterior sample
# can be a very bad approximation of the true posterior distribution. However, this should not
# stop one from making a preliminary exploration of the results, just keeping in mind
# that the results can be very misleading. Thus, even if the MCMC did not converge yet,
# you may proceed to the next scripts of exploring parameter estimates (S3) and making predictions (S4).
# Then re-run those scripts once you have a better approximatio of the posterior distribution.
