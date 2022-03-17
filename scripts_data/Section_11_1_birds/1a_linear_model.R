#LET US SIMULATE SOME DATA FROM THE LINEAR MODEL
set.seed(1)
n = 50
x = rnorm(n)
beta1 = 0 # intercept
beta2 = 1 # slope
L = beta1 + beta2*x
y = L + rnorm(n, sd = 1) # normally distribute error
plot(x,y)

#FITTING THE BASIC LINEAR MODEL to the data we just generated
basic.lm = lm(y ~ x) # basic linear model of R: the model is defined and fitted
abline(basic.lm)
summary(basic.lm) # you see that beta1 is
# approximated to 0.12 and beta2 to 0.96 while we know that their
# values are respectively 0 and 1


#FITTING THE LINEAR MODEL WITH Hmsc
library(Hmsc)
Y = as.matrix(y) # community matrix with 50 plots for 1 species
XData = data.frame(x = x) # env dataframe with 1 covariate
# Be careful to specify Y as a matrix and X as a dataframe !!!!!

# Construct the model using Hmsc. Will create a Hmsc object giving the number of
# species, covariates, traits and random levels
hmsc.lm = Hmsc(Y = Y, XData = XData, XFormula = ~x)
# You have a lot of different objects in it
hmsc.lm$XData # will give the data furnished
hmsc.lm$X # will add the intercept to the db furnished

# The Hmsc object says that we have a trait but actually it is just because
# you specify the intercept becausse the trait is implemented in the basic formula


# Now fit the model by adddin the posterior distribution to the previous obj
hmsc.lm = sampleMcmc(hmsc.lm, thin = 1, samples = 1000,
                     transient = 500, nChains = 2)

# Look what it added
hmsc.lm

# You need to first convert the Hmsc object to a coda object
mpost = convertToCodaObject(hmsc.lm)
# And plot it
plot(mpost$Beta) # give the trait plots of the 2 parameters
# You see the 1500 samplings. Each color correspond to a chain of the MCMC sampling
# Both colors seem to look the same, there's no trend in any of the chain
# Up: the intercept (C1)
# Down: covariate number 2 (C2) meaning the X specified



# Check the effective sample size. The actual sample size is the 2 chains * 1000 sample so
# 2000 actual sample size. If the effective sample size is low, it means that there is a very bad mixing
effectiveSize(mpost$Beta)


# Gelman diagnostic with $psrf (potential scale of reduction factor)
gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
# One number per parameter. Close to 1, the chains have converged: gooooood

# Look at the parameter estimates
summary(mpost$Beta) # best estimate for the estimate is .12. For the slope is .95. Then you have the SD
# The quantiles shouldn't be consider as similar to the confidence interval (but I doooooo)

# Compare them to the ones of the basic linear model
summary(basic.lm)


# Get the model fit
preds = computePredictedValues(hmsc.lm) # compute the predicted values
# you have 2000 posterior predictions

# You can use this function to resume the info
MF = evaluateModelFit(hM=hmsc.lm, predY=preds)

# You have R² (variance explained) and RMSE (root mean squared error)
MF$R2

#THIS IS HOW ONE COULD FIT E.G. A PROBIT MODEL
#m.probit = Hmsc(Y=Y, XData=XData, XFormula=~x, distr="probit")
