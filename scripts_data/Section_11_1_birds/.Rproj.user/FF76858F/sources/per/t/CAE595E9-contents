# setwd("C:/LocalData/ovaskain/all stuff/TEACHING/HMSC/2020 November/R-demonstration-1")
rm(list = ls())
library(Hmsc)
library(abind)
library(ggplot2)
library(gridExtra)

set.seed(1)

#READING, SELECTING AND FORMATTING THE DATA
da = read.csv("bird data\\data.csv")
# Data only for 2014
da = droplevels(subset(da,Year==2014))
#♦Habitat as a factor
da$Habitat = as.factor(da$Habitat)
XData = data.frame(hab=da$Habitat, clim = da$AprMay)
Y = as.matrix(1*(da$Corvus_monedula>0))
colnames(Y) = "Corvus monedula"
xy = as.matrix(cbind(da$x,da$y))

# Take a look at data
head(XData)
head(Y)


# Original data of Y
(da$Corvus_monedula)
# You can see that they modify it to change this abundance into presence/absence data


#PLOTTING THE DATA
mapData=data.frame(x=xy[,1],y=xy[,2],O=Y[,1],H=XData$hab, C=XData$clim)
ggplot(data = mapData, aes(x=x, y=y))+geom_point(size=1.5,aes(color=H)) + ggtitle("Habitat") + scale_color_discrete(name="") + theme(legend.position="bottom")
ggplot(data = mapData, aes(x=x, y=y, color=C))+geom_point(size=1.5) + ggtitle("Climate") + scale_color_gradient(low="blue", high="red", name="") + theme(legend.position="bottom")+labs(y="")
ggplot(data = mapData, aes(x=x, y=y, color=O))+geom_point(size=1.5) +  ggtitle("Corvus monedula") + scale_color_gradient(low="blue", high="red", name="")+theme(legend.position="bottom")+labs(y="")


#######################################################################
#SETTING UP THE HMSC-MODEL (actually 3 different models: m1, m2, m3)

# the study design is the a df. You need the study design if you have random effect
studyDesign = data.frame(route = as.factor(da$Route))

# You had the studydesign as the rowname of the xy coordinates
rownames(xy)=studyDesign[,1]

# the random level for the spatial coordinates: spatial data is the given by the xy coordinates
# You can add the temporal coordinates you can add a 3 dimension
rL = HmscRandomLevel(sData=xy)
# now you have the smapling unit level random effect

# This is the fixed effect: linear effet of the habitat and the 2nd order polynomial effect
# of the climatic conditions
XFormula = ~ hab + poly(clim,degree = 2,raw = TRUE) # raw relates to how polynomial function scales the first and second order term f the polynomial
# For Hmsc put always raw = T.


#######################################################################
# Create the models
# You implement Y, X data, the XFormula. Model 1 and 3 have the random effect.
# M1 is full model: both XFormula = fixed level and then random level
m1 = Hmsc(Y=Y, XData = XData, XFormula=XFormula,
          distr="probit", studyDesign=studyDesign,
          ranLevels=list(route=rL))


# Take a look at m1
m1

# About the covariates:
# you have 7 columns which are the intercept + 6 covariates.
# You have the habitat types. You have the first and second order of the climatic
# conditions . The 1 is the mean of the temperature (climatic conditions itself).
# The 2nd is the just the (first order)².

# To access the random part:
m1$ranLevels



# No random effect for m2 = fixed effect only
m2 = Hmsc(Y=Y, XData = XData, XFormula=XFormula,
          distr="probit", studyDesign=studyDesign)


# Intercept only model, implemented in XFormula
# Random effect only model
m3 = Hmsc(Y=Y, XData = XData, XFormula=~1,
          distr="probit", studyDesign=studyDesign,
          ranLevels=list(route=rL))


#################################################
#FITTING THE MODEL
# The number of iteration per chain is computed like this: (thin*samples + transient)*nChain
# Set the MCMC parameters (those one are ridiculously low for testing):
nChains = 2
samples = 20
thin = 1 # means that you take only one of the iterations when they're too
transient = round(0.5*samples*thin)



# Using the MCMC method to fit the model for the 3 models:
m1 = sampleMcmc(m1, thin = thin, samples = samples,
           transient = transient, nChains = nChains, )

# m2 runs faster cause no random effect. If you don't have random effect,
# it is not necessary to specify the study design
m2 = sampleMcmc(m2, thin = thin, samples = samples,
                transient = transient, nChains = nChains)
m3 = sampleMcmc(m3, thin = thin, samples = samples,
                transient = transient, nChains = nChains)
models = list(m1,m2,m3)
filename =paste0("models/model_chains_",as.character(nChains),
                 "_samples_",as.character(samples),
                 "_thin_",as.character(thin))
save(models,file = filename)

#READING IN THE FITTED MODEL
nChains = 4
samples = 250
thin = 1
filename =paste0("models/model_chains_",as.character(nChains),
                 "_samples_",as.character(samples),
                 "_thin_",as.character(thin))
load(filename)
models

#EXAMINING MCMC CONVERGENCE
# let's take the model 1
m = models[[1]]
mpost = convertToCodaObject(m)

# The effective size must be close to the number of chain*sample
effectiveSize(mpost$Beta)

# Compute the Gelman diagnosis (when close to 1 it's good)
gelman.diag(mpost$Beta,multivariate=FALSE)$psrf


# The alpha parameter is the spatial sclae parameter (better to lookat)
effectiveSize(mpost$Alpha[[1]])
gelman.diag(mpost$Alpha[[1]],multivariate=FALSE)$psrf



################################################
#EVALUATING MODEL FIT
m = models[[1]]
preds = computePredictedValues(m)
## You evalute and obtain 3 values: RMSE (root mean square error),
# AUC (Area Under the Curve, good to evaluate how
# good a presence/absence model is) and TjurR² (other good way for the ppresence/absence)
# If AUC >.5 the model runs better than random
# AUC is a mesure of discrimination.
# If R² = 1 perfect model and 0 is like random
# TjurR² is lower than the classical AUC and R² so 0.36 is quite good
evaluateModelFit(hM=m, predY=preds)


###############################################
#COMPUTING VARIANCE PARTITIONING
m = models[[1]]
VP = computeVariancePartitioning(m)
# Gives how each parameters explains the variance of the data
VP$vals
# you can see that the spatial rnadom effect is almost not involved in the variance partitioning


#EXAMINING PARAMETER VALUES
m = models[[3]]
mpost = convertToCodaObject(m)

# Beta is for the covariates.
# For each covariates, you have either the positive or negative response to each
summary(mpost$Beta, quantiles = c(0.025, 0.5, 0.975))[[2]]

# Alpha is for the spatial random effect. It's the spatial scale parameter
# If the posterior median scale parameter is 0 km, it means that we don't have signal
# It gives the spatial scale at which the variation in the data is occuring is under 50%
summary(mpost$Alpha[[1]], quantiles = c(0.025, 0.5, 0.975))[[2]]



## The partition of variance is explained in part 5.5 of the book

#############################################################################
#PREDICTIONS OVER ENVIRONMENTAL GRADIENTS
m = models[[2]]

# For this, we need to create the gradient.
# The gradient object has 3 objects.
# XDataNew which contains the gradient where the gradient goes from the lowest of the data to
# the highest of the data.You can also the habitat type, that the data must know if you want to make prediction
# (this is the non.focalVarialbes arg). The number 1 indicates that we want that you want the main
# effect of the habitat (i.e. it will take the habitat the main present).
Gradient = constructGradient(m,focalVariable = "clim",
                             non.focalVariables = list(hab = 1))
Gradient$XDataNew


## Here, we predict how the presence of the species is constrained by the climate
## We predict the occurence probability.
# Expected = T means that you compute a mean of the prediction (posterior proba).
# If you put expected = F, you will only obtain 1 realisation which has
# a internal variation
predY = predict(m, Gradient=Gradient,expected = TRUE)
## It gives the posterior distribution for the 1000 samples
## We have prediction for each species for the 20 climatic conditions
## So it calculates the incertitude


## The points are the actual data . Black line is the posterior mean predictions
## The title of the    posterior proba that the prediciton of this plot for the
# case where the climatic predictor is at the smaller value (left of the axis)
# is smaller than the prediction on the right (clim = 8)
plotGradient(m, Gradient, pred=predY, measure="Y", index = 1, showData = TRUE)


## This time, gradient for the habitat type
Gradient = constructGradient(m,focalVariable = "hab",
                             non.focalVariables = list(clim = 1))
Gradient$XDataNew

predY = predict(m, Gradient=Gradient,expected = TRUE)

## you get the mean occurence proba and the 95% interval.
# You can see that the habitat Urb is more important and it means that you
# will find more probably Corvus in urban habitat
plotGradient(m, Gradient, pred=predY, measure="Y", index = 1, showData = TRUE)



#PREDICTIONS OVER SPACE
grid = read.csv("bird data\\grid_1000.csv")

## You have to remove marine habitat cause we hadn't it in the learning
grid = droplevels(subset(grid,!(Habitat=="Ma")))

## Take onle the xy coordinates
xy.grid = as.matrix(cbind(grid$x,grid$y))

## Create a df gathering the clim and hab column from the data we want to
# to predict on.
XData.grid = data.frame(hab=as.factor(grid$Habitat), clim=grid$AprMay)

m = models[[1]]
## You take your model and the newXdata. The sDataNew indicates the study design
Gradient = prepareGradient(m, XDataNew = XData.grid, sDataNew = list(route=xy.grid))

## Then predict.
predY = predict(m, Gradient=Gradient, predictEtaMean = TRUE, expected = TRUE)
length(predY)
length(predY[[1]])

# Now compute the posterior mean of occurence for the species COrvus
EpredY = apply(abind(predY,along=3),c(1,2),mean)
length(EpredY)

mapData=data.frame(x=xy.grid[,1],y=xy.grid[,2], EpredY,H=XData.grid$hab, C=XData.grid$clim)

## First predict the habitat type
ggplot(data = mapData, aes(x=x, y=y, color=H))+geom_point(size=1.5) +
  ggtitle("Habitat") + scale_color_discrete(name="") +
  theme(legend.position="bottom")

## Then the climatic values
ggplot(data = mapData, aes(x=x, y=y, color=C))+geom_point(size=1.5) +
  ggtitle("Climate") + scale_color_gradient(low="blue", high="red", name = "") +
  theme(legend.position="bottom") +labs(y="")

## And now the prediction
ggplot(data = mapData, aes(x=x, y=y, color=Corvus.monedula))+geom_point(size=1.5) +
  ggtitle("Corvus monedula")+ scale_color_gradient(low="blue", high="red", name = "") +
  theme(legend.position="bottom")+labs(y="")
