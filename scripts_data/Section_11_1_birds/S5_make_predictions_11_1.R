# THIS SCRIPT MAKES PREDICTIONS OVER ENVIRONMENTAL GRADIENTS AND OVER SPACE FOR HMSC MODELS FOR THE BIRD EXAMPLE (SECTION 11.1) OF THE BOOK
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

# We first read in the model object. You may change the thin parameter to the highest
# thin for which you have fitted the model, to get as reliable results as possible

nChains = 2
samples = 100
thin = 1
filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

# We start by making gradient plots that visualise how the communities vary among the 
# environmental variables.

Gradient = constructGradient(m,focalVariable = "clim")
predY = predict(m, Gradient=Gradient, expected = TRUE)

# Occurrence probability of Corvus monedula (to see other species, change index)

plotGradient(m, Gradient, pred=predY, measure="Y", index = 50, showData = TRUE)

# Occurrence probability of Numenius arquata with a peak on the
plotGradient(m, Gradient, pred=predY, measure="Y", index = 47, showData = TRUE)
# We can draw a vertical line for the estimated
# location of the peak from the posterior mean of Beta coefficients
postBeta = getPostEstimate(m, parName="Beta")
abline(v = -postBeta$mean[6,47]/2/postBeta$mean[7,47], col="red")

# Species richness

plotGradient(m, Gradient, pred=predY, measure="S", showData = TRUE)

# Proportion of species with resident migratory strategy

plotGradient(m, Gradient, pred=predY, measure="T", index = 2, showData = TRUE)

# Community mean weighted log-transformed body mass

plotGradient(m, Gradient, pred=predY, measure="T", index = 4, showData = TRUE)

# We observe that the occurrence probability of the species C. monedula increases with the climatic variable.
# which result is consistent with the single species analysis of Section 5.7 of the book.
# The species richness increases as well with increasing temperature, reflecting the
# fact that most species respond positively to temperature. Concerning the traits,
# we observe that the mean proportion of resident species increases with increasing
# temperature, whereas the average body mass decreases with it.

# We next plot the same features over the habitat gradient

Gradient = constructGradient(m,focalVariable = "hab")
predY = predict(m, Gradient=Gradient, expected = TRUE)
plotGradient(m, Gradient, pred=predY, measure="Y", index = 50, showData = TRUE,  jigger = 0.1)
plotGradient(m, Gradient, pred=predY, measure="S", showData = TRUE, jigger = 0.1)
plotGradient(m, Gradient, pred=predY, measure="T", index = 2, showData = TRUE, jigger = 0.1)
plotGradient(m, Gradient, pred=predY, measure="T", index = 4, showData = TRUE, jigger = 0.1)

# Consistently with the single species analysis of Section 5.7, the occurrence probability of the species
# C. monedula is the highest in urban habitats.. Species richness is lowest in open habitats and wetlands,
# and these habitats also have the lowest proportions of resident species.
# These results partially reflect the responses of the species to climatic conditions,
# as these two habitats are the most common in Northern Finland with the coldest climate.
# The average body size is essentially independent of habitat type.

# We next perform spatial predictions. To do so, we proceed as in the single-species case sudy in Chapter 4. Thus, we import a grid of spatial coordinates, habitat types and climatic conditions for 1000 locations. We then apply the `prepareGradient` function to these data to prepare a spatial gradient for which the predictions are to be made
# While in the book the predictions are presented for 10000 prediction locations, we recommend running this
# script for 1000 prediction locations to make the running time faster. To choose which one to do,
# read either the file "grid_1000.csv" to "grid_10000.csv"

grid = read.csv(file.path(data.directory, "grid_1000.csv"), stringsAsFactors=TRUE)

# Let's look at what is included in the grid

head(grid)

# We see that there  are the same environmental variables as in the data used to fit the model,
# as well as the spatial coordinates. This is what we need to make the predictions.

# The habitat type "Ma" is present in the prediction data but it is not part of training data,
# so we can't make predictions for it

grid = droplevels(subset(grid,!(Habitat=="Ma")))

# We next construct the objects xy.grid and XData.grid that have the coordinates and
# environmental predictors, named similarly (hab and clim) as for the original data matrix (see m$XData) 

xy.grid = as.matrix(cbind(grid$x,grid$y))
XData.grid = data.frame(hab=grid$Habitat, clim=grid$AprMay, stringsAsFactors = TRUE)

# We next use the prepareGradient function to convert the environmental and spatial
# predictors into a format that can be used as input for the predict function

Gradient = prepareGradient(m, XDataNew = XData.grid, sDataNew = list(Route=xy.grid))

# We are now ready to compute the posterior predictive distribution (takes a minute to compute it)
nParallel=2
predY = predict(m, Gradient=Gradient, expected = TRUE, nParallel=nParallel)

# Note that we used expected = TRUE to predict occurrence probabilities (e.g. 0.2) instead of occurrences (0 or 1) 
# Note also that if you have very large prediction grid, you can use the predictEtaMean = TRUE option to speed up the computations
# predY = predict(m, Gradient=Gradient, predictEtaMean = TRUE, expected = TRUE)

# Let's explore the prediction object.

class(predY)

# It is a list... 

length(predY)

# ...of length 200, if you fitted two chains with 100 samples from each

dim(predY[[1]])

# Each prediction is a matrix with dimensions 951 x 50, as there are 951 prediction locations and 50 species.

head(predY[[1]])

# Each matrix is filled in with occurrence probabilities
# We may simply by ignoring parameter uncertainty and just looking at 
# the posterior mean prediction. 

EpredY=Reduce("+",predY)/length(predY)
dim(EpredY)

# EpredY is a 951 x 50 matrix of posterior mean occurrence probabilities
# The next step is to post-process the predictions to those community features
# that we wish to illustrate over the prediction space. With the script below,
# we derive from the predictions the occurrence probability of C. monedula (species number 50),
# the species richness, and community-weighted mean traits.
# We also include data on habitat type and climatic conditions to the dataframe mapData that
# includes all the information we need to visualize the predictions as maps

Cm = EpredY[,50]
S=rowSums(EpredY)
CWM = (EpredY%*%m$Tr)/matrix(rep(S,m$nt),ncol=m$nt)
xy = grid[,1:2]
H = XData.grid$hab
C = grid$AprMay
mapData=data.frame(xy,C,S,Cm,CWM,H, stringsAsFactors=TRUE)

# We will use the ggplot function from the ggplot2 package, so let's load the data
library(ggplot2)

# We first plot variation in the habitat and climatic conditions on which the predictions are based on.

ggplot(data = mapData, aes(x=x, y=y, color=H))+geom_point(size=2) + ggtitle("Habitat") + scale_color_discrete() + coord_equal()
ggplot(data = mapData, aes(x=x, y=y, color=C))+geom_point(size=2) + ggtitle("Climate") + scale_color_gradient(low="blue", high="red") + coord_equal()

# We then exemplify prediction for one focal species, here C. monedula over Finland

ggplot(data = mapData, aes(x=x, y=y, color=Cm))+geom_point(size=2) + ggtitle(expression(italic("Corvus monedula")))+ scale_color_gradient(low="blue", high="red") + coord_equal()

# This prediction is reassuringly very similar to that based on the single-species model of Chapter 5.7
# We next plot predicted species richness, which is highest in Southern Finland

ggplot(data = mapData, aes(x=x, y=y, color=S))+geom_point(size=2) + ggtitle("Species richness")+ scale_color_gradient(low="blue", high="red") + coord_equal()

# We next plot the proportion of resident species, also highest in Southern Finland 

ggplot(data = mapData, aes(x=x, y=y, color=MigrationR))+geom_point(size=2) + ggtitle("Proportion of resident species")+ scale_color_gradient(low="blue", high="red") + coord_equal()

# We next plot the community-weighted mean log-transformed body size, which is highest in Northern Finland

ggplot(data = mapData, aes(x=x, y=y, color=LogMass))+geom_point(size=2) + ggtitle("Mean log-transformed body mass") + scale_color_gradient(low="blue", high="red") + coord_equal()
