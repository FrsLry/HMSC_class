# THIS SCRIPT EXPLORES THE PARAMETER ESTIMATES OF HMSC MODELS FOR THE BIRD EXAMPLE (SECTION 11.1) OF THE BOOK
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

# We first perform a variance partitioning of the model.
# To be able to group the environmental variables, we look at the design matrix X that Hmsc has constructed
# by applying the XFormula to the XData.

head(m$X)

# We observe that the columns 2-5 relate to habitat variation and columns 6-7 to climatic variation.
# Arbitrarily, we include the intercept in the habitat variables, and thus perform the variance partitioning
# with the following grouping of the columns of the X matrix.

groupnames = c("habitat","climate")
group = c(1,1,1,1,1,2,2)
VP = computeVariancePartitioning(m, group = group, groupnames = groupnames)
plotVariancePartitioning(m,VP)
# Note that alternatively you may set 'las=2' to draw species titles vertically 
# But then you need to set margins wide enough to hold them
op = par(mar=c(9,4,1,1)+.1)
plotVariancePartitioning(m,VP, las=2)
par(op)

# The results of the variance partitioning show that the model (which is Model FULL of the book)
# relies only little on the spatial random effect.
# The climatic variable (including the squared term) explains roughly twice much variance as the
# habitat variable. We note that the numbers in the box in the figure show average variance proportions
# over the species, and that for some individual species (shown by the bars)
# the habitat matters more than climatic conditions. To examine variation among the species, we visualise their niches in Figure XX, which plot we generated for Model FULL with the `plotBeta` function.

# We next construct a beta-plot showing the estimates of species niche parameters

postBeta = getPostEstimate(m, parName="Beta")
plotBeta(m, post=postBeta, plotTree = TRUE, spNamesNumbers=c(FALSE, FALSE))

# In this figure, the habitat preferences of the species are compared to the reference level
# of broadleaved forests (Br). We observe predominatly negative responses to the habitat type Wetland (We),
# meaning that most species prefer broadleaved forests over wetland. Concerning the other
# habitat types (coniferous forests, open habitats and urban habitats), there are several
# species that show either a higher or a lower preference than for broadleaved forests.
# For most species, the linear response to spring temperature is primarily positive whereas
# the quadratic effect is primarily negative. As we shall she in the prediction plots,
# this means that the occurrence probability of most species peaks at warm climatic conditions.

# Negative quadratic to temperature actually means that the response
# to 'climate' has a peak at the 'climate' gradient (response is
# unimodal). Positive quadratic effect would indicate that the
# response is bimodal and has a minimum at 'climate' increasing to
# both directions. In both cases this only is true if the extreme
# point (peak, valley) is within the studied range of 'climate'. If
# the extreme is outside the studied 'climate' range, then we have a
# non-linear increasing or decreasing response over the studied
# range. If you wish to get an estimate of where the peak is, a
# rough estimate can be found using posterior mean betas from 'postBeta'.
# Note that to do this more properly, one should evaluate the full posterior distribution 
# of the estimate of the peak; see e.g. text around Fig. 10.16 of the book for an 
# example of how that can be done.

# First look at the range of 'clim':
range(m$XData$clim)

# The location of the extreme point (peak, valley) is defined by the
# first degree and second degree polynomial coefficients beta6 and
# beta7 as -beta6/2/beta7

-postBeta$mean[6,]/2/postBeta$mean[7,]

# For instance, Corvus_monedula is outside the range of clim and
# increases towards south, but Phylloscopus_trochilus peaks somewhere
# in central Finland. The possible
# location of the peak/valley and the increase towards S/N can only
# be deduced by the coefficients together and meaningfully only
# within the studied range of environment.

# We examine next if the species niches are linked to their traits with a Gamma-plot

postGamma = getPostEstimate(m, parName="Gamma")
plotGamma(m, post=postGamma, supportLevel = 0.95)

# Assuming that the results look the same as in the book (which may not be the case if
# e.g. sufficient mixing is not achieved), the plot suggests that the evironmental covariates
# are not tightly linked to the trait data.

# Let us next relax the level of required statistical support from 0.95 to 0.85.

plotGamma(m, post=postGamma, supportLevel = 0.85)

# Now resident migratory behaviour appears to be positively correlated with coniferous
# habitat use and short migratory behavior negatively correlated with wetland use.

# Another way of examing the influence of traits is to see how much of the variation they
# explain among the responses of the species to their covariates.

VP$R2T$Beta

# These results are consistent with the above figures: the traits explain only a very minor
# part of the variation. The same negative result is obtained also in the sense that
# traits explain only a negligible proportion of variation in species occurrence:

VP$R2T$Y

# We next evaluate the posterior distribution of the phylogenetic signal in species niches

mpost = convertToCodaObject(m)
round(summary(mpost$Rho, quantiles = c(0.025, 0.5, 0.975))[[2]],2)

# There is no evidence of phylogenetic signal in species niches

# In summary, the results reported above indicate that the 50 species included in our
# analyses respond rather individualistically to environmental variation,
# without a strong bearing to the included traits and phylogenetic relatioships.

# We next illustrate the species associations revealed by the random effects with the corrplot function.

library(corrplot)
OmegaCor = computeAssociations(m)
supportLevel = 0.95
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
toPlot = sign(toPlot)
plotOrder = corrMatOrder(OmegaCor[[1]]$mean,order="AOE")
corrplot(toPlot[plotOrder,plotOrder], method = "color", tl.cex=0.5,
         col=colorRampPalette(c("blue", "white", "red"))(255))

# The red and blue colours indicate those species pairs for which the support for
# either a positive or negative association is at least 0.95.
# Note that we have ordered the species (with the function corrMatOrder) so that the cluster of
# co-occurring species is most easily seen. To keep the original species order, write e.g. plotOrder = 1:m$ns

# We next examine at which spatial scale the variation captured by the random effect occurs.

mpost = convertToCodaObject(m) # mpost was calculated above
summary(mpost$Alpha[[1]], quantiles = c(0.025, 0.5, 0.975))

# There is no support for spatial signal in the residual, as the posterior distribution
# of the alphas overlap with zero. Thus, the variation is independent among the survey routes.
# Note that in case of the model without environmental covariates (not shown here but see the book),
# the variation in the leading factor occurs at the scale of ca. 150 km, reflecting the scale
# at which the relevant environmental conditions vary.
