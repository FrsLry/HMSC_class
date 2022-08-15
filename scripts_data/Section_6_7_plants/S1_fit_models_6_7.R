# THIS SCRIPT CONSTRUCTS AND FITS HMSC MODELS FOR THE PLANT EXAMPLE (SECTION 6.7) OF THE BOOK
# Ovaskainen, O. and Abrego, N. 2020. Joint Species Distribution Modelling - With Applications in R. Cambridge University Press.

# We first set the working directory to the directory that has all the
# files for this case study. We recommend that you create always the
# subdirectories "data" and "models" to the main directory. If running
# the script in another computer, all that needs to be changed in the
# script below is the main directory. Working directory below works
# only in our local computer system, and you need to modify this in
# your own system. For instance, in a MacOS system this could be
# "~/Data/HMSC_book/Section_6_7_plants".

## see your current working directory first
getwd()

# replace the lines below so that they correspond your working directory
# wd = "C:/HMSC_book/R-scripts/Section_6_7_plants"
# setwd(wd)

localDir = getwd()
data.directory = file.path(localDir, "scripts_data/Section_6_7_plants/data/")
model.directory = file.path(localDir, "scripts_data/Section_6_7_plants/models/")

# We load the Hmsc package, and set the random seed, so that all
# results from this script are reproducible.

library(Hmsc)
set.seed(1)

# The factor variables are coded as character strings in our data
# set, and we need them be handled as factors, but this depends on
# the version of R and local settings of options, and therefore we
# use 'stringsAsFactors=TRUE' in read.csv
# We also set the numerical variable site as a factor

data = read.csv(file=file.path(data.directory,"whittaker revisit data.csv"),
                stringsAsFactors=TRUE)
data$site = factor(data$site)
head(data)

# The data are in the long format: each row corresponds to one plant
# species in one site.  The column `value` is plant abundance.  The
# column `env` is Whittaker's index describing the site position on
# the topographic moisture gradient (TMG).

# Sites on mesic, north-facing slopes receive lower numbers of TMG
# and sites on warmer, south-facing slopes receive higher numbers of
# TMG (Damschen et al. 2010, Miller et al. 2018).

# The column `trait` is the functional trait that Miller et
# al. (2018) selected for their analyses: leaf tissue
# carbon-to-nitrogen ratio (C:N).

# The C:N ratio can be considered as a surrogate of competitive
# ability, plants with low C:N growing faster but having lower stress
# tolerance than plants with high C:N (Miller et al. 2018, Cornelissen
# et al. 2003, Poorter & Bongers 2006).

# It can be thus expected that species occurring on drier and warmer
# sites have on average higher C:N ratio, leading to positive
# relationship between TMG and C:N ratio.

# Miller et al. (2018) applied several statistical methods to examine
# if there is an association with the C:N ratio and the environmental
# gradient.

# Here we re-analyse this question with HMSC.

# We first reformat the data so that it works as input for Hmsc.To do
# so, in the script below we construct

# the matrix `Y` of species abundances

# the dataframe `XData` of the environmental variable TMG

# and the dataframe `TrData` of the trait C:N ratio.

sites = levels(data$site)
species = levels(data$species)
n = length(sites)
ns = length(species)

Y = matrix(NA, nrow = n, ncol = ns)
env = rep(NA,n)
trait = rep(NA,ns)
for (i in 1:n){
  for (j in 1:ns){
    row = data$site==sites[i] & data$species==species[j]
    Y[i,j] = data[row,]$value
    env[i] = data[row,]$env
    trait[j] = data[row,]$trait
  }
}

colnames(Y) = species
XData = data.frame(TMG = env)
TrData = data.frame(CN = trait)

head(Y)
head(XData)
head(TrData)

# Let us explore the raw data by plotting histograms of species
# prevalences and abundances.

# We define prevalence as the  fraction of occupied sampling units

P = colMeans(Y > 0)
hist(P, xlim=c(0,1), xlab = "Prevalence (P)")

# As usual for community ecology data, most species are rare in the
# sense of being present only in the minority of study sites.

# We define abundance as the mean number of individuals over sites
# where the species is present.

A = colSums(Y)/colSums(Y>0)
hist(A, xlab = "Abundance (A)")

# To account for relatedness among the species, we use a taxonomy as a
# proxy for phylogeny.

# To do so, we read in a classification of the species into families
# and genera.

# We then use the function `as.phylo` from the `ape` package to
# construct a taxonomical tree.

# We assume equal branch lengths among families, among genera within a
# family, and among species within a genus.
library(dplyr)
taxonomy = read.csv(file=file.path(data.directory, "taxonomy.csv"))
taxonomy <- taxonomy %>% mutate_all(as.factor)
library(ape)
plant.tree <- as.phylo(~family/genus/species, data=taxonomy, collapse = FALSE)
plant.tree$edge.length <- rep(1,length(plant.tree$edge))
plot(plant.tree)

# After this initial exploration of the data, we are ready to set up
# the HMSC models.

# To examine the robustness of the results, we fit two models:

# a probit model for presence-absence data, and a lognormal Poisson
# model for the full count data.

# We first define the presence-absence model

# We will use TMG as the only environmental covariate, and assume a
# linear effect

XFormula = ~TMG

# We will use CN as the only trait covariate, and assume a linear
# effect

TrFormula = ~CN

# We next define the presence-absence model.

# We convert the count data to presence-absences by Y=(Y>0)

# We include environmental covariates through XData and XFormula

# We include trait covariates by TrData and TrFormula

# We include the phylogenetic relationships by phyloTree

# We assume the probit distribution as appropriate for
# presence-absence data

# The data have no specific structure and do not include any random
# effects, thus we do not define studyDesign nor ranLevels

model.pa = Hmsc(Y=(Y>0),
                XData = XData,  XFormula=XFormula,
                TrData = TrData, TrFormula = TrFormula,
                phyloTree = plant.tree,
                distr="probit")

# It is always a good idea to look at the model object.

model.pa

# Evaluating the above line should give "Hmsc object with 52 sampling
# units, 75 species, 2 covariates, 2 traits and 0 random levels"

getCall(model.pa)

# Evaluating the above line shows the call by which the model object was constructed

# It is alway a good idea to look at also the X-matrix that Hmsc has
# generated from XData and XFormula

# This will show e.g. if continuous variables are accidentally treated
# as factors or vice versa

head(model.pa$X)

# Evaluating the above line shows that there are two columns in the X
# matrix: the intercept and the linear effect of TMG.

# Note that the Hmsc object says that there are two covariates, thus
# here the intercept is counted as well.

# It is alway a good idea to look at also the Tr-matrix that Hmsc has
# generated from TrData and TrFormula

# This will show e.g. if continuous variables are accidentally treated
# as factors or vice versa

head(model.pa$Tr)

# Evaluating the above line shows that there are two columns in the Tr
# matrix: the intercept and the linear effect of CN.

# Note that the Hmsc object says that there are two traits, thus here
# the intercept is counted as well.

# It is a good idea to explore also the other objects that are
# included in the model object, to ensure that they are as they
# should, and to gain a better understanding of the model structure

# As one example, let us plot the correlation strcuture (matrix C)
# derived from the phylogenetic tree

library(Matrix)
image(Matrix(model.pa$C))

# The abundance model is defined exactly in the same way,

# except that now the data are not trunccated to presence-absence,

# and now we assume a lognormal Poisson model rather than the probit model.

model.abu = Hmsc(Y=Y,
                 XData = XData,  XFormula=XFormula,
                 TrData = TrData, TrFormula = TrFormula,
                 phyloTree = plant.tree,
                 distr="lognormal poisson")

# Alternatively, we may define model.abu as an updated version of model.pa
model.abu = update(model.pa, Y = Y, distr="lognormal poisson")

# both will give the same model version
getCall(model.abu)

# We will fit both models so that we store 100 posterior samples for
# each of two chains
# We will set these below as nChains = 2, samples = 100

# We note that for more "final" results, one might wish to have
# e.g. 1000 samples for each of four chains

# In the script below, we fit both models (presence-absence and
# abundance) one after each other

# After fitting the model, we save the ftted model object to a file

# Loading the fitted model object then serves as the starting point
# for exploring the results

# The script runs over a loop where thin is first 1, then 10, then
# 100, and so on

# Thin is the thinning parameter of the MCMC chain.

# The transient (also called burn-in) is set to 50*thin

# When thin = 1, there will be 50 burn-in and 100 actual
# iterations. All actual iterations are stored.

# When thin = 10, there will be 500 burn-in and 1000 actual
# iterations. The actual iterations are thinned by 10, so 100 are
# stored.

# When thin = 100, there will be 5000 burn-in and 10000 actual
# iterations. The actual iterations are thinned by 100, so 100 are
# stored.

# A long MCMC chain is needed to achieve convergence

# Thinning is applied to avoid storing model objects of very large
# size

# Even if in the end thin = 1000 is required to achieve converge, We
# recommend to run the loop thin = 1, 10, 100, 1000

# This is for several reasons.

# First of all, it will not be known beforehand how much thinning is
# needed to achieve satisfactory convergence

# Second, thin = 1 will run very fast, whereas thin = 1000 will take
# very long (1000 times longer)

# After thin = 1 is completed, it is already possible to develop all
# the remaining scripts that explore the fitted model

# When exploring the fitted model, often one realizes changes that
# need to be made, even if the fitting has not converged

# Third, running the model fitting for thin = 1, 10, 100, 1000 does
# not take much longer than running it just for thin = 1000 (it takes
# ca. 12% longer)

# Thus, in summary, running the model fitting for thin = 1, 10, 100,
# 1000 typically saves a lot of time, as it allows one to proceed fast
# in writing (and revising) all the scripts that are needed from
# defining the model to producing the result tables and figures

# The idea is not to run the entire loop in one go, as that would take
# a lot of time. Just run thin = 1, and then move to develop the next
# scripts.

# Most modern computers have several CPUs and you can run chains in
# parallel to further save time. Below we set nParallel=2 so that we
# use one  CPU for each of the two chains.
# Setting nParallel > 1 has two consequences: tracing info
# vanishes, and random sequences will change from nParallel=1 and
# results are no longer reproducible compared to that choice

# You may then leave the remaining part of the loop (e.g. thin = 10,
# 100, 1000) to run e.g. overnight

## We set up basic model parameters

samples = 100
nChains = 2
nParallel = 2

# for (thin in c(1,10,100,1000)){
for (thin in c(1,10)) {
  transient = 50*thin
  model.pa = sampleMcmc(model.pa, thin = thin, samples = samples, transient = transient, nChains = nChains, nParallel = nParallel)
  filename = file.path(model.directory, paste0("model_pa_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
  save(model.pa, file=filename)
  model.abu = sampleMcmc(model.abu, thin = thin, samples = samples, transient = transient, nChains = nChains, nParallel = nParallel)
  filename=file.path(model.directory, paste0("model_abu_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
  save(model.abu, file=filename)
}

# The script saves the finalized model after each thin. This means
# that you will not lose finished models if you interrupt the script,
# but they are in directory models.
