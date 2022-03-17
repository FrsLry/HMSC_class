# THIS SCRIPT CONSTRUCTS AND FITS HMSC MODELS FOR THE FUNGAL EXAMPLE (SECTION 7.9) OF THE BOOK
# Ovaskainen, O. and Abrego, N. 2020. Joint Species Distribution Modelling - With Applications in R. Cambridge University Press.

# We first set the working directory to the directory that has all the files for this case study
# We recommend that you create always the subdirectories "data" and "models" to the main directory
# If running the script in another computer, all that needs to be changed in the script below is the main directory

# wd = "C:/HMSC_book/R-scripts/Section_7_9_fungi"
# setwd(wd)

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")

# We load the Hmsc package, and set the random seed, so that all results from this script are reproducable

library(Hmsc)
set.seed(1)

# We next read in the data and have a look at it
data = read.csv(file=file.path(data.directory, "data.csv"),stringsAsFactors=TRUE)
n = dim(data)[1]
head(data[,1:6])

# For each of n=99 logs, the decay class (DC) of the log is classified as 1-4,
# with 1 represening freshly fallen logs and 4 much decayed logs.
# The column readcount is the total number of sequences obtained for each sample,
# which number can be viewed to represent observation effort and is commonly called the sequencing depth.
# The remaing columns include the sequence counts for each species.
# The sequences were identified (in a data pre-processing step) using the probabilistic taxonomic placement
# algorithm (PROTAX-Fungi, for reference see the book section).
# The sequences that could not be assigned with at least 50% probability to any species are pooled in the category `unk`.
# This category contains typically the majority of the samples. This is because species-level identification
# is challenging e.g. due to incompleteness of the reference databases.
# The remaining columns consist of those 413 fungal species that were detected in the data at least once with 
# at least 50% identification probability. Using 50% as a threshold means that many of the species-level names may be wrong,
# in which the correct species name likely to belong to some closely related species. 

# We first organise the data. We construct the dataframe XData that includes the decay class
# and the total sequence count, as these will be used as fixed effects.
# We then construct the dataframe YData of the species data, including the unknown class.

XData = data.frame(DC = as.factor(data$DC), readcount = data$readcount)
YData = data[,4:dim(data)[2]]

# We are primarily interested in co-occurrences among the species. Detecting co-occurrences requires a substantial
# amount of information, for which reason the above script selects only those species that were found in at least 10 logs.
# This leaves us with 31 species, as you can verify by checking dim(YData) before and after evaluating the script below

sel.sp = colSums(YData>0)>=10
YData = YData[,sel.sp]

# We next explore the raw data by plotting histograms of species prevalence (P) and abundance (A).
# While in the model we will use the raw sequence counts, here we illustrate the data as relative abundance,
# and thus below we normalize abundance by "sum(YData)". Further, the relative abundances
# involve many orders of magnitude, some being e.g. 0.0001 and other as 0.1. For this reason,
# below we take the log-10 transform of the relative abundances when making the histogram.
# This "spreads" the variation so that it more visible. The value of -4 corresponds to relative
# abundance 0.0001, the value of -3 to relative abundance 0.001, the value of -2 to relative abundance 0.01,
# the value of -1 to relative abundance 0.1, and the value of 0 to relative abudance of 1 (which would mean
# that the species dominates the entire sample.

P = colMeans(YData>0)
A = colSums(YData)/sum(YData)

par(mfrow=c(1,2))
hist(P,xlim=c(0,1),breaks = seq(from=0,to=1,by=0.1), col = "grey", xlab = "Prevalence")
hist(log(A,base=10),breaks=10, col = "grey", xlab = "Abundance")

# As usually, there is great variation among the species, some being common and others rare.
# By exploring the variables P and A, you can see e.g. that excluding the unknown class that is present in all sampling units,
# the most prevalent species are  Fomitopsis pinicola (prevalence=0.54), Capronia kleinmondensis (0.42) and Phellophilus nigrolimitatus (0.37).
# In terms of abundance, the unknown class represents the vast majority (86%) of all the 
# sequences included in our analyses. The most abundant species are F. pinicola (5%), P. nigrolimitatus (1.8%) and A. serialis (1.5%). 

# We will fit six different HMSC models to these data, consisting of three model types and two choices of the
# explanatory variables. Below, and in the book, all of these models are fitted in a single loop.
# However, it is not recommended to construct multiple models straight away in a loop, as that comes with
# the risk of not checking each of the models carefully. For example, with these data the simplest starting point
# would be to fit a presence-absence model with decay class and log-transformed read-count as explanator variables,
# and the random effect of the log as a random effect. This model can be defined as follows:

Y = as.matrix(YData)

# For this model We truncate the data to presence-absence
Y = Y>0

# Even if the model definition HMSC accepts logical values (TRUE/FALSE) as the response matrix Y for presence-absence data,
# some of the post-processing functions assume that Y has numerical 0/1 values. So it is safer to 
# convert Y to numerical either as Y = 1*Y or as

Y = apply(Y,MARGIN = 2,FUN = as.numeric)

XFormula = ~DC + log(readcount)

studyDesign = data.frame(sample = data$LogID)
rL = HmscRandomLevel(units = studyDesign$sample)

m = Hmsc(Y = Y.pa, XData = XData, XFormula = XFormula, 
         studyDesign = studyDesign, ranLevels = list(sample = rL),
         distr="probit")

# To understand what was done with the above script, let us look at the objects:
# The study design simply lists the id of each log.

head(studyDesign)

# The motivation for including the studyDesign when defining the model object m was that we wished
# to include a random effect at the log level. The random effects are included by the ranLevels option of Hmsc(...)
# The list includes all levels at which a random effect is to be included. Here the random effect is included as
# "sample = rL". The left-hand side (sample) corresponds to the relevant column of the studyDesign, for which
# the log-id column (the only column) was named as sample. The right-hand side (rL) describes the structure of 
# random effect. To see the structure, evaluate

rL

# This gives the information that "Hmsc random level object with 99 units. Spatial dimensionality is 0 and number of covariates is 0."
# Thus, rL is an "unstructured" random effect, with no reference to space or covariates. This is the simplest type of a random effect.
# There are 99 units, as each log is a separate unit. We note in passing that often the random effects are not defined at the sampling unit level,
# but at a higher hierarchical level that includes many sampling units. Let us assume e.g. that the studyDesign
# would include another column describing the plot to which the log would belong, and that we would like to set up
# another random effect at the level of the plot. In such a case, this should NOT be defined as
# rL.plot = HmscRandomLevel(units = studyDesign$plot), but rather as 
# rL.plot = HmscRandomLevel(units = levels(studyDesign$plot)). The levels() picks each level (here each plot)
# only once, and thus they should not be repeated when defining the random effect strurcture.
# This was a side remark that is not relevant to this case study, see examples on hieararchical study designs for more information on this topic. 

# When defining the object m, we imported the data as Y = 1*(Y>0) to truncate it to presence-absence
# Corresponding to this, we assumed the probit link function

# It is always a good idea to look at how the XData and XFormula are translate to the design matrix X. 
# These can be seen as e.g. by evaluating the following

m$XFormula
head(m$XData)
head(m$X)

# We observed that DC is treated as factor (as it should), with DC=1 being the baseline level and thus not visible in the matrix X
# Further, log(readcount) is treated as continuous covariate (as it should)

# It is always a good idea to look at the model object as well:

m

# This should give "Hmsc object with 99 sampling units, 31 species, 5 covariates, 1 traits and 1 random levels"
# Note that "5 covariates" corresponds to the columns of the X matrix. We have two environmental variables,
# but with adding the intercept and expanind the factor of DC into three columns, there are in total 5 "covariates"
# Note further that there is "1 trait" even if we did not define any traits. To see why this is the case, evaluate

head(m$Tr)

# This shows that the trait matrix contains the intercept, which is the "1 trait". The intercept models the 
# mean response of the species to the environmental covariates.

# Above we exemplified how to construct a single model, and we recommend the users of Hmsc to typically define
# their models one by one. But when the collection of the models is "ready", it is often convenient
# to treat all of them within a single loop. In this way, e.g. the code for model fitting or evaluating
# MCMC convergence needs to be written only once.
# Thus, i the scripts below, the object models will include all six models that were discussed in the book chapter, organised as a
# list of lists, so the models[[i]][[j]] is the model type i=1,2,3 for the choice of explanatory variables j=1,2.

# The model type i=1 is a lognormal Poisson model that is fitted to the sequence count data as they are.
# So this model will simultaneously account both for whether the species is present or not as well as how abundant
# it is when it is present. The model types i=2 and i=3 form together a hurdle model that separates presence-absence
# variation from abundance variation. Thus, model i=2 is a probit model that is fitted to sequences counts
# truncated to presence-absence, whereas model i=3 is a normal model that is fitted to log-transformed sequence
# counts conditional on presence. This means that for i=3, sampling units where the species is not present are shown in the
# Y matrix as missing data (NA) rather than as zeros.

# Concerning the choices on the explanatory variables, we always include log-transformed sequencing depth,
# as this variable measures the total observation effort. Sequencing depth is the only variable included with j=1,
# and thus in these models we will estimate raw co-occurrences. With j=2, we additionally include the categorical
# variable of the decay class, and thus in these models we will estimate residual co-occurrences.
# We note that there are also many other properties of the log than the decay class that are likely to
# influence fungal occurrences, such as the diameter of the log. However, we ignore the other variables,
# as the data was specifically sampled to include much variation in decay class but less variation in other properties such as diameter.

# In all models, we also a random effect at the sampling unit level. The random effect models associations among the species,
# which is what we are primarily interested about.

models = list()
for (i in 1:3){
  Y = as.matrix(YData)
  if (i==2) {Y = 1*(Y>0)}
  if (i==3) {
    Y[Y==0] = NA
    Y = log(Y)
  }
  tmp = list()
  for (j in 1:2){
    XFormula = switch(j, ~1 + log(readcount), ~DC + log(readcount))
    m = Hmsc(Y = Y, XData = XData, XFormula = XFormula, 
             studyDesign = studyDesign, ranLevels = list(sample = rL),
             distr=switch(i,"lognormal poisson","probit","normal"),
             YScale = TRUE)
    tmp[[j]] = m
  }
  models[[i]] = tmp
}

# In the above script, we have used the option Yscale = TRUE to scale the response data to zero mean and unit variance.
# is discussed in more detail in Section 8.3 of the book, this scaling influences only the normal model,
# and it is done to make the default priors of Hmsc compatible with the data.

# We will fit each of the six models so that we store 100 posterior samples for each of two chains
# We note that for more "final" results, one might wish to have e.g. 1000 samples for each of four chains

nChains = 2
samples = 100

# We next loop over both the model types as well as the selections of explanatory variables to fit all the six models.
# After fitting all models, we save the models object (including the six fitted model objects) to a file
# Loading the fitted models then serves as the starting point for exploring the results
# The script runs over a loop where thin is first 1, then 10, then 100, and so on
# Thin is the thinning parameter of the MCMC chain.
# The transient (also called burn-in) is set to 50*thin
# When thin = 1, there will be 50 burn-in and 100 actual iterations. All actual iterations are stored.
# When thin = 10, there will be 500 burn-in and 1000 actual iterations. The actual iterations are thinned by 10, so 100 are stored.
# When thin = 100, there will be 5000 burn-in and 10000 actual iterations. The actual iterations are thinned by 100, so 100 are stored.
# A long MCMC chain is needed to achieve convergence
# Thinning is applied to avoid storing model objects of very large size
# Even if in the end thin = 1000 is required to achieve converge, We recommend to run the loop thin = 1, 10, 100, 1000
# This is for several reasons.
# First of all, it will not be known beforehand how much thinning is needed to achieve satisfactory convergence
# Second, thin = 1 will run very fast, whereas thin = 1000 will take very long (1000 times longer)
# After thin = 1 is completed, it is already possible to develop all the remaining scripts that explore the fitted model
# When exploring the fitted model, often one realizes changes that need to be made, even if the fitting has not converged
# Third, running the model fitting for thin = 1, 10, 100, 1000 does not take much longer than running it just for thin = 1000 (it takes ca. 12% longer)
# Thus, in summary, running the model fitting for thin = 1, 10, 100, 1000 typically saves a lot of time,
# as it allows one to proceed fast in writing (and revising) all the scripts that are needed from defining the model to producing the result tables and figures
# The idea is not to run the entire loop in one go, as that would take a lot of time. Just run thin = 1, and then move to develop the next scripts. 
# You may then leave the remaining part of the loop (e.g. thin = 10, 100, 1000) to run e.g. overnight

for (thin in c(1,10,100,1000)){
  transient = 50*thin
  for (i in 1:3){
    for (j in 1:2){
      cat("model = ",i, ", modeltype = ",j,"\n",sep="")
      models[[i]][[j]] = sampleMcmc(models[[i]][[j]], thin = thin, samples = samples, transient = transient,
                                    nChains = nChains, nParallel = nChains, initPar = if(i==3) {NULL} else {"fixed effects"})
    }
  }
  filename=file.path(model.directory, paste0("models_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
  save(models,file=filename)
}

# MCMC convergence can be difficult to achieve especially in those models that are not base on normal distribution
# For this reason, in the script above we initialize the "lognormal poisson" (i=1) and "probit" (i=2) models with
# initPar="fixed effects", with which option the MCMC chains are not started from locations randomized from the prior
# but from a maximum likelihood solution to the fixed-effects part of the model
