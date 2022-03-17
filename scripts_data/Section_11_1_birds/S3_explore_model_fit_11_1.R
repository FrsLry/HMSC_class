# THIS SCRIPT COMPUTES THE MODEL FIT OF HMSC MODELS FOR THE BIRD EXAMPLE (SECTION 11.1) OF THE BOOK
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

# We next evaluate model fit in terms of explanatory power

preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds)

# We also evaluate model fit in terms of predictive power based on two-fold cross-validation.
# Doing the latter takes a lot of time, as the model needs to be re-fitted twice. For this reason,
# we store the results from cross-validation.
# Set run.cross.validation = TRUE to perform the cross-validation and to save the results.
# Set run.cross.validation = FALSE to read in results from cross-validation that you have run previously.

nParallel = 2
run.cross.validation = TRUE # start with TRUE when you introduce the script
filename=file.path(model.directory, paste0("CV_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(m$thin)))
if(run.cross.validation){
  partition = createPartition(m, nfolds = 2, column = "route")
  preds = computePredictedValues(m,partition=partition, nParallel = nParallel)
  MFCV = evaluateModelFit(hM=m, predY=preds)
  save(partition,MFCV,file=filename)
} else {
  load(filename)
}

tmp = c(mean(MF$AUC), mean(MFCV$AUC), mean(MF$TjurR2), mean(MFCV$TjurR2))
names(tmp)=c("AUC","AUC (CV)","TjurR2","TjurR2 (CV)")
tmp

# We observe that (with thin = 1) the average (over the species) AUC for explanatory power is 0.88,
# and  the average AUC for predictive power is 0.79.
# The average TjurR2 for explanatory power is 0.42,
# and  the average TjurR2  for predictive power is 0.31.
# These results are expected in the sense that explanatory power is typically higher than predictive power,
# and that AUC and TjurR2 are based on different units, AUC being typically higher.
# See the book for more explanation of both of these topics!

# While the code above yields the average results over the species, we may also look at the species-specific results.
# For example, let us plot the explanatory and predictive AUC measures with respect to each other.

plot(MF$AUC,MFCV$AUC)
abline(0,1)

# The call to abline(0,1) adds the identity line (y=x) to the plot. For points (=species) below the line,
# explanatory power is greater than predictive power.

# If we would have fitted multiple models, we could compare them either based on the cross-validation
# (as was done by above) or with WAIC. WAIC can be computed as follows:

WAIC = computeWAIC(m)
WAIC

# As we fitted only one model, WAIC is not informative, as the absolute value is not relevant, just a comparison between models
