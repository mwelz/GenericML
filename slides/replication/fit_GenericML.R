# ------------------------------------------------------------------------------
# This script applies GenericML() on the preprocessed data of
# a study by Crépon et al. (2015, AER:AE, DOI: 10.1257/app.20130535).
# This randomized study investigated the effects of microcredit availability on
# household economic behavior in rural Morocco.
#
# The script outputs a file "GenericML_object.Rdata".
#
# The dataset was kindly made available by Esther Duflo.
#
# Note: your machine requires at least 6 cores to run this script.
# You may change "num_cores" to less cores (or even opt for no parallelization),
# but the results in the slides are only replicable for 6 cores on Unix machines
# (by virtue of parallel seeding)
#
# This script was written by M. WELZ (welz@ese.eur.nl)
# ------------------------------------------------------------------------------


## install (if applicable) and load relevant packages
required_packages <- c("ranger", "glmnet", "e1071", "xgboost", "GenericML")
missing <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing) > 0L) install.packages(missing)

## load GenericML
library(GenericML, quietly = TRUE)

## load data
load("slides/data/morocco_preprocessed.Rdata")

## prepare learners
learners <- c("random_forest",
              'mlr3::lrn("cv_glmnet", s = "lambda.min", alpha = 0.5)',
              "mlr3::lrn('svm')",
              "mlr3::lrn('xgboost')")


## setup_X1() customizes inclusion of controls and fixed effects
# include BCA and CATE controls
# add fixed effects along variable "vil_pair"
X1 <- setup_X1(funs_Z = c("B", "S"),
               fixed_effects = vil_pair)

## setup_vcov() customizes covariance estimation
# calls functions from the "sandwich" package
# cluster standard errors along "demi_paire"
vcov <- setup_vcov(estimator = "vcovCL",
                   arguments = list(cluster = demi_paire))

## run GenericML()
x <- GenericML(
  Z = Z, D = D, Y = Y,                      # observed data
  learners_GenericML = learners,            # learners
  learner_propensity_score = "constant",    # = 0.5 (RCT)
  num_splits = 100L,                        # number splits
  quantile_cutoffs = c(0.2, 0.4, 0.6, 0.8), # grouping
  significance_level = 0.05,                # significance level
  X1_BLP = X1, X1_GATES = X1,               # regression setup
  vcov_BLP = vcov, vcov_GATES = vcov,       # covariance setup
  parallel = TRUE, num_cores = 6L,          # parallelization
  seed = 20220621)                          # RNG seed

## save GenericML() object
save(x, file = "slides/replication/GenericML_object.Rdata")
