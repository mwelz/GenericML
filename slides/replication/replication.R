# -------------------------------------------------------------------------------
# This script replicates all results in the slides of our talk at useR! 2022
# You can find the slides as "useR2022.pdf" in the parent folder of this script
#
# Note: your machine requires at least 6 cores to run this script.
# You may change "num_cores" to less cores (or even opt for no parallelization),
# but the results in the slides are only replicable for 6 cores on Unix machines
# (by virtue of parallel seeding)
#
# This script was written by M. WELZ (welz@ese.eur.nl)
# -------------------------------------------------------------------------------

# install (if applicable) and load relevant packages
required_packages <- c("ranger", "glmnet", "e1071", "xgboost", "GenericML", "devtools")
missing <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
# if(length(missing) > 0L) install.packages(missing) # uncomment to install

# install version 0.2.3 which is not yet on CRAN
# devtools::install_github("mwelz/GenericML")
library("GenericML")

# load data, available in GitHub repo mwelz/GenericML
url_data <-
  url(paste0(
    "https://github.com/mwelz/GenericML/blob/main/slides",
    "/data/morocco_preprocessed.Rdata?raw=true"
  ))
load(url_data)

# specify learners
learners <-
  c("random_forest",
    "mlr3::lrn('cv_glmnet', s = 'lambda.min', alpha = 0.5)",
    "mlr3::lrn('svm')",
    "mlr3::lrn('xgboost')")


# include BCA and CATE controls
# add fixed effects along variable "vil_pair"
X1 <- setup_X1(funs_Z = c("B", "S"),
               fixed_effects = vil_pair)


# calls functions from the "sandwich" package
# cluster standard errors along "demi_paire"
vcov <- setup_vcov(estimator = "vcovCL",
                   arguments = list(cluster = demi_paire))


# run GenericML()
# load("slides/replication/GenericML_object.Rdata") # uncomment if you want to load the object below
genML <- GenericML(
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

# save GenericML() object to the repo
# save(genML, file = "slides/replication/GenericML_object.Rdata")

# BLP
results_BLP <- get_BLP(genML, plot = TRUE)
results_BLP       # print method
plot(results_BLP) # plot method

# GATES
results_GATES <- get_GATES(genML, plot = TRUE)
results_GATES
plot(results_GATES)

# CLAN
results_CLAN <- get_CLAN(genML, variable = "head_age_bl", plot = TRUE)
results_CLAN
plot(results_CLAN)

# best learners
get_best(genML)
