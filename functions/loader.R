# temporary function to load everything in the workspace. Will be superfluous once repo is a package

library(ggplot2)
library(mlr3)
library(mlr3learners) # potentially a bug in mlr3; if this is not loaded, mlr::lrn won't recognize the learner class

source(paste0(getwd(), "/functions/GenericML.R"))
source(paste0(getwd(), "/functions/helpers.R"))
source(paste0(getwd(), "/functions/clan.R"))
source(paste0(getwd(), "/functions/blp.R"))
source(paste0(getwd(), "/functions/gates.R"))
source(paste0(getwd(), "/functions/misc.R"))
source(paste0(getwd(), "/functions/mlr3.R"))
source(paste0(getwd(), "/functions/plot.R"))
source(paste0(getwd(), "/functions/S3.R"))

