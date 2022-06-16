# ------------------------------------------------------------------------------
# This script loads and preprocesses the raw dataset of
# a study by Crépon et al. (2015, AER:AE, DOI: 10.1257/app.20130535).
# This randomized study investigated the effects of microcredit availability on
# household economic behavior in rural Morocco.
#
# The script outputs a preprocessed file, "morocco_preprocessed.Rdata".
#
# The dataset was kindly made available by Esther Duflo.
#
# This script was written by M. WELZ (welz@ese.eur.nl) and is based on a script
# of V. CHERNOZHUKOV, M. DEMIRER, E. DUFLO, and I. FERNANDEZ-VAL.
# ------------------------------------------------------------------------------


### 0. install (if applicable) the 'readstata13' package ----
if(isFALSE("readstata13" %in% installed.packages()[,"Package"])){
  install.packages("readstata13")
} # IF


### 1. preprocessing ----
## load and encode data
data        <- readstata13::read.dta13("slides/data/morocco_raw.dta")
meta        <- attributes(data)
data$paire  <- factor(data$paire)
a           <- as.data.frame(model.matrix(~data$paire-1L))
colnames(a) <- (substring( names(a), 6L, 12L))
data        <- cbind(data, a)
colnames(data)[which(colnames(data)=="paire")] <- "vil_pair"

## prepare control variables
controls_vars  <- c("members_resid_bl", "nadults_resid_bl",
                    "head_age_bl", "act_livestock_bl", "act_business_bl",
                    "borrowed_total_bl", "members_resid_d_bl",
                    "nadults_resid_d_bl", "head_age_d_bl",
                    "act_livestock_d_bl", "act_business_d_bl",
                    "borrowed_total_d_bl", "ccm_resp_activ",
                    "other_resp_activ", "ccm_resp_activ_d",
                    "other_resp_activ_d")
controls_paire <- names(data)[(substring( names(data), 1L, 5L) == "paire")]
controls_description <- cbind(meta$names[meta$names %in% controls_vars],
                              description = meta$var.labels[meta$names %in% controls_vars])
controls_vars <- controls_description[, 1L] # reorder
controls      <- c(controls_vars, controls_paire)

## vector of outcome variables
Y_nams <- c("loansamt_total", "output_total", "profit_total", "consumption")
Y_description <- cbind(Y_nams, description = meta$var.labels[meta$names %in% Y_nams])

## 'paire' is a one-hot encoded ID variable with 81 levels
vil_paire_description  <- meta$var.labels[meta$names %in% "paire"]
demi_paire_description <- meta$var.labels[meta$names %in% "demi_paire"]
group_description <- cbind(variable = c("vil_pair", "demi_paire"),
                           description = c(vil_paire_description, demi_paire_description))

## column names
colnames(Y_description) <- colnames(controls_description) <-
  c("variable", "description")

## vector of treatment indicators
D_nam <- "treatment"

## ID variables; useful for fixed effects and clustering
group_vars <- c("vil_pair", "demi_paire")

## drop NAs for reduced data
data_red <- na.omit(data[, c(Y_nams, D_nam, controls, group_vars)])


### 2. print variable description ----
controls_description
Y_description
group_description


### 3. gather preprocessed data ----
## dependent variable: total borrowing here
Y <- data_red[, "loansamt_total"]

## treatment assignment
D <- data_red$treatment

## covariates
# drop last 'paire' one-hot category to avoid singularity
Z <- as.matrix(data_red[, controls[-length(controls)]])

## group variables
vil_pair   <- as.integer(data_red$vil_pair)
demi_paire <- as.integer(data_red$demi_paire)

## description
description <- list(Y = Y_description[1L,,drop = FALSE],
                    Z = controls_description,
                    group = group_description,
                    D = "1 if household has access to microcredits, 0 otherwise")

## save
save(Z, Y, D, vil_pair, demi_paire, description,
     file = "slides/data/morocco_preprocessed.Rdata")
