rm(list = ls()) ; cat("\014")
library(glmnet) # apparently needs to be loaded (?) TODO
library(ggplot2)

# load the functions
source(paste0(getwd(), "/functions/generic-ml-estimation-funs.R"))
source(paste0(getwd(), "/functions/generic-ml-auxiliary-funs.R"))

set.seed(1)
num.obs  <- 500
num.vars <- 5

# ate parameter
theta <- -2

# treatment assignment (assume RCT)
D <- rbinom(num.obs, 1, 0.5) 

# covariates
Z <- mvtnorm::rmvnorm(num.obs, mean = rep(0, num.vars), sigma = diag(num.vars))
colnames(Z) <- paste0("z", 1:num.vars)

# coefficients (including an intercept)
beta <- c(1, 2, -3, 0, 0, 2)

# compute Pr(Y = 1 | X) for each individual (with noise)
#eps     <-  rnorm(n, mean = 0, sd = 0.5)
y0.star <- as.numeric(cbind(1, Z) %*% beta) #+ eps
y1.star <- theta + y0.star

# experimental: use standardized logit link
mu <- 0 # mean(y0.star)
s <- 1 #sqrt(var(y0.star))
pi0 <- plogis(y0.star, location = mu, scale = s)
pi1 <- plogis(y1.star, location = mu, scale = s)

# hte
hte <- pi1 - pi0

# true ATE
ate <- mean(pi1) - mean(pi0)

# create binary outcomes
Y0 <- rbinom(num.obs, 1, pi0)
Y1 <- rbinom(num.obs, 1, pi1)
Y  <- ifelse(D == 1, Y1, Y0) # observed outcome


#######################

# arguments: 
quantile.cutoffs         <- c(0.2, 0.4, 0.6, 0.8) # for the GATES grouping of S (argument)
proportion.in.main.set   <- 0.5 # argument
Z.clan                   <- NULL # argument. The matrix of variables that shall be considered in CLAN
learners.genericML       <- c('glm', 'mlr3::lrn("ranger", num.trees = 100)')
learner.propensity.score <- 'mlr3::lrn("glmnet", lambda = 0, alpha = 1)' # non-penalized logistic regression
num.splits               <- 3
significance.level       <- 0.05
store.splits             <- FALSE
store.learners           <- FALSE

# TODO: The TODOs in genericML()
genML <- genericML(Z = Z, D = D, Y = Y, 
                   learner.propensity.score = learner.propensity.score, 
                   learners.genericML = learners.genericML,
                   num.splits = num.splits,
                   Z.clan = Z.clan,
                   quantile.cutoffs = quantile.cutoffs,
                   proportion.in.main.set = proportion.in.main.set, 
                   significance.level = significance.level,
                   store.splits = store.splits,
                   store.learners = store.learners)

# save object
# save(genML, hte, ate, Y, D, Z, file = paste0(getwd(), "/tests-ml/gen-ml-test.Rdata"))

# analyze
genML$VEIN$best.learners$GATES # difference is insignificant, so no hetero
genML$VEIN$best.learners$BLP  # beta2 is insignificant, so no hetero
genML$VEIN$best.learners$CLAN$z1 # there seems to be hetero along z1

# GATES
genericML.plot(genML, type = "GATES") # no hetero
genericML.plot(genML, type = "BLP")   # no hetero

genericML.plot(genML, type = "CLAN", CLAN.variable = "z1")   # no hetero
genericML.plot(genML, type = "CLAN", CLAN.variable = "z2")   # hetero
genericML.plot(genML, type = "CLAN", CLAN.variable = "z3")   # slight hetero
genericML.plot(genML, type = "CLAN", CLAN.variable = "z4")   # no hetero
genericML.plot(genML, type = "CLAN", CLAN.variable = "z5")   # no hetero
