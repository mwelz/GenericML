rm(list = ls()) ; gc(); cat("\014")

# Install the required packages if they are not already installed
required.packages <- c("ggplot2", "mlr3", "mlr3learners", 
                       "mvtnorm", "glmnet", "ranger", "e1071",
                       "kknn", "MASS", "nnet", "xgboost")
new.packages      <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
rm(required.packages, new.packages)


# load the functions
source(paste0(getwd(), "/functions/generic-ml-estimation-funs.R"))
source(paste0(getwd(), "/functions/generic-ml-auxiliary-funs.R"))

### 1. Data Generation (linear, no treatment effect heterogeneity) ----
set.seed(1)
num.obs  <- 10000
num.vars <- 5

# ATE parameter
ATE <- 2

# random treatment assignment 
D <- rbinom(num.obs, 1, 0.5) 

# covariates
Z <- mvtnorm::rmvnorm(num.obs, mean = rep(0, num.vars), sigma = diag(num.vars))
colnames(Z) <- paste0("z", 1:num.vars)

# counterfactual outcomes
Y0 <- as.numeric(Z %*% c(2, -3, 0, 1, 2)) + rnorm(num.obs)
Y1 <- ATE + Y0

# observed outcome
Y  <- ifelse(D == 1, Y1, Y0) 


### 2. Prepare the arguments for genericML() ---- 

# quantile cutoffs for the GATES grouping of the estiimated CATEs 
quantile.cutoffs         <- c(0.2, 0.4, 0.6, 0.8) # 20%, 40%, 60%, 80% quantiles

# specify the learner of the propensity score (observe the mlr3 syntax)
# any mlr3 learner can be specified (list: https://mlr3learners.mlr-org.com/)
learner.propensity.score <- "mlr3::lrn('glmnet', lambda = 0, alpha = 1)" # non-penalized logistic regression

# specify the considered learners of the BCA and the CATE.  Any mlr3 learner can be specified (list: https://mlr3learners.mlr-org.com/)
# three learners here: elastic net, random forest, and SVM
learners.genericML       <- c("elastic.net", "mlr3::lrn('ranger', num.trees = 100)", "mlr3::lrn('svm')")

# specify the data that shall be used for the CLAN
# here, we use all variables of Z and uniformly distributed noise
Z.clan <- cbind(Z, random = runif(num.obs))

# specify the number of splits 
num.splits               <- 100

# specify the significance level
significance.level       <- 0.05

# specify the proportion of samples that shall be selected in the main set
proportion.in.main.set   <- 0.5

# specify whether or not the splits and auxiliary results of the learners shall be stored
store.splits             <- FALSE
store.learners           <- FALSE


### step 1: compute propensity scores ----
propensity.scores.obj <- propensity.score(Z = Z, D = D, 
                                          learner = make.mlr3.string(learner.propensity.score, 
                                                                     regr = FALSE))
propensity.scores     <- propensity.scores.obj$propensity.scores

### step 2: for each ML method, do the generic ML analysis ----
learners = learners.genericML

# initialize
generic.targets <- initializer.for.splits(Z = Z, Z.clan = Z.clan, 
                                          learners = learners, num.splits = num.splits, 
                                          quantile.cutoffs = quantile.cutoffs)

num.vars.in.Z.clan <- ifelse(is.null(Z.clan), ncol(Z), ncol(Z.clan))
genericML.by.split <- list()
N     <- length(Y)
N.set <- 1:N

if(store.splits) splits.mat <- matrix(NA_character_, N, num.splits)

# loop over the sample splits
for(s in 1:num.splits){
  
  # perform sample splitting into main set and auxiliary set
  M.set <- sort(sample(x = N.set, size = floor(proportion.in.main.set * N), replace = FALSE),
                decreasing = FALSE)
  A.set <- setdiff(N.set, M.set)
  
  if(store.splits){
    
    splits.mat[M.set, s] <- "M"
    splits.mat[A.set, s] <- "A"
    
  } # IF
  
  
  # loop over the learners
  for(i in 1:length(learners)){
    
    learner = learners[i]
    
    
    ### step 1: input checks ---- 
    if(is.null(Z.clan)) Z.clan <- Z # if no input provided, set it equal to Z
    
    ### step 2a: learn proxy predictors by using the auxiliary set ----
    
    # get the proxy baseline estimator for the main sample
    proxy.baseline.obj <- baseline.proxy.estimator(Z = Z, D = D, Y = Y, 
                                                   auxiliary.sample = A.set, 
                                                   learner = make.mlr3.string(learner, regr = TRUE))
    proxy.baseline     <- proxy.baseline.obj$baseline.predictions.main.sample
    
    # get the proxy estimator of the CATE for the main sample
    proxy.cate.obj <- 
      CATE.proxy.estimator(Z = Z, D = D, Y = Y,
                           auxiliary.sample = A.set, 
                           learner = make.mlr3.string(learner, regr = TRUE),
                           proxy.baseline.estimates = proxy.baseline.obj$baseline.predictions.full.sample)
    proxy.cate <- proxy.cate.obj$CATE.predictions.main.sample
    
    
    ### step 2b: estimate BLP parameters by OLS (TODO: HT transformation!) ----
    blp.obj <- get.BLP.params.classic(D = D[M.set], Y = Y[M.set],
                                      propensity.scores = propensity.scores[M.set],
                                      proxy.baseline = proxy.baseline, 
                                      proxy.cate = proxy.cate, 
                                      significance.level = significance.level)
    
    
    ### step 2c: estimate GATES parameters by OLS (TODO: HT transformation!) ----
    # group the proxy estimators for the CATE in the main sample by quantiles. TODO: intervals need to be [) instead of (]
    group.membership.main.sample <- quantile.group(proxy.cate, 
                                                   cutoffs = quantile.cutoffs, 
                                                   quantile.nam = TRUE) 
    
    gates.obj <- get.GATES.params.classic(D = D[M.set], Y = Y[M.set],
                                          propensity.scores = propensity.scores[M.set],
                                          group.membership.main.sample = group.membership.main.sample, 
                                          proxy.baseline = proxy.baseline, proxy.cate = proxy.cate,
                                          significance.level = significance.level)
    
    
    ### step 2d: estimate CLAN parameters in the main sample
    clan.obj <- get.CLAN.parameters(Z.clan.main.sample = Z.clan[M.set,], 
                                    group.membership.main.sample = group.membership.main.sample)
    
    
    ### step 2e: get parameters over which we maximize to find the "best" ML method ----
    best.obj <- best.ml.method.parameters(BLP.obj = blp.obj, 
                                          GATES.obj = gates.obj, 
                                          proxy.cate.main.sample = proxy.cate,
                                          group.membership.main.sample = group.membership.main.sample)
    
    
    
    
    generic.targets[[i]]$BLP[,,s]   <- blp.obj$generic.targets
    generic.targets[[i]]$GATES[,,s] <- gates.obj$generic.targets
    generic.targets[[i]]$best[,,s]  <- c(best.obj$lambda, best.obj$lambda.bar)
    
    if(store.learners){
      
      genericML.by.split[[learners[i]]][[s]] <- generic.ml.obj
      
    }
    
    for(j in 1:num.vars.in.Z.clan){
      generic.targets[[i]]$CLAN[[j]][,,s] <- clan.obj$generic.targets[[j]]
    }
    
  } # FOR learners
} # FOR num.splits





# extract the best learners
best.learners <- get.best.learners(generic.targets)

### step 3: perform VEIN analysis ---- 
#vein <- VEIN(generic.targets, best.learners)
gen.ml.ls <- initialize.gen.ml(generic.targets)
generic.ml.across.learners.obj = generic.targets

for(learner in learners){
  
  # BLP and GATES
  for(type in c("BLP", "GATES")){
    
    gen.ml.ls[[type]][[learner]][,"Estimate"] <- 
      apply(generic.ml.across.learners.obj[[learner]][[type]][,"Estimate",], 1, function(z) Med(z)$Med)
    gen.ml.ls[[type]][[learner]][,"CB lower"] <- 
      apply(generic.ml.across.learners.obj[[learner]][[type]][,"CB lower",], 1, function(z) Med(z)$upper.median)
    gen.ml.ls[[type]][[learner]][,"CB upper"] <- 
      apply(generic.ml.across.learners.obj[[learner]][[type]][,"CB upper",], 1, function(z) Med(z)$lower.median)
    pval <- 
      apply(generic.ml.across.learners.obj[[learner]][[type]][,"Pr(>|z|)",], 1, function(z) Med(z)$lower.median)
    gen.ml.ls[[type]][[learner]][,"p-value raw"] <- pval
    pval.adj <- 2 * pval
    pval.adj[pval.adj > 1] <- 1 # p-values cannot exceed 1
    gen.ml.ls[[type]][[learner]][,"p-value adjusted"] <- pval.adj
    
  } # FOR type
  
  

  
  
  # CLAN
  z.clan.nam <- names(generic.ml.across.learners.obj[[1]]$CLAN)
  
  for(z.clan in z.clan.nam){
    
    gen.ml.ls$CLAN[[learner]][[z.clan]][,"Estimate"] <- 
      apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"Estimate",], 1, function(z) Med(z)$Med)
    gen.ml.ls$CLAN[[learner]][[z.clan]][,"CB lower"] <- 
      apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"CB lower",], 1, function(z) Med(z)$upper.median)
    gen.ml.ls$CLAN[[learner]][[z.clan]][,"CB upper"] <- 
      apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"CB upper",], 1, function(z) Med(z)$upper.median)
    pval <- 
      apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"Pr(>|z|)",], 1, function(z) Med(z)$lower.median)
    gen.ml.ls$CLAN[[learner]][[z.clan]][,"p-value raw"] <- pval
    pval.adj <- 2 * pval
    pval.adj[pval.adj > 1] <- 1 # p-values cannot exceed 1
    gen.ml.ls$CLAN[[learner]][[z.clan]][,"p-value adjusted"] <- pval.adj
    
  } # FOR Z.clan
  
} # FOR learners

