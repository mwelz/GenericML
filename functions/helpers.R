make.mlr3.string <- function(learner.str, regr = TRUE){
  # helper function. Requires input of type 'mlr3::lrn("cv_glmnet", s = "lambda.min")' (note the absence of classif and regr)
  
  if(substr(learner.str, start = 1, stop = 6) != "mlr3::"){
    
    learner <- learner.str
    
  } else{
    
    learner <- paste0(substr(learner.str, start = 1, stop = 11), ifelse(regr, "regr.", "classif."), 
                      substr(learner.str, start = 12, stop = 1e8))
    
    learner <- eval(parse(text = learner))
    
  } # END IF
  
  return(learner)
  
} # END FUN



initializer.for.splits <- function(Z, Z_CLAN, learners,
                                   num.splits, quantile.cutoffs,
                                   differences.control_GATES,
                                   differences.control_CLAN){
  
  # helper function that initializes object in the generic ML splitting procedure
  
  if(is.null(Z_CLAN)){
    d <- ncol(Z)
    Z_CLAN.nam <- colnames(Z)
  } else{
    d <- ncol(Z_CLAN)
    Z_CLAN.nam <- colnames(Z_CLAN)
  }
  
  K <- length(quantile.cutoffs) + 1
  CLAN_group.base  <- ifelse(differences.control_CLAN$group.to.subtract.from == "least", 1, K)
  GATES_group.base <- ifelse(differences.control_GATES$group.to.subtract.from == "least", 1, K)
  CLAN_groups.to.be.subtracted  <- differences.control_CLAN$groups.to.be.subtracted
  GATES_groups.to.be.subtracted <- differences.control_GATES$groups.to.be.subtracted
  
  clan <- array(NA_real_, dim = c(K + length(CLAN_groups.to.be.subtracted), 7, num.splits), 
                dimnames = list(c(paste0("delta.", 1:K),
                                  paste0(
                                    "delta.", CLAN_group.base, "-",
                                    "delta.", CLAN_groups.to.be.subtracted)), 
                                c("Estimate", "CB lower", "CB upper", "Std. Error", 
                                  "z value", "Pr(<z)", "Pr(>z)"),
                                NULL))
  
  gates <- array(NA_real_, dim = c(K + length(GATES_groups.to.be.subtracted), 7, num.splits),
                 dimnames = list(c(paste0("gamma.", 1:K),
                                   paste0(
                                     "gamma.", GATES_group.base, "-",
                                     "gamma.", GATES_groups.to.be.subtracted)), 
                                 c("Estimate", "CB lower", "CB upper", "Std. Error", 
                                   "z value", "Pr(<z)", "Pr(>z)"),
                                 NULL))
  
  blp <- array(NA_real_, dim = c(2, 7, num.splits),
               dimnames = list(c("beta.1", "beta.2"), 
                               c("Estimate", "CB lower", "CB upper", "Std. Error", 
                                 "z value", "Pr(<z)", "Pr(>z)"),
                               NULL))
  
  best <- array(NA_real_, dim = c(1, 2, num.splits), 
                dimnames = list(NULL, c("lambda", "lambda.bar"), NULL))
  
  clan.lists <- lapply(1:d, function(...) clan )
  names(clan.lists) <- Z_CLAN.nam
  
  out.ls <- lapply(1:length(learners), function(...){
    
    list(BLP = blp, GATES = gates, CLAN = clan.lists, best = best)
    
  })
  
  names(out.ls) <- learners
  return(out.ls)
  
} # END FUN



initialize.gen.ml <- function(generic.ml.across.learners.obj){
  
  # helper function for initialization of the final returned object
  
  gates.nam  <- rownames(generic.ml.across.learners.obj[[1]]$GATES[,,1])
  z.clan.nam <- names(generic.ml.across.learners.obj[[1]]$CLAN)
  clan.nam   <- rownames(generic.ml.across.learners.obj[[1]]$CLAN[[1]][,,1])
  learners   <- names(generic.ml.across.learners.obj)
  
  blp.mat <- matrix(NA_real_, nrow = 2, ncol = 5, 
                    dimnames = list(c("beta.1", "beta.2"), 
                                    c("Estimate", "CB lower", "CB upper",
                                      "Pr(<z) adjusted", "Pr(>z) adjusted")))
  
  gates.mat <- matrix(NA_real_, nrow = length(gates.nam), ncol = 5, 
                      dimnames = list(gates.nam, 
                                      c("Estimate", "CB lower", "CB upper",
                                        "Pr(<z) adjusted", "Pr(>z) adjusted")))
  
  clan.mat <- matrix(NA_real_, nrow = length(clan.nam), ncol = 5, 
                     dimnames = list(clan.nam, 
                                     c("Estimate", "CB lower", "CB upper", 
                                       "Pr(<z) adjusted", "Pr(>z) adjusted")))
  
  
  gates.ls <- lapply(learners, function(...) gates.mat)
  names(gates.ls) <- learners
  
  blp.ls <- lapply(learners, function(...) blp.mat)
  names(blp.ls) <- learners
  
  z.clan.ls <- lapply(z.clan.nam, function(...) clan.mat)
  names(z.clan.ls) <- z.clan.nam
  
  clan.ls <- lapply(learners, function(...) z.clan.ls)
  names(clan.ls) <- learners
  
  return(list(BLP = blp.ls, GATES = gates.ls, CLAN = clan.ls))
  
} # END FUN



generic.ml.across.learners <- function(Z, D, Y, 
                                       propensity.scores, 
                                       learners, 
                                       num.splits = 50,
                                       Z_CLAN = NULL, 
                                       X1.variables_BLP           = list(functions_of_Z = c("B"),
                                                                         custom_covariates = NULL,
                                                                         fixed_effects = NULL),
                                       X1.variables_GATES         = list(functions_of_Z = c("B"),
                                                                         custom_covariates = NULL,
                                                                         fixed_effects = NULL),
                                       HT.transformation = FALSE,
                                       vcov.control_BLP           = list(estimator = "vcovHC",
                                                                         arguments = list(type = "const")),
                                       vcov.control_GATES         = list(estimator = "vcovHC",
                                                                         arguments = list(type = "const")),
                                       equal.group.variances_CLAN = FALSE,
                                       proportion.in.main.set     = 0.5, 
                                       quantile.cutoffs           = c(0.25, 0.5, 0.75),
                                       differences.control_GATES  = list(group.to.subtract.from = "most",
                                                                         groups.to.be.subtracted = 1),
                                       differences.control_CLAN   = list(group.to.subtract.from = "most",
                                                                         groups.to.be.subtracted = 1),
                                       significance.level = 0.05, 
                                       minimum.variation = 1e-05,
                                       store.learners = FALSE,
                                       store.splits = FALSE){
  
  # initialize
  generic.targets <- initializer.for.splits(Z = Z, Z_CLAN = Z_CLAN, 
                                            learners = learners, num.splits = num.splits, 
                                            quantile.cutoffs = quantile.cutoffs, 
                                            differences.control_GATES = differences.control_GATES,
                                            differences.control_CLAN = differences.control_CLAN)
  
  num.vars.in.Z_CLAN <- ifelse(is.null(Z_CLAN), ncol(Z), ncol(Z_CLAN))
  genericML.by.split <- list()
  N     <- length(Y)
  N.set <- 1:N
   
  # make the custom covariates a matrix to prevent bug later
  if(!is.null(X1.variables_BLP$custom_covariates)){
    X1.variables_BLP$custom_covariates <- as.matrix(X1.variables_BLP$custom_covariates)
  } # IF
  if(!is.null(X1.variables_GATES$custom_covariates)){
    X1.variables_GATES$custom_covariates <- as.matrix(X1.variables_GATES$custom_covariates)
  } # IF
  
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
      
      generic.ml.obj <- 
        get.generic.ml.for.given.learner(Z = Z, D = D, Y = Y, 
                                         propensity.scores = propensity.scores,
                                         learner = learners[i],
                                         M.set = M.set, A.set = A.set,
                                         Z_CLAN                       = Z_CLAN, 
                                         X1.variables_BLP             = X1.variables_BLP,
                                         X1.variables_GATES           = X1.variables_GATES,
                                         HT.transformation            = HT.transformation,
                                         vcov.control_BLP             = vcov.control_BLP,
                                         vcov.control_GATES           = vcov.control_GATES,
                                         equal.group.variances_CLAN   = equal.group.variances_CLAN,
                                         quantile.cutoffs             = quantile.cutoffs,
                                         differences.control_GATES    = differences.control_GATES,
                                         differences.control_CLAN     = differences.control_CLAN,
                                         significance.level           = significance.level,
                                         minimum.variation            = minimum.variation)
        
      generic.targets[[i]]$BLP[,,s]   <- generic.ml.obj$BLP$generic.targets
      generic.targets[[i]]$GATES[,,s] <- generic.ml.obj$GATES$generic.targets
      generic.targets[[i]]$best[,,s]  <- c(generic.ml.obj$best$lambda, generic.ml.obj$best$lambda.bar)
      
      if(store.learners){
        
        genericML.by.split[[learners[i]]][[s]] <- generic.ml.obj
        
      }
      
      for(j in 1:num.vars.in.Z_CLAN){
        generic.targets[[i]]$CLAN[[j]][,,s] <- generic.ml.obj$CLAN$generic.targets[[j]]
      }
      
    } # FOR learners
  } # FOR num.splits
  
  if(!store.learners) genericML.by.split <- NULL
  if(!store.splits)   splits.mat <- NULL
  
  return(list(generic.targets = generic.targets, 
              genericML.by.split = genericML.by.split,
              splits = splits.mat))
  
} # END FUN



get.best.learners <- function(generic.ml.across.learners.obj){
  
  learners <- names(generic.ml.across.learners.obj)  
  
  # for each learner, take medians over the number of splits
  best.analysis <- sapply(learners, 
                          function(learner) apply(generic.ml.across.learners.obj[[learner]]$best, c(1,2), median))
  rownames(best.analysis) <- c("lambda", "lambda.bar")
  
  return(list(best.learner.for.CATE  = learners[which.max(best.analysis["lambda", ])],
              best.learner.for.GATES = learners[which.max(best.analysis["lambda.bar", ])],
              lambda.overview = t(best.analysis)))
  
} # END FUN



# TODO: make documentation; generic.ml.across.learners.obj is NOT the correct description for the input anymore (check genericML function!)
VEIN <- function(generic.ml.across.learners.obj, best.learners.obj){
  
  gen.ml.ls <- initialize.gen.ml(generic.ml.across.learners.obj)
  learners  <- names(generic.ml.across.learners.obj)
  
  for(learner in learners){
    
    # BLP and GATES
    for(type in c("BLP", "GATES")){
      
      gen.ml.ls[[type]][[learner]][,"Estimate"] <- 
        apply(generic.ml.across.learners.obj[[learner]][[type]][,"Estimate",], 1, function(z) Med(z)$Med)
      gen.ml.ls[[type]][[learner]][,"CB lower"] <- 
        apply(generic.ml.across.learners.obj[[learner]][[type]][,"CB lower",], 1, function(z) Med(z)$upper.median)
      gen.ml.ls[[type]][[learner]][,"CB upper"] <- 
        apply(generic.ml.across.learners.obj[[learner]][[type]][,"CB upper",], 1, function(z) Med(z)$lower.median)
      p.left.raw <- 
        apply(generic.ml.across.learners.obj[[learner]][[type]][,"Pr(<z)",], 1, function(z) Med(z)$lower.median)
      p.right.raw <- 
        apply(generic.ml.across.learners.obj[[learner]][[type]][,"Pr(>z)",], 1, function(z) Med(z)$lower.median)
      p.left.adj  <- 2 * p.left.raw
      p.right.adj <- 2 * p.right.raw
      p.left.adj[p.left.adj > 1]   <- 1 # p-values cannot exceed 1
      p.right.adj[p.right.adj > 1] <- 1
      gen.ml.ls[[type]][[learner]][,"Pr(<z) adjusted"] <- p.left.adj
      gen.ml.ls[[type]][[learner]][,"Pr(>z) adjusted"] <- p.right.adj
      
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
      p.left.raw <- 
        apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"Pr(<z)",], 1, function(z) Med(z)$lower.median)
      p.right.raw <- 
        apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"Pr(>z)",], 1, function(z) Med(z)$lower.median)
      p.left.adj  <- 2 * p.left.raw
      p.right.adj <- 2 * p.right.raw
      p.left.adj[p.left.adj > 1]   <- 1 # p-values cannot exceed 1
      p.right.adj[p.right.adj > 1] <- 1
      gen.ml.ls$CLAN[[learner]][[z.clan]][,"Pr(<z) adjusted"] <- p.left.adj
      gen.ml.ls$CLAN[[learner]][[z.clan]][,"Pr(>z) adjusted"] <- p.right.adj
      
    } # FOR Z_CLAN
    
  } # FOR learners
  
  
  return(list(best.learners = list(BLP = gen.ml.ls$BLP[[best.learners.obj$best.learner.for.CATE]],
                                   GATES = gen.ml.ls$GATES[[best.learners.obj$best.learner.for.GATES]],
                                   CLAN = gen.ml.ls$CLAN[[best.learners.obj$best.learner.for.GATES]]),
              all.learners = gen.ml.ls))
  
} # END FUN


# helper function that calculates an error covariance matrix estimator of a linear model
#
## @param x a linear model object
## @param vcov.control a list with two elements called 'estimator' and 'arguments'. The argument 'estimator' is a string specifying the covariance matrix estimator to be used; specifies a covariance estimator function in the sandwich package (https://cran.r-project.org/web/packages/sandwich/sandwich.pdf). Recommended estimators are "vcovBS", "vcovCL", "vcovHAC", and "vcovHC". Default is 'vcovHC'. The element 'arguments' is a list of arguments that shall be passed to the function specified in the element 'estimator'. Default leads to the (homoskedastic) ordinary least squares covariance matrix estimator. See the reference manual of the sandwich package for details (https://cran.r-project.org/web/packages/sandwich/vignettes/sandwich.pdf).
get.vcov <- function(x,
                     vcov.control = list(estimator = "vcovHC",
                                         arguments = list(type = "const"))){
  
  # append the model so that we can pass this list to do.call
  arguments   <- vcov.control$arguments
  arguments$x <- x
  
  # return the estimate
  do.call(what = eval(parse(text = paste0("sandwich::", vcov.control$estimator))),
          args = arguments)

} # FUN


# helper function to prepare the custom part of the regressor matrix in BLP and GATES
get.df.from.X1.variables <- function(functions.of.Z_mat,
                                     X1.variables){
  
  custom        <- X1.variables$custom_covariates
  fixed.eff     <- X1.variables$fixed_effects
  out           <- data.frame(functions.of.Z_mat[, X1.variables$functions_of_Z])
  colnames(out) <- X1.variables$functions_of_Z
  
  if(!is.null(fixed.eff)) out$fixed.effects <- factor(fixed.eff)
  if(!is.null(custom))    out <- data.frame(out, custom)
  
  out
  
} # FUN