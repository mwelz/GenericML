# helper function
# the learners need to be mlr3 objects!
generic.ml.across.learners <- function(Z, D, Y, 
                                       propensity.scores, 
                                       learners, # need to be mlr3 objects!
                                       learners.names,
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
  

  # if no input provided, set Z_CLAN it equal to Z
  if(is.null(Z_CLAN)) Z_CLAN <- Z 
  
  num.vars.in.Z_CLAN <- ncol(Z_CLAN)
  genericML.by.split <- list()
  N     <- length(Y)
  N.set <- 1:N
  
  # set variable names fo CLAN
  if(is.null(colnames(Z_CLAN))) colnames(Z_CLAN) <- paste0("V", 1:num.vars.in.Z_CLAN)
  
  # make the custom covariates a matrix to prevent bug later
  if(!is.null(X1.variables_BLP$custom_covariates)){
    X1.variables_BLP$custom_covariates <- as.matrix(X1.variables_BLP$custom_covariates)
  } # IF
  if(!is.null(X1.variables_GATES$custom_covariates)){
    X1.variables_GATES$custom_covariates <- as.matrix(X1.variables_GATES$custom_covariates)
  } # IF
  
  # initialize
  if(store.splits) splits.mat <- matrix(NA_character_, N, num.splits)
  
  # initialize
  generic.targets <- initializer.for.splits(Z = Z, Z_CLAN = Z_CLAN, 
                                            learners = learners.names, num.splits = num.splits, 
                                            quantile.cutoffs = quantile.cutoffs, 
                                            differences.control_GATES = differences.control_GATES,
                                            differences.control_CLAN = differences.control_CLAN)
  
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
        get.generic.ml.for.given.learner_NoChecks(Z = Z, D = D, Y = Y, 
                                                  propensity.scores = propensity.scores,
                                                  learner = learners[[i]],
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
        
        genericML.by.split[[learners.names[i]]][[s]] <- generic.ml.obj 
        
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
