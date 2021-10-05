# helper functions
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
                                       parallel = .Platform$OS.type == "unix",
                                       num.cores = parallel::detectCores(), 
                                       seed = NULL,
                                       store.learners = FALSE,
                                       store.splits = FALSE){
  
  # call correct main function
  do.call(what = get(ifelse(parallel, 
                            "generic.ml.across.learners_parallel", 
                            "generic.ml.across.learners_serial")),
          args = list(Z = Z, D = D, Y = Y, 
                      propensity.scores          = propensity.scores, 
                      learners                   = learners, 
                      learners.names             = learners.names,
                      num.splits                 = num.splits,
                      Z_CLAN                     = Z_CLAN, 
                      X1.variables_BLP           = X1.variables_BLP,
                      X1.variables_GATES         = X1.variables_GATES,
                      HT.transformation          = HT.transformation,
                      vcov.control_BLP           = vcov.control_BLP,
                      vcov.control_GATES         = vcov.control_GATES,
                      equal.group.variances_CLAN = equal.group.variances_CLAN,
                      proportion.in.main.set     = proportion.in.main.set, 
                      quantile.cutoffs           = quantile.cutoffs,
                      differences.control_GATES  = differences.control_GATES,
                      differences.control_CLAN   = differences.control_CLAN,
                      significance.level         = significance.level, 
                      minimum.variation          = minimum.variation,
                      num.cores                  = num.cores, 
                      seed                       = seed,
                      store.learners             = store.learners,
                      store.splits               = store.splits))
  
} # FUN



generic.ml.across.learners_serial <- function(Z, D, Y, 
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
                                              num.cores = parallel::detectCores(), # dead argument here
                                              seed = NULL,
                                              store.learners = FALSE,
                                              store.splits = FALSE){
  

  # if no input provided, set Z_CLAN it equal to Z
  if(is.null(Z_CLAN)) Z_CLAN <- Z 
  
  # set seed
  if(is.null(seed)){
    
    RNGkind("L'Ecuyer-CMRG")
    
  } else{
    
    set.seed(seed, "L'Ecuyer")
    
  } # IF
  
  num.vars.in.Z_CLAN <- ncol(Z_CLAN)
  genericML.by.split <- list()
  N     <- length(Y)
  N.set <- 1:N
  prop <- proportion.in.main.set * N 
  
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
  if(store.splits) splits.mat <- matrix(NA_character_, N, num.splits, 
                                        dimnames = c(NULL, paste0("split_", 1:num.splits)))
  
  # initialize
  generic.targets <- initializer.for.splits(Z = Z, Z_CLAN = Z_CLAN, 
                                            learners = learners.names, num.splits = num.splits, 
                                            quantile.cutoffs = quantile.cutoffs, 
                                            differences.control_GATES = differences.control_GATES,
                                            differences.control_CLAN = differences.control_CLAN)
  
  # loop over the sample splits
  for(s in 1:num.splits){
    
    # perform sample splitting into main set and auxiliary set
    M.set <- sort(sample(x = N.set, size = floor(prop), replace = FALSE),
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




generic.ml.across.learners_parallel <- function(Z, D, Y,
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
                                              num.cores = parallel::detectCores(),
                                              seed = NULL,
                                              store.learners = FALSE,
                                              store.splits = FALSE){

  # if no input provided, set Z_CLAN it equal to Z
  if(is.null(Z_CLAN)) Z_CLAN <- Z

  # set seed
  if(is.null(seed)){

    RNGkind("L'Ecuyer-CMRG")

  } else{

    set.seed(seed, "L'Ecuyer")

  } # IF

  num.vars.in.Z_CLAN <- ncol(Z_CLAN)
  K <- length(quantile.cutoffs) + 1
  num.generic.targets.gates <- K + length(differences.control_GATES$groups.to.be.subtracted)
  num.generic.targets.clan <- K + length(differences.control_CLAN$groups.to.be.subtracted)
  N     <- length(Y)
  N.set <- 1:N
  num.learners <- length(learners)
  prop <- proportion.in.main.set * N

  # set variable names fo CLAN
  if(is.null(colnames(Z_CLAN))) colnames(Z_CLAN) <- paste0("V", 1:num.vars.in.Z_CLAN)

  # make the custom covariates a matrix to prevent bug later
  if(!is.null(X1.variables_BLP$custom_covariates)){
    X1.variables_BLP$custom_covariates <- as.matrix(X1.variables_BLP$custom_covariates)
  } # IF
  if(!is.null(X1.variables_GATES$custom_covariates)){
    X1.variables_GATES$custom_covariates <- as.matrix(X1.variables_GATES$custom_covariates)
  } # IF


  # loop over the sample splits
  out <- parallel::mclapply(mc.cores = num.cores, X = 1:num.splits, FUN = function(s){

    # initialize
    blp.3d     <- get_blp.3d(num.learners, learners.names)
    gates.3d   <- get_gates.3d(num.learners, learners.names, num.generic.targets.gates)
    best.3d    <- get_best.3d(num.learners, learners.names)
    clan.3d.ls <- get_clan.3d.ls(num.learners, learners.names, num.generic.targets.clan, num.vars.in.Z_CLAN, colnames(Z_CLAN))


    if(store.learners){
      genericML.by.split.ls <- rep(list(NA_real_), num.learners)
    } else{
      genericML.by.split.ls <- NULL
    } # IF


    # perform sample splitting into main set and auxiliary set
    M.set <- sort(sample(x = N.set, size = floor(prop), replace = FALSE),
                  decreasing = FALSE)
    A.set <- setdiff(N.set, M.set)

    if(store.splits){

      splits.vec        <- rep("M", N)
      splits.vec[A.set] <- "A"

    } else{
      splits.vec <- NULL
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

      blp.3d[,,i]   <- generic.ml.obj$BLP$generic.targets
      gates.3d[,,i] <- generic.ml.obj$GATES$generic.targets
      best.3d[,,i]  <- c(generic.ml.obj$best$lambda, generic.ml.obj$best$lambda.bar)

      for(j in 1:num.vars.in.Z_CLAN){
        clan.3d.ls[[j]][,,i] <- generic.ml.obj$CLAN$generic.targets[[j]]
      }

      if(store.learners) genericML.by.split.ls[[i]] <- generic.ml.obj

    } # FOR learners

    # return
    return(
    list(BLP = blp.3d, GATES = gates.3d, best = best.3d,
         CLAN = clan.3d.ls, GML = genericML.by.split.ls,
         splits = splits.vec))

  }) # MCLAPPLY


  # now bring it in desired form:
  generic.targets <- initializer.for.splits(Z = Z, Z_CLAN = Z_CLAN,
                                            learners = learners.names, num.splits = num.splits,
                                            quantile.cutoffs = quantile.cutoffs,
                                            differences.control_GATES = differences.control_GATES,
                                            differences.control_CLAN = differences.control_CLAN)

  for(s in 1:num.splits){
    for(i in 1:length(learners)){

      generic.targets[[i]]$BLP[,,s]   <- out[[s]]$BLP[,,i]
      generic.targets[[i]]$GATES[,,s] <- out[[s]]$GATES[,,i]
      generic.targets[[i]]$best[,,s]  <- out[[s]]$best[,,i]

      for(j in 1:num.vars.in.Z_CLAN){
        generic.targets[[i]]$CLAN[[j]][,,s] <- out[[s]]$CLAN[[j]][,,i]
      } # FOR

    } # FOR
  } # FOR


  if(store.splits){
    splits <- sapply(1:num.splits, function(s) out[[s]]$splits)
    colnames(splits) <- paste0("split_", 1:num.splits)
  } else{
    splits <- NULL
  } # IF

  if(store.learners){

    stored.learners <- list()

    for(s in 1:num.splits){
      for(i in 1:length(learners)){
        stored.learners[[learners.names[i]]][[s]] <- out[[s]]$GML[[i]]
      }
    }

  } else{
    stored.learners <- NULL
  } # IF

  # return
  return(list(generic.targets = generic.targets,
              genericML.by.split = stored.learners,
              splits = splits))

} # FUN
