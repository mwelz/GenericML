# helper functions
generic.ml.across.learners <- function(Z, D, Y,
                                       propensity_scores,
                                       learners, # need to be mlr3 objects!
                                       learners.names,
                                       num_splits           = 100,
                                       Z_CLAN               = Z,
                                       X1_BLP               = setup_X1(),
                                       X1_GATES             = setup_X1(),
                                       HT                   = FALSE,
                                       vcov_BLP             = setup_vcov(),
                                       vcov_GATES           = setup_vcov(),
                                       equal_variances_CLAN = FALSE,
                                       prop_aux             = 0.5,
                                       quantile_cutoffs     = c(0.25, 0.5, 0.75),
                                       diff_GATES           = setup_diff(),
                                       diff_CLAN            = setup_diff(),
                                       significance_level   = 0.05,
                                       min_variation        = 1e-05,
                                       parallel             = TrueIfUnix(),
                                       num_cores            = parallel::detectCores(),
                                       seed                 = NULL,
                                       store_learners       = FALSE,
                                       store_splits         = FALSE){

  # call correct main function
  do.call(what = get(ifelse(parallel,
                            "generic.ml.across.learners_parallel",
                            "generic.ml.across.learners_serial")),
          args = list(Z = Z, D = D, Y = Y,
                      propensity_scores          = propensity_scores,
                      learners                   = learners,
                      learners.names             = learners.names,
                      num_splits                 = num_splits,
                      Z_CLAN                     = Z_CLAN,
                      X1_BLP                     = X1_BLP,
                      X1_GATES                   = X1_GATES,
                      HT                         = HT,
                      vcov_BLP                   = vcov_BLP,
                      vcov_GATES                 = vcov_GATES,
                      equal_variances_CLAN       = equal_variances_CLAN,
                      prop_aux                   = prop_aux,
                      quantile_cutoffs           = quantile_cutoffs,
                      diff_GATES                 = diff_GATES,
                      diff_CLAN                  = diff_CLAN,
                      significance_level         = significance_level,
                      min_variation              = min_variation,
                      num_cores                  = num_cores,
                      seed                       = seed,
                      store_learners             = store_learners,
                      store_splits               = store_splits))

} # FUN



generic.ml.across.learners_serial <- function(Z, D, Y,
                                              propensity_scores,
                                              learners, # need to be mlr3 objects!
                                              learners.names,
                                              num_splits           = 100,
                                              Z_CLAN               = Z,
                                              X1_BLP               = setup_X1(),
                                              X1_GATES             = setup_X1(),
                                              HT                   = FALSE,
                                              vcov_BLP             = setup_vcov(),
                                              vcov_GATES           = setup_vcov(),
                                              equal_variances_CLAN = FALSE,
                                              prop_aux             = 0.5,
                                              quantile_cutoffs     = c(0.25, 0.5, 0.75),
                                              diff_GATES           = setup_diff(),
                                              diff_CLAN            = setup_diff(),
                                              significance_level   = 0.05,
                                              min_variation        = 1e-05,
                                              num_cores            = parallel::detectCores(), # dead argument here
                                              seed                 = NULL,
                                              store_learners       = FALSE,
                                              store_splits         = FALSE){


  # set seed
  if(is.null(seed)){

    # ensure reproducibility
    rng <- RNGkind()
    RNGkind("L'Ecuyer-CMRG")
    on.exit(RNGkind(kind = rng[1], normal.kind = rng[2], sample.kind = rng[3]))

  } else{

    # ensure reproducibility
    rng <- RNGkind()
    set.seed(seed, "L'Ecuyer")
    on.exit(RNGkind(kind = rng[1], normal.kind = rng[2], sample.kind = rng[3]))

  } # IF

  num.vars.in.Z_CLAN <- ncol(Z_CLAN)
  genericML.by.split <- list()
  N     <- length(Y)
  N_set <- 1:N
  prop  <- floor(prop_aux * N)


  # initialize
  if(store_splits) splits.mat <- matrix(NA_character_, N, num_splits,
                                        dimnames = list(NULL, paste0("split_", 1:num_splits)))

  # initialize
  generic_targets <- initializer.for.splits(Z = Z, Z_CLAN = Z_CLAN,
                                            learners = learners.names, num_splits = num_splits,
                                            quantile_cutoffs = quantile_cutoffs,
                                            diff_GATES = diff_GATES,
                                            diff_CLAN = diff_CLAN)

  # loop over the sample splits
  for(s in 1:num_splits){

    # perform sample splitting into main set and auxiliary set
    split.ls <- sample_split(D = D, N = N, N_set = N_set, prop = prop)
    M_set    <- split.ls$M_set
    A_set    <- split.ls$A_set

    if(store_splits){

      splits.mat[M_set, s] <- "M"
      splits.mat[A_set, s] <- "A"

    } # IF


    # loop over the learners
    for(i in 1:length(learners)){

      generic.ml.obj <-
        GenericML_single_NoChecks(Z = Z, D = D, Y = Y,
                                  propensity_scores = propensity_scores,
                                  learner = learners[[i]],
                                  M_set = M_set, A_set = A_set,
                                  Z_CLAN                       = Z_CLAN,
                                  X1_BLP                       = X1_BLP,
                                  X1_GATES                     = X1_GATES,
                                  HT                           = HT,
                                  vcov_BLP                     = vcov_BLP,
                                  vcov_GATES                   = vcov_GATES,
                                  equal_variances_CLAN         = equal_variances_CLAN,
                                  quantile_cutoffs             = quantile_cutoffs,
                                  diff_GATES                   = diff_GATES,
                                  diff_CLAN                    = diff_CLAN,
                                  significance_level           = significance_level,
                                  min_variation                = min_variation)

      generic_targets[[i]]$BLP[,,s]   <- generic.ml.obj$BLP$generic_targets
      generic_targets[[i]]$GATES[,,s] <- generic.ml.obj$GATES$generic_targets
      generic_targets[[i]]$best[,,s]  <- c(generic.ml.obj$best$lambda, generic.ml.obj$best$lambda.bar)

      if(store_learners){

        genericML.by.split[[learners.names[i]]][[s]] <- generic.ml.obj

      }

      for(j in 1:num.vars.in.Z_CLAN){
        generic_targets[[i]]$CLAN[[j]][,,s] <- generic.ml.obj$CLAN$generic_targets[[j]]
      }

    } # FOR learners
  } # FOR num_splits

  if(!store_learners) genericML.by.split <- NULL
  if(!store_splits)   splits.mat <- NULL

  return(list(generic_targets = generic_targets,
              genericML.by.split = genericML.by.split,
              splits = splits.mat))

} # END FUN




generic.ml.across.learners_parallel <- function(Z, D, Y,
                                              propensity_scores,
                                              learners, # need to be mlr3 objects!
                                              learners.names,
                                              num_splits           = 100,
                                              Z_CLAN               = Z,
                                              X1_BLP               = setup_X1(),
                                              X1_GATES             = setup_X1(),
                                              HT                   = FALSE,
                                              vcov_BLP             = setup_vcov(),
                                              vcov_GATES           = setup_vcov(),
                                              equal_variances_CLAN = FALSE,
                                              prop_aux             = 0.5,
                                              quantile_cutoffs     = c(0.25, 0.5, 0.75),
                                              diff_GATES           = setup_diff(),
                                              diff_CLAN            = setup_diff(),
                                              significance_level   = 0.05,
                                              min_variation        = 1e-05,
                                              num_cores            = parallel::detectCores(),
                                              seed                 = NULL,
                                              store_learners       = FALSE,
                                              store_splits         = FALSE){


  # set seed
  if(is.null(seed)){

    # ensure reproducibility
    rng <- RNGkind()
    RNGkind("L'Ecuyer-CMRG")
    on.exit(RNGkind(kind = rng[1], normal.kind = rng[2], sample.kind = rng[3]))

  } else{

    # ensure reproducibility
    rng <- RNGkind()
    set.seed(seed, "L'Ecuyer")
    on.exit(RNGkind(kind = rng[1], normal.kind = rng[2], sample.kind = rng[3]))

  } # IF

  num.vars.in.Z_CLAN <- ncol(Z_CLAN)
  K <- length(quantile_cutoffs) + 1
  num.generic_targets.gates <- K + length(diff_GATES$subtracted)
  num.generic_targets.clan <- K + length(diff_CLAN$subtracted)
  N     <- length(Y)
  N_set <- 1:N
  num.learners <- length(learners)
  prop <- floor(prop_aux * N)


  # loop over the sample splits
  out <- parallel::mclapply(mc.cores = num_cores, X = 1:num_splits, FUN = function(s){

    # initialize
    blp.3d     <- get_blp.3d(num.learners, learners.names)
    gates.3d   <- get_gates.3d(num.learners, learners.names, num.generic_targets.gates)
    best.3d    <- get_best.3d(num.learners, learners.names)
    clan.3d.ls <- get_clan.3d.ls(num.learners, learners.names, num.generic_targets.clan, num.vars.in.Z_CLAN, colnames(Z_CLAN))


    if(store_learners){
      genericML.by.split.ls <- rep(list(NA_real_), num.learners)
    } else{
      genericML.by.split.ls <- NULL
    } # IF


    # perform sample splitting into main set and auxiliary set
    split.ls <- sample_split(D = D, N = N, N_set = N_set, prop = prop)
    M_set    <- split.ls$M_set
    A_set    <- split.ls$A_set

    if(store_splits){

      splits.vec        <- rep("M", N)
      splits.vec[A_set] <- "A"

    } else{
      splits.vec <- NULL
    } # IF


    # loop over the learners
    for(i in 1:length(learners)){

      generic.ml.obj <-
        GenericML_single_NoChecks(Z = Z, D = D, Y = Y,
                                  propensity_scores = propensity_scores,
                                  learner = learners[[i]],
                                  M_set = M_set, A_set = A_set,
                                  Z_CLAN                       = Z_CLAN,
                                  X1_BLP                       = X1_BLP,
                                  X1_GATES                     = X1_GATES,
                                  HT                           = HT,
                                  vcov_BLP                     = vcov_BLP,
                                  vcov_GATES                   = vcov_GATES,
                                  equal_variances_CLAN         = equal_variances_CLAN,
                                  quantile_cutoffs             = quantile_cutoffs,
                                  diff_GATES                   = diff_GATES,
                                  diff_CLAN                    = diff_CLAN,
                                  significance_level           = significance_level,
                                  min_variation                = min_variation)

      blp.3d[,,i]   <- generic.ml.obj$BLP$generic_targets
      gates.3d[,,i] <- generic.ml.obj$GATES$generic_targets
      best.3d[,,i]  <- c(generic.ml.obj$best$lambda, generic.ml.obj$best$lambda.bar)

      for(j in 1:num.vars.in.Z_CLAN){
        clan.3d.ls[[j]][,,i] <- generic.ml.obj$CLAN$generic_targets[[j]]
      }

      if(store_learners) genericML.by.split.ls[[i]] <- generic.ml.obj

    } # FOR learners

    # return
    return(
    list(BLP = blp.3d, GATES = gates.3d, best = best.3d,
         CLAN = clan.3d.ls, GML = genericML.by.split.ls,
         splits = splits.vec))

  }) # MCLAPPLY


  # now bring it in desired form:
  generic_targets <- initializer.for.splits(Z = Z, Z_CLAN = Z_CLAN,
                                            learners = learners.names, num_splits = num_splits,
                                            quantile_cutoffs = quantile_cutoffs,
                                            diff_GATES = diff_GATES,
                                            diff_CLAN = diff_CLAN)

  for(s in 1:num_splits){
    for(i in 1:length(learners)){

      generic_targets[[i]]$BLP[,,s]   <- out[[s]]$BLP[,,i]
      generic_targets[[i]]$GATES[,,s] <- out[[s]]$GATES[,,i]
      generic_targets[[i]]$best[,,s]  <- out[[s]]$best[,,i]

      for(j in 1:num.vars.in.Z_CLAN){
        generic_targets[[i]]$CLAN[[j]][,,s] <- out[[s]]$CLAN[[j]][,,i]
      } # FOR

    } # FOR
  } # FOR


  if(store_splits){
    splits <- sapply(1:num_splits, function(s) out[[s]]$splits)
    colnames(splits) <- paste0("split_", 1:num_splits)
  } else{
    splits <- NULL
  } # IF

  if(store_learners){

    stored.learners <- list()

    for(s in 1:num_splits){
      for(i in 1:length(learners)){
        stored.learners[[learners.names[i]]][[s]] <- out[[s]]$GML[[i]]
      }
    }

  } else{
    stored.learners <- NULL
  } # IF

  # return
  return(list(generic_targets = generic_targets,
              genericML.by.split = stored.learners,
              splits = splits))

} # FUN
