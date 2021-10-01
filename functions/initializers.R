initializer.for.splits <- function(Z, Z_CLAN, learners,
                                   num.splits, quantile.cutoffs,
                                   differences.control_GATES,
                                   differences.control_CLAN){
  
  # helper function that initializes object in the generic ML splitting procedure
  
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
  
  clan.lists <- lapply(1:ncol(Z_CLAN), function(...) clan )
  names(clan.lists) <- colnames(Z_CLAN)
  
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


# for parallelized across learner function:
get_blp.3d <- function(num.learners, learners.names){
  
  array(NA_real_, dim = c(2, 7, num.learners), dimnames = list(c("beta.1", "beta.2"), c("Estimate", "CB lower", "CB upper", "Std. Error", "z value", "Pr(<z)", "Pr(>z)"), learners.names))
  
} # FUN


get_gates.3d <- function(num.learners, learners.names, num.generic.targets.gates){
  
  array(NA_real_, dim = c(num.generic.targets.gates, 7, num.learners),
        dimnames = list(paste0("gamma.", 1:num.generic.targets.gates), c("Estimate", "CB lower", "CB upper", "Std. Error", "z value", "Pr(<z)", "Pr(>z)"), learners.names))
  
} # FUN


get_best.3d <- function(num.learners, learners.names){
  
  array(NA_real_, dim = c(1, 2, num.learners), 
        dimnames = list(NULL, c("lambda", "lambda.bar"), learners.names))
  
} # FUN


get_clan.3d.ls <- function(num.learners, learners.names, num.generic.targets.clan, num.vars.in.Z_CLAN, Z_CLAN.names){
  
  clan.3d.ls <- lapply(1:num.vars.in.Z_CLAN, 
                       function(...) array(NA_real_, dim = c(num.generic.targets.clan, 7, num.learners), dimnames = list(paste0("gamma.", 1:num.generic.targets.clan), c("Estimate", "CB lower", "CB upper", "Std. Error", "z value", "Pr(<z)", "Pr(>z)"), learners.names)) )
  names(clan.3d.ls) <- Z_CLAN.names
  return(clan.3d.ls)
  
} # FUN


