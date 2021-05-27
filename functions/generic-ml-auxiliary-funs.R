library(ggplot2)

#' Estimates the BLP parameters based on the main sample M. 
#' 
#' @param D a binary vector of treatment status of length _|M|_
#' @param Y a vector of responses of length _|M|_
#' @param proxy.baseline a vector of proxy baseline estimates of length _M_
#' @param proxy.cate a vector of proxy CATE estimates of length _M_
#' @param propensity.scores a vector of propensity scores of length _|M|_
#' @return BLP coefficients with inference statements
#' 
#' @export
#' 
#' TODO: implement same with HT transformation! 
get.BLP.params.classic <- function(D, Y, propensity.scores, 
                                   proxy.baseline, proxy.cate, 
                                   significance.level = 0.05){
  
  # prepare weights
  weights <- 1 / (propensity.scores * (1 - propensity.scores))
  
  # prepare covariate matrix
  X <- data.frame(B = proxy.baseline, 
                  S = proxy.cate,
                  beta.1 = D - propensity.scores, 
                  beta.2 = (D - propensity.scores) * (proxy.cate - mean(proxy.cate))) 
  
  # fit weighted linear regression by OLS
  blp.obj <- lm(Y ~., data = data.frame(Y, X), weights = weights)
  
  # extract coefficients
  coefficients     <- summary(blp.obj)$coefficients
  
  # inference on beta2: test the null that it is 1) = 0, 2) = 1.
  beta2.inference <- matrix(NA_real_, 2, 2)
  rownames(beta2.inference) <- c("H0: beta.2 = 0", "H0: beta.2 = 1")
  colnames(beta2.inference) <- c("t value", "Pr(>|t|)")
  beta2.inference[1,] <- coefficients["beta.2", c("t value", "Pr(>|t|)")]
  beta2.inference[2, "t value"] <- 
    (coefficients["beta.2", "Estimate"] - 1) / coefficients["beta.2", "Std. Error"]  
  beta2.inference[2, "Pr(>|t|)"] <- 
    2 * pt(abs(beta2.inference[2, "t value"]), df = blp.obj$df.residual, lower.tail = FALSE)
  
  # generic targets
  generic.targets <- coefficients[c("beta.1", "beta.2"), ]
  colnames(generic.targets) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  generic.targets[,"Pr(>|z|)"] <- 2 * pnorm(abs(generic.targets[,"z value"]), lower.tail = FALSE)
  ci.lo <- generic.targets[,"Estimate"] - qnorm(1-significance.level/2) * generic.targets[,"Std. Error"]
  ci.up <- generic.targets[,"Estimate"] + qnorm(1-significance.level/2) * generic.targets[,"Std. Error"]
  generic.targets <- cbind(generic.targets, ci.lo, ci.up)
  colnames(generic.targets) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "CB lower", "CB upper")
  generic.targets <- generic.targets[,c("Estimate", "CB lower", "CB upper", 
                                        "Std. Error", "z value", "Pr(>|z|)")]
  

  return(list(lm.obj = blp.obj, 
              blp.coefficients = blp.obj$coefficients[c("beta.1", "beta.2")],
              generic.targets = matrix(generic.targets, nrow = 1,
                                       dimnames = list("beta.2", names(generic.targets))),
              coefficients = coefficients,
              beta2.inference = beta2.inference))
  
} # END FUN


#' Estimates the GATES parameters based on the main sample M. 
#' 
#' @param D a binary vector of treatment status of length _|M|_
#' @param Y a vector of responses of length _|M|_
#' @param propensity.scores a vector of propensity scores of length _|M|_
#' @param proxy.baseline a vector of proxy baseline estimates of length _M_
#' @param proxy.cate a vector of proxy CATE estimates of length _M_
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate 
#' the group memberships (such a matrix is returned by the function quantile.group())
#' @return GATES coefficients 
#' 
#' @export
#' 
#' TODO: implement same with HT transformation! 
get.GATES.params.classic <- function(D, Y, 
                                     propensity.scores, 
                                     proxy.baseline, proxy.cate,
                                     group.membership.main.sample,
                                     significance.level = 0.05){
  
  # make the group membership a binary matrix
  groups <- 1 * group.membership.main.sample
  
  # prepare weights
  weights <- 1 / (propensity.scores * (1 - propensity.scores))
  
  # prepare covariate matrix
  X <- data.frame(B = proxy.baseline, 
                  S = proxy.cate,
                  (D - propensity.scores) * groups)
  colnames(X) <- c(colnames(X)[c(1,2)], paste0("gamma.", 1:ncol(groups)))
  
  # fit weighted linear regression by OLS
  gates.obj <- lm(Y ~., data = data.frame(Y, X), weights = weights)
  
  # extract coefficients
  coefficients                 <- summary(gates.obj)$coefficients
  gates.coefficients           <- coefficients[paste0("gamma.", 1:ncol(groups)), 1]
  gates.coefficients.quantiles <- colnames(groups)
  names(gates.coefficients.quantiles) <- paste0("gamma.", 1:ncol(groups))
  
  # prepare generic target parameters
  coefficients.temp <- coefficients[-c(1,2,3),]
  colnames(coefficients.temp) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  coefficients.temp[,"Pr(>|z|)"] <- 2 * pnorm(abs(coefficients.temp[,"z value"]), lower.tail = FALSE)
  ci.lo <- coefficients.temp[,"Estimate"] - qnorm(1-significance.level/2) * coefficients.temp[,"Std. Error"]
  ci.up <- coefficients.temp[,"Estimate"] + qnorm(1-significance.level/2) * coefficients.temp[,"Std. Error"]
  coefficients.temp <- cbind(coefficients.temp, ci.lo, ci.up)
  colnames(coefficients.temp) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "CB lower", "CB upper")
  coefficients.temp <- coefficients.temp[,c("Estimate", "CB lower", "CB upper", 
                                            "Std. Error", "z value", "Pr(>|z|)")]
  
  # prepare generic target parameters for the difference
  covmat  <- stats::vcov(gates.obj)
  diff    <- coefficients[paste0("gamma.", ncol(groups)), "Estimate"] - 
    coefficients["gamma.1", "Estimate"]
  diff.se <- sqrt(covmat[paste0("gamma.", ncol(groups)), paste0("gamma.", ncol(groups))] +
                    covmat["gamma.1", "gamma.1"] - 2 * covmat[paste0("gamma.", ncol(groups)), "gamma.1"])
  ci.lo   <- diff - qnorm(1-significance.level/2) * diff.se
  ci.up   <- diff + qnorm(1-significance.level/2) * diff.se
  zstat   <- diff / diff.se
  pval    <- 2 * pnorm(abs(zstat), lower.tail = FALSE)
  
  return(list(lm.obj = gates.obj, 
              gates.coefficients = gates.coefficients,
              gates.coefficients.quantiles = gates.coefficients.quantiles,
              generic.targets = rbind(coefficients.temp, 
                                      matrix(c(diff, ci.lo, ci.up, diff.se, zstat, pval), nrow = 1, 
                                             dimnames = list("gamma.K-gamma.1", NULL)) ),
              coefficients = coefficients))

} # END FUN


#' Estimates the CLAN parameters in the main sample
#' 
#' @param Z.clan.main.sample a matrix with _|M|_ rows. Each column represents a variable for which CLAN shall be performed.
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate 
#' the group memberships (such a matrix is returned by the function quantile.group())
#' @return The two CLAN parameters ("most" affected and "least" affected) for each variable in Z.clan.main.sample
#' 
#' @export
get.CLAN.parameters <- function(Z.clan.main.sample, 
                                group.membership.main.sample, 
                                significance.level = 0.05){
  
  K <- ncol(group.membership.main.sample)
  
  # initialize
  generic.targets   <- list()
  clan.coefficients <- matrix(NA_real_, nrow = 3, ncol = ncol(Z.clan.main.sample))
  
  # loop over the CLAN variables
  for(j in 1:ncol(Z.clan.main.sample)){
    
    # initialize matrix
    out.mat <- matrix(NA_real_, nrow = 3, ncol = 6)
    
    # get summary statistics for least affected group
    ttest.delta1 <- stats::t.test(Z.clan.main.sample[group.membership.main.sample[, 1], j])
    ci.lo        <- ttest.delta1$estimate - qnorm(1-significance.level/2) * ttest.delta1$stderr 
    ci.up        <- ttest.delta1$estimate + qnorm(1-significance.level/2) * ttest.delta1$stderr 
    pval         <- 2 * pnorm(abs(ttest.delta1$statistic), lower.tail = FALSE)
    out.mat[1,]  <- c(ttest.delta1$estimate, ci.lo, ci.up,
                      ttest.delta1$stderr, ttest.delta1$statistic, pval)
    
    # get summary statistics for most affected group
    ttest.deltaK <- stats::t.test(Z.clan.main.sample[group.membership.main.sample[, K], j])
    ci.lo        <- ttest.deltaK$estimate - qnorm(1-significance.level/2) * ttest.deltaK$stderr 
    ci.up        <- ttest.deltaK$estimate + qnorm(1-significance.level/2) * ttest.deltaK$stderr 
    pval         <- 2 * pnorm(abs(ttest.deltaK$statistic), lower.tail = FALSE)
    out.mat[2,]  <- c(ttest.deltaK$estimate, ci.lo, ci.up,
                      ttest.deltaK$stderr, ttest.deltaK$statistic, pval)
    
    # get summary statistics for difference between most and least affected group
    diff        <- ttest.deltaK$estimate - ttest.delta1$estimate
    diff.se     <- sqrt( ttest.deltaK$stderr^2 + ttest.delta1$stderr^2 )
    ci.lo       <- diff - qnorm(1-significance.level/2) * diff.se
    ci.up       <- diff + qnorm(1-significance.level/2) * diff.se
    diff.ttest  <- diff / diff.se
    pval        <- 2 * pnorm(abs(diff.ttest), lower.tail = FALSE)
    out.mat[3,] <- c(diff, ci.lo, ci.up, diff.se, diff.ttest, pval)
    
    colnames(out.mat)     <- c("Estimate", "CB lower", "CB upper", "Std. Error", "z value", "Pr(>|z|)")
    rownames(out.mat)     <- c("delta.1", "delta.K", "delta.K-delta.1")
    generic.targets[[j]]  <- out.mat
    clan.coefficients[,j] <- out.mat[,1] 

  } # END FOR
  
  names(generic.targets) <- colnames(clan.coefficients) <- colnames(Z.clan.main.sample)
  rownames(clan.coefficients) <- c("delta.1", "delta.K", "delta.K-delta.1")
  
  return(list(clan.coefficients = clan.coefficients,
              generic.targets   = generic.targets))
  
} # END FUN


#' returns the two parameters that are used to find the best ML method
#' 
#' @param BLP.obj an object as returned by get.CLAN.parameters()
#' @param GATES.obj an object as returned by get.BLP.parameters()
#' @param proxy.cate.main.sample Proxy CATE estimators for the main sample
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate 
#' the group memberships (such a matrix is returned by the function quantile.group())
#' @return lambda and lambda.bar parameters
#' 
#' @export
best.ml.method.parameters <- function(BLP.obj,
                                      GATES.obj, 
                                      proxy.cate.main.sample, 
                                      group.membership.main.sample){
  
  return(list(lambda = as.numeric(BLP.obj$blp.coefficients["beta.2"]^2 * var(proxy.cate.main.sample)),
              lambda.bar = as.numeric(colSums(group.membership.main.sample) %*%  GATES.obj$gates.coefficients^2)))
  
} # END FUN


make.mlr3.string <- function(learner.str, regr = TRUE){
  # helper function. Requires input of type 'mlr3::lrn("cv_glmnet", s = "lambda.min")' (note the absence of classif and regr)
  
  if(substr(learner.str, start = 1, stop = 6) != "mlr3::"){
    
    learner <- learner.str
    
  } else{
    
    learner <- paste0(substr(learner.str, start = 1, stop = 11), ifelse(regr, "regr.", "classif."), 
                      substr(learner.str, start = 12, stop = 1e8))
    
    learner <-eval(parse(text = learner))
    
  } # END IF
  
  return(learner)
  
} # END FUN


initializer.for.splits <- function(Z, Z.clan, learners,
                                   num.splits, quantile.cutoffs){
  
  # helper function that initializes object in the generic ML splitting procedure
  
  if(is.null(Z.clan)){
    d <- ncol(Z)
    Z.clan.nam <- colnames(Z)
  } else{
    d <- ncol(Z.clan)
    Z.clan.nam <- colnames(Z.clan)
  }
  
  clan <- array(NA_real_, dim = c(3, 6, num.splits), 
                dimnames = list(c("delta.1", "delta.K", "delta.K-delta.1"), 
                                c("Estimate", "CB lower", "CB upper", "Std. Error", "z value", "Pr(>|z|)"),
                                NULL))
  
  gates <- array(NA_real_, dim = c(length(quantile.cutoffs)+2, 6, num.splits),
                 dimnames = list(
                   c(paste0("gamma.", 1:(length(quantile.cutoffs)+1)), "gamma.K-gamma.1"),
                   c("Estimate", "CB lower", "CB upper", "Std. Error", "z value", "Pr(>|z|)"), 
                   NULL))
  
  blp <- array(NA_real_, dim = c(2, 6, num.splits),
               dimnames = list(c("beta.1", "beta.2"), 
                               c("Estimate", "CB lower", "CB upper", "Std. Error", "z value", "Pr(>|z|)"),
                               NULL))
  
  best <- array(NA_real_, dim = c(1, 2, num.splits), 
                dimnames = list(NULL, c("lambda", "lambda.bar"), NULL))
  
  clan.lists <- lapply(1:d, function(...) clan )
  names(clan.lists) <- Z.clan.nam
  
  out.ls <- lapply(1:length(learners), function(...){
    
    list(BLP = blp, GATES = gates, CLAN = clan.lists, best = best)
    
  })
  
  names(out.ls) <- learners
  return(out.ls)
  
} # END FUN


generic.ml.across.learners <- function(Z, D, Y, 
                                       propensity.scores, 
                                       learners, 
                                       num.splits = 50,
                                       Z.clan = NULL, 
                                       proportion.in.main.set = 0.5, 
                                       quantile.cutoffs = c(0.25, 0.5, 0.75),
                                       significance.level = 0.05, 
                                       store.learners = FALSE,
                                       store.splits = FALSE){
  
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
      
      generic.ml.obj <- 
        get.generic.ml.for.given.learner(Z = Z, D = D, Y = Y, 
                                         propensity.scores = propensity.scores, 
                                         learner = learners[i], 
                                         M.set = M.set, A.set = A.set,
                                         Z.clan = Z.clan, 
                                         proportion.in.main.set = proportion.in.main.set, 
                                         quantile.cutoffs = quantile.cutoffs,
                                         significance.level = significance.level)
      
      generic.targets[[i]]$BLP[,,s]   <- generic.ml.obj$BLP$generic.targets
      generic.targets[[i]]$GATES[,,s] <- generic.ml.obj$GATES$generic.targets
      generic.targets[[i]]$best[,,s]  <- c(generic.ml.obj$best$lambda, generic.ml.obj$best$lambda.bar)
      
      if(store.learners){
        
        genericML.by.split[[learners[i]]][[s]] <- generic.ml.obj

      }
      
      for(j in 1:num.vars.in.Z.clan){
        generic.targets[[i]]$CLAN[[j]][,,s] <- generic.ml.obj$CLAN$generic.targets[[j]]
      }
      
    } # FOR learners
  } # FOR num.splits
  
  if(!store.learners) genericML.by.split <- NULL
  if(!store.splits)   splits.mat <- NULL
  
  return(list(generic.targets = generic.targets, 
         genericML.by.split = genericML.by.split,
         splits = splits.mat)) # TODO: also extract A.set, M.set, generic.ml.obj
  
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


# get lower and upper median as in comment 4.2 in Chernozhukov et al. (2021)
Med <- function(x){
  
  # get the empirical CDF of X
  ecdf.x <- ecdf(x)
  
  # evaluate the ensuing probabilities
  F.x <- ecdf.x(x)
  
  # get lower median and upper median
  lower       <- min(x[F.x >= 0.5])
  upper.array <- x[(1 - F.x) >= 0.5]
  upper       <- ifelse(length(upper.array) == 0, 
                        lower,
                        max(upper.array)) # account for case where upper.array is empty
  
  return(list(lower.median = lower, 
              upper.median = upper, 
              Med = mean(c(lower, upper))))
  
} # END FUN



initialize.gen.ml <- function(generic.ml.across.learners.obj){
  
  # helper function for initialization of the final returned object
  
  gates.nam  <- rownames(generic.ml.across.learners.obj[[1]]$GATES[,,1])
  z.clan.nam <- names(generic.ml.across.learners.obj[[1]]$CLAN)
  clan.nam   <- rownames(generic.ml.across.learners.obj[[1]]$CLAN[[1]][,,1])
  learners   <- names(generic.ml.across.learners.obj)
  
  blp.mat <- matrix(NA_real_, nrow = 2, ncol = 5, 
                    dimnames = list(c("beta.1", "beta.2"), 
                                    c("Estimate", "CB lower", "CB upper", "p-value adjusted", "p-value raw")))
  
  gates.mat <- matrix(NA_real_, nrow = length(gates.nam), ncol = 5, 
                      dimnames = list(gates.nam, 
                                      c("Estimate", "CB lower", "CB upper", "p-value adjusted", "p-value raw")))
  
  clan.mat <- matrix(NA_real_, nrow = length(clan.nam), ncol = 5, 
                     dimnames = list(clan.nam, 
                                     c("Estimate", "CB lower", "CB upper", "p-value adjusted", "p-value raw")))
  
  
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


#' Makes plot for generic ML
#' 
#' @param genericML.obj an object as returned by genericML()
#' @param learner the learner whose results are to be returned. Default is 'best'
#' @param type the analysis to be plotted. Either 'GATES', 'BLP', or 'CLAN'. Default is 'GATES'.
#' @param CLAN.variable CLAN variable to be plotted. Only applicable if type = 'CLAN'.
#' @param limits the limits of the y-axis of the plot.
#' @param title the title of the plot.
genericML.plot <- function(genericML.obj,
                           learner = "best",
                           type = "GATES",
                           CLAN.variable = NULL,
                           limits = NULL,
                           title = NULL){
  
  ## 0. input check ----
  if(type == "CLAN"){
    
    if(is.null(CLAN.variable)) stop("No CLAN variable specified")
    if(!CLAN.variable %in% names(genericML.obj$VEIN$best.learners$CLAN)) stop("This variable was not used for CLAN.")
    
  } # IF
  
  if(!(type %in% c("BLP", "CLAN", "GATES"))) stop("Type needs to be either 'BLP', 'CLAN', or 'GATES'.")
  
  if(learner == "best" & type != "CLAN"){
    
    data <- genericML.obj$VEIN$best.learners[[type]]
    
    
  } else if(learner == "best" & type == "CLAN"){
    
    data <- genericML.obj$VEIN$best.learners[[type]][[CLAN.variable]]
    
  } else if(!(learner %in% genericML.obj$arguments$learners.genericML)){
    
    stop("This learner is not used in the generic ML procedure.")
    
  } else if(learner != "best" & type == "CLAN" ){
    
    data <- genericML.obj$VEIN$all.learners[[type]][[learner]][[CLAN.variable]]
    
  } else{
    
    data <- genericML.obj$VEIN$all.learners[[type]][[learner]]
    
  } # IF
  
  
  if(learner == "best"){
    data.blp <- genericML.obj$VEIN$best.learners$BLP
  } else{
    data.blp <- genericML.obj$VEIN$all.learners$BLP[[learner]]
  } # IF
  
  
  if(is.null(limits) & type != "BLP"){
    
    limits <- c(min(c(data[, "CB lower"], data.blp["beta.1", "CB lower"])),
                max(c(data[, "CB upper"], data.blp["beta.1", "CB upper"])))
    
  } else if(is.null(limits) & type == "BLP"){
    
    limits <- c(min(data[, "CB lower"]),
                max(data[, "CB upper"]))
    
  } # IF
  
  # adjusted confidence level
  confidence.level <- 1 - 2 * genericML.obj$arguments$significance.level
  
  
  ## 1.1 make plot for GATES or CLAN ----
  
  if(type != "BLP"){
    
    K  <- nrow(data) - 1
    df <- data.frame(point.estimate = data[, "Estimate"],
                     ci.lower = data[, "CB lower"],
                     ci.upper = data[, "CB upper"],
                     group = c(paste0("G", 1:K), paste0("G", K, "-G1")))
    
    p <- ggplot(mapping = aes(x = group,
                              y = point.estimate), data = df) +
      geom_hline(aes(yintercept = 0),
                 color = "black", linetype = "dotted") +
      geom_hline(aes(yintercept = data.blp["beta.1", "Estimate"],
                     color = "ATE"),
                 linetype = "dashed") +
      geom_hline(aes(yintercept = data.blp["beta.1", "CB lower"],
                     color = paste0(100*confidence.level, "% CB (ATE)")),
                 linetype = "dashed")  +
      geom_hline(yintercept = data.blp["beta.1", "CB upper"],
                 linetype = "dashed", color = "red") +
      geom_point(aes(color = paste0(type, " with ",  100*confidence.level, "% CB")), size = 3) +
      geom_errorbar(mapping = aes(ymin = ci.lower,
                                  ymax = ci.upper)) +
      theme_light() +
      ylab("Treatment Effect") +
      xlab("Group by HTE Score") +
      scale_colour_manual(values = c("red","blue", "black")) +
      theme(legend.title = element_blank(),
            legend.position = "bottom") + 
      guides(color = guide_legend(override.aes = list(
        linetype = 0, size = 4, shape = 15, alpha = 1))) +
      ylim(limits[1], limits[2]) +
      ggtitle(title)
    
    return(p)
    
  } else{
    
    ## 1.2 make plot for BLP ----
    df <- data.frame(point.estimate = data[, "Estimate"],
                     ci.lower = data[, "CB lower"],
                     ci.upper = data[, "CB upper"],
                     group = c("beta.1 (ATE)", "beta.2 (HTE)"))
    
    p <- ggplot(mapping = aes(x = group,
                              y = point.estimate), data = df) +
      geom_hline(aes(yintercept = 0),
                 color = "black", linetype = "dotted") +
      geom_point(size = 3) +
      geom_errorbar(mapping = aes(ymin = ci.lower,
                                  ymax = ci.upper)) +
      theme_light() +
      ylab("Treatment Effect") +
      xlab(paste0(type, " with ",  100*confidence.level, "% CB")) +
      ylim(limits[1], limits[2]) +
      ggtitle(title)
    
    return(p)
    
  } # IF
} # FUN


#' Partition a vector x into groups based on its quantiles
#' 
#' @param x the vector to be partitioned
#' @param cutoffs the quantile cutoffs for the partition. Default are the quartiles.
#' @param quantile.nam logical. Shall the cutoff values be included?
quantile.group <- function(x,
                           cutoffs = c(0.25, 0.5, 0.75),
                           quantile.nam = TRUE){
  # cutoffs are the quantile cutoffs (like c(0.25, 0.5, 0.75))
  q         <- quantile(x, cutoffs)
  q         <- c(-Inf, q, Inf)
  groups    <- as.character(cut(x, breaks = q, include.lowest = TRUE, right = FALSE, dig.lab = 3))
  group.nam <- unique(groups)
  group.nam <- group.nam[order(
    as.numeric(substr(sub("\\,.*", "", group.nam), 2, stop = 1e8L)), 
    decreasing = FALSE)] # ensure the order is correct
  group.mat <- matrix(NA, length(x), length(group.nam))
  nam       <- rep(NA, length(group.nam))
  
  for(j in 1:length(group.nam)){
    if(j == 1){
      nam[j] <- paste0("<", 100*cutoffs[j], "% quantile")
    } else if (j == length(group.nam)){
      nam[j] <- paste0(">=", 100*cutoffs[j-1], "% quantile")
    } else{
      nam[j] <- paste0("[", 100*cutoffs[j-1], ",", 100*cutoffs[j], ")% quantile")
    }
    group.mat[,j] <- groups == group.nam[j]
  }
  
  if(quantile.nam){
    colnames(group.mat) <- nam
  } else{
    colnames(group.mat) <- group.nam
  }
  return(group.mat)
} # FUN
