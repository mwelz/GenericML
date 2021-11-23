#' Plot method for a \code{GenericML} object
#'
#' Visualizes the estimates of the generic targets of interest: plots the point estimates as well as the corresponding confidence bounds. The generic targets of interest can be (subsets of) the parameters of the BLP, GATES, or CLAN analysis.
#'
#' @param x An instance of \code{\link{GenericML}}.
#' @param type The analysis whose parameters shall be plotted. Either \code{"GATES"}, \code{"BLP"}, or \code{"CLAN"}. Default is \code{"GATES"}.
#' @param learner The learner whose results are to be returned. Default is \code{"best"} for the best learner as measured by the \eqn{Lambda} parameters.
#' @param CLAN_variable Name of the CLAN variable to be plotted. Only applicable if \code{type = "CLAN"}.
#' @param groups Character vector indicating the per-group parameter estimates that shall be plotted in GATES and CLAN analyses. Default is \code{"all"} for all parameters. If there are \eqn{K} groups, this variable is a subset of \code{c("G1", "G2",...,"GK", "G1-G2", "G1-G2",..., "G1-GK", "GK-G1", "GK-G2",...)}, where Gk denotes the k-th group. Note that this set depends on the choices of the arguments \code{"diff_GATES"} and \code{"diff_CLAN"} of the \code{\link{GenericML}} object.
#' @param limits The limits of the y-axis of the plot.
#' @param title The title of the plot.
#' @param ... Additional arguments to be passed down.
#'
#' @return
#' An object of class \code{"ggplot"} (see \code{\link[ggplot2]{ggplot}}).
#'
#' @seealso
#' \code{\link{GenericML}},
#' \code{\link{get_BLP}},
#' \code{\link{get_GATES}},
#' \code{\link{get_CLAN}},
#' \code{\link{setup_diff}}
#'
#' @import ggplot2
#'
#' @examples
#' if(require("glmnet")) {
#'
#' ## generate data
#' set.seed(1)
#' n  <- 150                                  # number of observations
#' p  <- 5                                    # number of covariates
#' D  <- rbinom(n, 1, 0.5)                    # random treatment assignment
#' Z  <- matrix(runif(n*p), n, p)             # design matrix
#' Y0 <- as.numeric(Z %*% rexp(p) + rnorm(n)) # potential outcome without treatment
#' Y1 <- 2 + Y0                               # potential outcome under treatment
#' Y  <- ifelse(D == 1, Y1, Y0)               # observed outcome
#'
#' ## name the columns of Z
#' colnames(Z) <- paste0("V", 1:p)
#'
#' ## specify learners
#' learners <- c("lasso")
#'
#' ## specify quantile cutoffs (the 4 quartile groups here)
#' quantile_cutoffs <- c(0.25, 0.5, 0.75)
#'
#' ## specify the differenced generic targets of GATES and CLAN
#' diff_GATES <- setup_diff(subtract_from = "most",
#'                          subtracted = c(1,2,3))
#' diff_CLAN  <- setup_diff(subtract_from = "least",
#'                          subtracted = c(3,2))
#'
#' ## perform generic ML inference
#' # small number of splits to keep computation time low
#' x <- GenericML(Z, D, Y, learners, num_splits = 2,
#'                quantile_cutoffs = quantile_cutoffs,
#'                diff_GATES = diff_GATES,
#'                diff_CLAN = diff_CLAN,
#'                parallel = FALSE)
#'
#' ## plot BLP parameters
#' plot(x, type = "BLP")
#'
#' ## plot GATES parameters "G1", "G4", "G4-G1"
#' plot(x, type = "GATES", groups = c("G1", "G4", "G4-G1"))
#'
#' ## plot CLAN parameters "G1", "G2", "G2-G1" of variable "V1":
#' plot(x, type = "CLAN", CLAN_variable = "V1",
#'      groups = c("G1", "G2", "G1-G3"))
#' }
#'
#' @export
plot.GenericML <- function(x,
                           type = "GATES",
                           learner = "best",
                           CLAN_variable = NULL,
                           groups = "all",
                           limits = NULL,
                           title = NULL,
                           ...){

  if(class(x) != "GenericML"){
    stop("The object 'x' must be an instance of GenericML()")
  }

  # specify the title
  if(is.null(title)){

    if(type == "CLAN"){
      title <- paste0("VEIN of CLAN for variable '", CLAN_variable, "'")
    } else{
      title <- paste0("VEIN of ", type)
    }

  } # IF

  ## 0. input check ----
  if(type == "CLAN"){

    if(is.null(CLAN_variable)) stop("No CLAN variable specified")
    if(!CLAN_variable %in% names(x$VEIN$best_learners$CLAN)) stop("This variable was not used for CLAN.")

  } # IF

  if(!(type %in% c("BLP", "CLAN", "GATES"))) stop("Type needs to be either 'BLP', 'CLAN', or 'GATES'.")

  if(learner == "best" & type != "CLAN"){

    data <- x$VEIN$best_learners[[type]]


  } else if(learner == "best" & type == "CLAN"){

    data <- x$VEIN$best_learners[[type]][[CLAN_variable]]

  } else if(!(learner %in% x$arguments$learners_GenericML)){

    stop("This learner is not used in the generic ML procedure.")

  } else if(learner != "best" & type == "CLAN" ){

    data <- x$VEIN$all_learners[[type]][[learner]][[CLAN_variable]]

  } else{

    data <- x$VEIN$all_learners[[type]][[learner]]

  } # IF


  if(learner == "best"){
    data.blp <- x$VEIN$best_learners$BLP
  } else{
    data.blp <- x$VEIN$all_learners$BLP[[learner]]
  } # IF


  if(is.null(limits) & type != "BLP"){

    limits <- c(min(c(0, data[, "CB lower"], data.blp["beta.1", "CB lower"])),
                max(c(0, data[, "CB upper"], data.blp["beta.1", "CB upper"])))

  } else if(is.null(limits) & type == "BLP"){

    limits <- c(min(c(0, data[, "CB lower"])),
                max(c(0, data[, "CB upper"])))

  } # IF

  # adjusted confidence level
  confidence.level <- 1 - 2 * x$arguments$significance_level


  ## 1.1 make plot for GATES or CLAN ----
  # initialize variables (ugly hack for this R CMD check issue: https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when)
  group <- point.estimate <- ci.lower <- ci.upper <- NULL

  if(type != "BLP"){

    if(type == "GATES"){
      group.nam <- gsub("gamma.", "G", rownames(data))
      ylab      <- "Treatment Effect"
    } else{
      group.nam <- gsub("delta.", "G", rownames(data))
      ylab      <- paste0("Value of '", CLAN_variable, "'")
    } # IF

    # create data frame for ggplot
    df <- data.frame(point.estimate = data[, "Estimate"],
                     ci.lower = data[, "CB lower"],
                     ci.upper = data[, "CB upper"],
                     group = factor(group.nam, levels = group.nam))

    if(all(groups == "all")){

      groups_new <- group.nam

    } else if(!all(groups %in% group.nam)){

      stop(
      paste0("In 'groups', the input group(s) '",
             paste0(groups[!groups %in% group.nam], collapse = "','"),
             "' is/are not a subset of c('", paste(group.nam, collapse = "','"), "')."))

    } else{
      groups_new <- groups
    } # IF


    # subset data frame for the plot
    df <- df[group.nam %in% groups_new,,drop = FALSE]

    # make the plot
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
      ylab(ylab) +
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
                     group = c("beta.1", "beta.2"))

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
      ggtitle(title) +
      scale_x_discrete(breaks = c("beta.1", "beta.2"),
                       labels = c(expression(paste(beta[1], " (ATE)")),
                                  expression(paste(beta[2], " (HTE)"))))

    return(p)

  } # IF
} # FUN
