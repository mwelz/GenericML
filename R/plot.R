#' Set up information for a \code{GenericML()} plot
#'
#' Extract the relevant information for visualizing the point and interval estimates of the generic targets of interest. The generic targets of interest can be (subsets of) the parameters of the BLP, GATES, or CLAN analysis.
#'
#' @param x An object of the class \code{"\link{GenericML}"}, as returned by the function \code{\link{GenericML}()}.
#' @param type The analysis whose parameters shall be plotted. Either \code{"GATES"}, \code{"BLP"}, or \code{"CLAN"}. Default is \code{"GATES"}.
#' @param learner The learner whose results are to be returned. Default is \code{"best"} for the best learner as measured by the \eqn{\Lambda} parameters.
#' @param CLAN_variable Name of the CLAN variable to be plotted. Only applicable if \code{type = "CLAN"}.
#' @param groups Character vector indicating the per-group parameter estimates that shall be plotted in GATES and CLAN analyses. Default is \code{"all"} for all parameters. If there are \eqn{K} groups, this variable is a subset of \code{c("G1", "G2",...,"GK", "G1-G2", "G1-G2",..., "G1-GK", "GK-G1", "GK-G2",...)}, where Gk denotes the k-th group. Note that this set depends on the choices of the arguments \code{"diff_GATES"} and \code{"diff_CLAN"} of the \code{"\link{GenericML}"} object.
#'
#'
#' @details
#' This function is used internally by \code{\link{plot.GenericML}()}. It may also be useful for users who want to produce a similar plot, but who want more control over what information to display or how to display that information.
#'
#' @return
#' An object of class \code{"setup_plot"}, which is a list with the following elements.
#' \describe{
#'   \item{\code{data_plot}}{A data frame containing point and interval estimates of the generic target specified in the argument \code{type}.}
#'   \item{\code{data_BLP}}{A data frame containing point and interval estimates of the BLP analysis.}
#'   \item{\code{confidence_level}}{The confidence level of the confidence intervals. The confidence level is equal to  \code{1 - 2 * significance_level}, which is the adjustment proposed in the paper.}}
#'
#' @examples
#' if(require("ranger") && require("ggplot2")) {
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
#' learners <- c("random_forest")
#'
#' ## perform generic ML inference
#' # small number of splits to keep computation time low
#' x <- GenericML(Z, D, Y, learners,
#'                num_splits = 2,
#'                parallel = FALSE)
#'
#' ## the plot we wish to replicate
#' plot(x = x, type = "GATES")
#'
#' ## get the data to plot the GATES estimates
#' data <- setup_plot(x = x, type = "GATES")
#'
#' ## define variables to appease the R CMD check
#' group <- estimate <- ci_lower <- ci_upper <- NULL
#'
#' ## replicate the plot(x, type = "GATES")
#' # for simplicity, we skip aligning the colors
#' ggplot(mapping = aes(x = group,
#'                      y = estimate), data = data$data_plot) +
#'   geom_hline(aes(yintercept = 0),
#'              color = "black", linetype = "dotted") +
#'   geom_hline(aes(yintercept = data$data_BLP["beta.1", "estimate"],
#'                  color = "ATE"),
#'              linetype = "dashed") +
#'   geom_hline(aes(yintercept = data$data_BLP["beta.1", "ci_lower"],
#'                  color = paste0(100*data$confidence_level, "% CI (ATE)")),
#'              linetype = "dashed")  +
#'   geom_hline(yintercept = data$data_BLP["beta.1", "ci_upper"],
#'              linetype = "dashed", color = "red") +
#'   geom_point(aes(color = paste0("GATES with ",  100*data$confidence_level, "% CI")), size = 3) +
#'   geom_errorbar(mapping = aes(ymin = ci_lower,
#'                               ymax = ci_upper))
#' }
#'
#' @seealso
#' \code{\link{plot.GenericML}()}
#'
#' @export
setup_plot <- function(x,
                       type = "GATES",
                       learner = "best",
                       CLAN_variable = NULL,
                       groups = "all")
{

  ## 0.1 input check ----
  if(!inherits(x = x, what = "GenericML", which = FALSE)){
    stop("The object 'x' must be an instance of GenericML()")
  }

  if(type == "CLAN"){

    if(is.null(CLAN_variable)) stop("No CLAN variable specified")
    if(!CLAN_variable %in% names(x$VEIN$best_learners$CLAN)) stop("This variable was not used for CLAN.")

  } # IF

  if(!(type %in% c("BLP", "CLAN", "GATES"))) stop("Type needs to be either 'BLP', 'CLAN', or 'GATES'.")

  ## 0.2 prepare data ----
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
    df_blp <- as.data.frame(x$VEIN$best_learners$BLP)
  } else{
    df_blp <- as.data.frame(x$VEIN$all_learners$BLP[[learner]])
  } # IF

  # subset data frame of BLP to relevant information
  df_blp <- df_blp[,c("Estimate", "CB lower", "CB upper")]
  colnames(df_blp) <- c("estimate", "ci_lower", "ci_upper")

  # adjusted confidence level
  confidence_level <- 1.0 - 2.0 * x$arguments$significance_level

  if(type != "BLP"){

    ## 1.1 data for GATES of CLAN plot ----

    if(type == "GATES"){
      group.nam <- gsub("gamma.", "G", rownames(data))
    } else{
      group.nam <- gsub("delta.", "G", rownames(data))
    } # IF

    # create data frame for ggplot
    df <- data.frame(estimate = data[, "Estimate"],
                     ci_lower = data[, "CB lower"],
                     ci_upper = data[, "CB upper"],
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

  } else{

    ## 1.2 data for BLP plot ----
    df <- data.frame(estimate = data[, "Estimate"],
                     ci_lower = data[, "CB lower"],
                     ci_upper = data[, "CB upper"],
                     group = c("beta.1", "beta.2"))

  } # IF

  ## return
  return(
    structure(list(data_plot = df,
                   data_BLP = df_blp,
                   confidence_level = confidence_level),
              class = "setup_plot")
  )

} # FUN




#' Plot method for a \code{"GenericML"} object
#'
#' Visualizes the estimates of the generic targets of interest: plots the point estimates as well as the corresponding confidence intervals. The generic targets of interest can be (subsets of) the parameters of the BLP, GATES, or CLAN analysis.
#'
#' @param x An object of the class \code{"\link{GenericML}"}, as returned by the function \code{\link{GenericML}()}.
#' @param type The analysis whose parameters shall be plotted. Either \code{"GATES"}, \code{"BLP"}, or \code{"CLAN"}. Default is \code{"GATES"}.
#' @param learner The learner whose results are to be returned. Default is \code{"best"} for the best learner as measured by the \eqn{\Lambda} parameters.
#' @param CLAN_variable Name of the CLAN variable to be plotted. Only applicable if \code{type = "CLAN"}.
#' @param groups Character vector indicating the per-group parameter estimates that shall be plotted in GATES and CLAN analyses. Default is \code{"all"} for all parameters. If there are \eqn{K} groups, this variable is a subset of \code{c("G1", "G2",...,"GK", "G1-G2", "G1-G2",..., "G1-GK", "GK-G1", "GK-G2",...)}, where Gk denotes the k-th group. Note that this set depends on the choices of the arguments \code{"diff_GATES"} and \code{"diff_CLAN"} of the \code{"\link{GenericML}"} object.
#' @param ATE Logical. If \code{TRUE} (default), then the BLP estimate of the average treatment effect along with confidence intervals will be added to the plot. Only applicable if \code{type} is \code{"CLAN"} or \code{"GATES"}.
#' @param limits A numeric vector of length two holding the limits of the y-axis of the plot.
#' @param title The title of the plot.
#' @param ... Additional arguments to be passed down.
#'
#' @return
#' An object of class \code{"\link[ggplot2]{ggplot}"}.
#'
#' @details
#' If you wish to retrieve the data frame that this plot method visualizes, please use \code{\link{setup_plot}()}.
#'
#' @seealso
#' \code{\link{setup_plot}()},
#' \code{\link{GenericML}()},
#' \code{\link{get_BLP}()},
#' \code{\link{get_GATES}()},
#' \code{\link{get_CLAN}()},
#' \code{\link{setup_diff}()}
#'
#' @import ggplot2
#'
#' @examples
#' if(require("ranger")) {
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
#' learners <- c("random_forest")
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
                           ATE = TRUE,
                           limits = NULL,
                           title = NULL,
                           ...){

  ## input check
  stopifnot(is.logical(ATE) & length(ATE) == 1L)

  ## get the data to be plotted
  data <- setup_plot(x = x,
                     type = type,
                     learner = learner,
                     CLAN_variable = CLAN_variable,
                     groups = groups)
  df <- data$data_plot
  df_blp <- data$data_BLP
  confidence_level <- data$confidence_level


  ## specify the title
  if(is.null(title)){

    if(type == "CLAN"){
      title <- paste0("VEIN of CLAN for variable '", CLAN_variable, "'")
    } else{
      title <- paste0("VEIN of ", type)
    }

  } # IF


  ## specify the limits
  if(is.null(limits) & type != "BLP"){

    if(ATE){

      limits <- c(min(c(0.0, df[, "ci_lower"], df_blp["beta.1", "ci_lower"])),
                  max(c(0.0, df[, "ci_upper"], df_blp["beta.1", "ci_upper"])))

    } else{

      limits <- c(min(c(0.0, df[, "ci_lower"])),
                  max(c(0.0, df[, "ci_upper"])))

    } # IF

  } else if(is.null(limits) & type == "BLP"){

    limits <- c(min(c(0.0, df[, "ci_lower"])),
                max(c(0.0, df[, "ci_upper"])))

  } # IF


  ### make plot for GATES or CLAN
  # initialize variables (ugly hack for this R CMD check issue: https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when)
  group <- estimate <- ci_lower <- ci_upper <- NULL

  if(type != "BLP"){

    ## case 1: GATES or BLP

    if(type == "GATES"){
      ylab <- "Treatment Effect"
    } else{
      ylab <- paste0("Value of '", CLAN_variable, "'")
    } # IF

    # prepare the plot
    p <- ggplot(mapping = aes(x = group,
                              y = estimate), data = df) +
      geom_point(aes(color = paste0(type, " with ",  100*confidence_level, "% CI")), size = 3) +
      geom_errorbar(mapping = aes(ymin = ci_lower,
                                  ymax = ci_upper)) +
      geom_hline(aes(yintercept = 0),
                 color = "black", linetype = "dotted") +
      theme_light() +
      ylab(ylab) +
      xlab("Group by HTE Score") +
      theme(legend.title = element_blank(),
            legend.position = "bottom") +
      guides(color = guide_legend(override.aes = list(
        linetype = 0, size = 4, shape = 15, alpha = 1))) +
      ylim(limits[1], limits[2]) +
      ggtitle(title)


    if(ATE){

      p <- p +
        geom_hline(aes(yintercept = df_blp["beta.1", "estimate"],
                       color = "ATE"),
                   linetype = "dashed") +
        geom_hline(aes(yintercept = df_blp["beta.1", "ci_lower"],
                       color = paste0(100*confidence_level, "% CI (ATE)")),
                   linetype = "dashed")  +
        geom_hline(yintercept = df_blp["beta.1", "ci_upper"],
                   linetype = "dashed", color = "red") +
        scale_colour_manual(values = c("red","blue", "black"))

    } else{
      p <- p + scale_colour_manual(values = "black")
    } # IF


  } else{

    ## case 2: BLP
    p <- ggplot(mapping = aes(x = group,
                              y = estimate), data = df) +
      geom_hline(aes(yintercept = 0),
                     color = "black", linetype = "dotted") +
      geom_point(size = 3) +
      geom_errorbar(mapping = aes(ymin = ci_lower,
                                  ymax = ci_upper)) +
      theme_light() +
      ylab("Treatment Effect") +
      xlab(paste0(type, " with ",  100*confidence_level, "% CI")) +
      ylim(limits[1], limits[2]) +
      ggtitle(title) +
      scale_x_discrete(breaks = c("beta.1", "beta.2"),
                       labels = c(expression(paste(beta[1], " (ATE)")),
                                  expression(paste(beta[2], " (HTE)"))))

  } # IF

  return(p)

} # FUN


#' @export
plot.BLP_info <- function(x, ...){
  p <- x$plot
  if(is.null(p)){
    stop(paste0("This object does not contain a ggplot object.",
        " Try the argument plot = TRUE in get_BLP()."))
  } else{
    return(p)
  }
} # FUN


#' @export
plot.GATES_info <- function(x, ...){
  p <- x$plot
  if(is.null(p)){
    stop(paste0("This object does not contain a ggplot object.",
                " Try the argument plot = TRUE in get_GATES()."))
  } else{
    return(p)
  }
} # FUN


#' @export
plot.CLAN_info <- function(x, ...){
  p <- x$plot
  if(is.null(p)){
    stop(paste0("This object does not contain a ggplot object.",
                " Try the argument plot = TRUE in get_CLAN()."))
  } else{
    return(p)
  }
} # FUN
