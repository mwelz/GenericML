#' Makes plots for generic ML
#'
#' @param x An instance of \code{\link{GenericML}}.
#' @param learner The learner whose results are to be returned. Default is \code{"best"} for the best learner as measured by the lambda parameters.
#' @param type The analysis to be plotted. Either \code{"GATES"}, \code{"BLP"}, or \code{"CLAN"}. Default is \code{"GATES"}.
#' @param CLAN_variable Name of CLAN variable to be plotted. Only applicable if \code{type = "CLAN"}.
#' @param groups Character vector for the groups to be plotted for GATES and CLAN. Default is \code{"all"} for all groups. If there are \eqn{K} groups, this variable is a subset of \code{c("G1", "G2",...,"GK", "GK-G1", "GK-G2",...)}, but this set depends on the choices of the arguments \code{"diff_GATES"} and \code{"diff_CLAN"} of the \code{\link{GenericML}} function that generated the instance \code{x}.
#' @param limits The limits of the y-axis of the plot.
#' @param title The title of the plot.
#' @param ... Additional arguments.
#'
#' @import ggplot2
#'
#' @export
plot.GenericML <- function(x,
                           learner = "best",
                           type = "GATES",
                           CLAN_variable = NULL,
                           groups = "all",
                           limits = NULL,
                           title = NULL,
                           ...){

  # for better readability
  GenericML.obj <- x

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
    if(!CLAN_variable %in% names(GenericML.obj$VEIN$best_learners$CLAN)) stop("This variable was not used for CLAN.")

  } # IF

  if(!(type %in% c("BLP", "CLAN", "GATES"))) stop("Type needs to be either 'BLP', 'CLAN', or 'GATES'.")

  if(learner == "best" & type != "CLAN"){

    data <- GenericML.obj$VEIN$best_learners[[type]]


  } else if(learner == "best" & type == "CLAN"){

    data <- GenericML.obj$VEIN$best_learners[[type]][[CLAN_variable]]

  } else if(!(learner %in% GenericML.obj$arguments$learners_GenericML)){

    stop("This learner is not used in the generic ML procedure.")

  } else if(learner != "best" & type == "CLAN" ){

    data <- GenericML.obj$VEIN$all_learners[[type]][[learner]][[CLAN_variable]]

  } else{

    data <- GenericML.obj$VEIN$all_learners[[type]][[learner]]

  } # IF


  if(learner == "best"){
    data.blp <- GenericML.obj$VEIN$best_learners$BLP
  } else{
    data.blp <- GenericML.obj$VEIN$all_learners$BLP[[learner]]
  } # IF


  if(is.null(limits) & type != "BLP"){

    limits <- c(min(c(data[, "CB lower"], data.blp["beta.1", "CB lower"])),
                max(c(data[, "CB upper"], data.blp["beta.1", "CB upper"])))

  } else if(is.null(limits) & type == "BLP"){

    limits <- c(min(data[, "CB lower"]),
                max(data[, "CB upper"]))

  } # IF

  # adjusted confidence level
  confidence.level <- 1 - 2 * GenericML.obj$arguments$significance_level


  ## 1.1 make plot for GATES or CLAN ----
  # initialize variables (ugly hack for this R CMD check issue: https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when)
  group <- point.estimate <- ci.lower <- ci.upper <- NULL

  if(type != "BLP"){

    if(type == "GATES"){
      group.nam <- gsub("gamma.", "G", rownames(data))
    } else{
      group.nam <- gsub("delta.", "G", rownames(data))
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
