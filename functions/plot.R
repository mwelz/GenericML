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
