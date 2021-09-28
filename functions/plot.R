#' Makes plot for generic ML
#' 
#' @param GenericML.obj an object as returned by genericML()
#' @param learner the learner whose results are to be returned. Default is 'best'
#' @param type the analysis to be plotted. Either 'GATES', 'BLP', or 'CLAN'. Default is 'GATES'.
#' @param CLAN.variable CLAN variable to be plotted. Only applicable if type = 'CLAN'.
#' @param groups.to.be.plotted The groups to be plotted for GATES and CLAN. Default is 'all'. If there are K groups, this variable can be set to a subset of {"G1", "G2",...,"GK", "GK-G1", "GK-G2",...}, but this set depends on the choices of the arguments 'differences.control_GATES' and 'differences.control_CLAN' of the GenericML() function that generated GenericML.obj.
#' @param limits the limits of the y-axis of the plot.
#' @param title the title of the plot.
plot.GenericML <- function(GenericML.obj,
                           learner = "best",
                           type = "GATES",
                           CLAN.variable = NULL,
                           groups.to.be.plotted = "all",
                           limits = NULL,
                           title = NULL){
  
  ## 0. input check ----
  if(type == "CLAN"){
    
    if(is.null(CLAN.variable)) stop("No CLAN variable specified")
    if(!CLAN.variable %in% names(GenericML.obj$VEIN$best.learners$CLAN)) stop("This variable was not used for CLAN.")
    
  } # IF
  
  if(!(type %in% c("BLP", "CLAN", "GATES"))) stop("Type needs to be either 'BLP', 'CLAN', or 'GATES'.")
  
  if(learner == "best" & type != "CLAN"){
    
    data <- GenericML.obj$VEIN$best.learners[[type]]
    
    
  } else if(learner == "best" & type == "CLAN"){
    
    data <- GenericML.obj$VEIN$best.learners[[type]][[CLAN.variable]]
    
  } else if(!(learner %in% GenericML.obj$arguments$learners.genericML)){
    
    stop("This learner is not used in the generic ML procedure.")
    
  } else if(learner != "best" & type == "CLAN" ){
    
    data <- GenericML.obj$VEIN$all.learners[[type]][[learner]][[CLAN.variable]]
    
  } else{
    
    data <- GenericML.obj$VEIN$all.learners[[type]][[learner]]
    
  } # IF
  
  
  if(learner == "best"){
    data.blp <- GenericML.obj$VEIN$best.learners$BLP
  } else{
    data.blp <- GenericML.obj$VEIN$all.learners$BLP[[learner]]
  } # IF
  
  
  if(is.null(limits) & type != "BLP"){
    
    limits <- c(min(c(data[, "CB lower"], data.blp["beta.1", "CB lower"])),
                max(c(data[, "CB upper"], data.blp["beta.1", "CB upper"])))
    
  } else if(is.null(limits) & type == "BLP"){
    
    limits <- c(min(data[, "CB lower"]),
                max(data[, "CB upper"]))
    
  } # IF
  
  # adjusted confidence level
  confidence.level <- 1 - 2 * GenericML.obj$arguments$significance.level
  
  
  ## 1.1 make plot for GATES or CLAN ----
  
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
    
    if(all(groups.to.be.plotted == "all")){
      
      groups.to.be.plotted_new <- group.nam
      
    } else if(!all(groups.to.be.plotted %in% group.nam)){
      
      stop(
      paste0("In 'groups.to.be.plotted', the input group(s) '", 
             paste0(groups.to.be.plotted[!groups.to.be.plotted %in% group.nam], collapse = "','"),
             "' is/are not a subset of c('", paste(group.nam, collapse = "','"), "')."))
      
    } else{
      groups.to.be.plotted_new <- groups.to.be.plotted
    } # IF
      

    # subset data frame for the plot
    df <- df[group.nam %in% groups.to.be.plotted_new,,drop = FALSE] 
    
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
