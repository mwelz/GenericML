# void functions for input checks

InputChecks_D <- function(D){
  
  # input checks
  if((!all(c(0, 1) %in% unique(D))) | (length(unique(D)) != 2)) stop("Treatment assignment D is not binary", call. = FALSE)
  if(!is.numeric(D)) stop("D is not numeric", call. = FALSE)
  
} # FUN


InputChecks_Y <- function(Y){
  
  if(!is.numeric(Y)) stop("Y is not numeric", call. = FALSE)
  if(!(is.matrix(Y) | is.data.frame(Y) | is.vector(Y))) stop("Y must be either a vector, data frame, or matrix", call. = FALSE)
  
  if(is.matrix(Y) | is.data.frame(Y)){
    
    if(ncol(Y) != 1) stop("The matrix/data frame Y must have only one column", call. = FALSE)
    
  } # IF

} # FUN


InputChecks_Z <- function(Z){
  
  if(!is.numeric(Z)) stop("Z is not numeric", call. = FALSE)
  if(!(is.matrix(Z) | is.data.frame(Z))) stop("Z must be either a data frame or matrix", call. = FALSE)
  
} # FUN


InputChecks_Z_CLAN <- function(Z_CLAN){
  
  if(!is.null(Z_CLAN)){
    
    if(!(is.matrix(Z_CLAN) | is.data.frame(Z_CLAN) | is.vector(Z_CLAN))){
      stop("Y must be either a vector, data frame, or matrix", call. = FALSE)
    }
  }
  
} # FUN




# helper that throws error in case of illegal input in 'X1.variables'
input.checks.X1 <- function(X1.variables, num.obs){
  
  if(!all(c("functions_of_Z", "custom_covariates", "fixed_effects") %in% names(X1.variables))){
      stop(paste0("The list ", deparse(substitute(X1.variables)), "must consist of three elements called ",
                  "'functions_of_Z', 'custom_covariates', and 'fixed_effects"), call. = FALSE)
    } # IF
  
  
  if(!(is.matrix(X1.variables$custom_covariates) | is.vector(X1.variables$custom_covariates) | is.null(X1.variables$custom_covariates))){
    stop(paste0(deparse(substitute(X1.variables)), "$custom_covariates must be either NULL, a vector, or a matrix"), call. = FALSE)
  } # IF
  
  
  if(!is.null(X1.variables$custom_covariates)){
    if(!is.numeric(X1.variables$custom_covariates)) stop(paste0(deparse(substitute(X1.variables)), "$custom_covariates must be numeric or NULL"), call. = FALSE)
    
    if(is.vector(X1.variables$custom_covariates)){
      
      if(length(X1.variables$custom_covariates) != num.obs){
        stop(paste0(deparse(substitute(X1.variables)), "$custom_covariates must be NULL or of the same length as Y"), call. = FALSE)
      } # IF
      
    } else{
      
      if(nrow(X1.variables$custom_covariates) != num.obs){
        stop(paste0(deparse(substitute(X1.variables)), "$custom_covariates must be NULL or of the same length as Y"), call. = FALSE)
      } # IF
      
    } # IF
    
  } # IF
  
  
  if(!(is.vector(X1.variables$fixed_effects) | is.null(X1.variables$fixed_effects))){
    stop(paste0(deparse(substitute(X1.variables)), "$fixed_effects must be either NULL or a vector"), call. = FALSE)
  } # IF
  
  
  if(is.vector(X1.variables$fixed_effects)){
    
    if(!is.numeric(X1.variables$fixed_effects)){
      stop(paste0(deparse(substitute(X1.variables)), "$fixed_effects must be numeric"), call. = FALSE)
    } else{
      
      if(length(X1.variables$fixed_effects) != num.obs){
        stop(paste0(deparse(substitute(X1.variables)), "$fixed_effects must be NULL or of the same length as Y"),
             call. = FALSE)
      } # IF
      
    }# IF
  } # IF
    
  
  legalinput <- X1.variables$functions_of_Z %in% c("S", "B", "p")
  
  if(!all(legalinput)){
    
    stop(paste0("Entries '", 
                paste(X1.variables$functions_of_Z[!legalinput], collapse = "', '"), 
                "' of ", deparse(substitute(X1.variables)),
                "$functions_of_Z are not contained in c('S', 'B', 'p')!"), call. = FALSE)
    
  } # IF
  
  if(!is.null(X1.variables$fixed_effects) & !is.vector(X1.variables$fixed_effects)){
    stop("The fixed effects need to be a vector", call. = FALSE)
  } # IF
  
} # FUN


