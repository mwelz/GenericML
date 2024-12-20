% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accessors.R
\name{get_best}
\alias{get_best}
\title{Accessor function for the best learner estimates}
\usage{
get_best(x)
}
\arguments{
\item{x}{An object of the class \code{"\link{GenericML}"}, as returned by the function \code{\link{GenericML}()}.}
}
\value{
An object of class \code{"best"}, which consists of the following components:
\describe{
\item{\code{BLP}}{A string holding the name of the best learner for a BLP analysis.}
\item{\code{GATES}}{A string holding the name of the best learner for a GATES analysis.}
\item{\code{CLAN}}{A string holding the name of the best learner for a CLAN analysis (same learner as in \code{GATES}).}
\item{\code{overview}}{A numeric matrix of the estimates of the performance measures \eqn{\Lambda} and \eqn{\bar{\Lambda}} for each learner.}}
}
\description{
The best learner is determined by maximizing the criteria \eqn{\Lambda} and \eqn{\bar{\Lambda}}, see Sections 5.2 and 5.3 of the paper. This function accesses the estimates of these two criteria,
}
\examples{
if(require("rpart") && require("ranger")){
## generate data
set.seed(1)
n  <- 150                                  # number of observations
p  <- 5                                    # number of covariates
D  <- rbinom(n, 1, 0.5)                    # random treatment assignment
Z  <- matrix(runif(n*p), n, p)             # design matrix
Y0 <- as.numeric(Z \%*\% rexp(p) + rnorm(n)) # potential outcome without treatment
Y1 <- 2 + Y0                               # potential outcome under treatment
Y  <- ifelse(D == 1, Y1, Y0)               # observed outcome

## column names of Z
colnames(Z) <- paste0("V", 1:p)

## specify learners
learners <- c("tree", "mlr3::lrn('ranger', num.trees = 10)")

## perform generic ML inference
# small number of splits to keep computation time low
x <- GenericML(Z, D, Y, learners, num_splits = 2,
               parallel = FALSE)

## access best learner
get_best(x)

## access BLP generic targets for best learner w/o plot
get_BLP(x, learner = "best", plot = FALSE)

## access BLP generic targets for ranger learner w/o plot
get_BLP(x, learner = "mlr3::lrn('ranger', num.trees = 10)", plot = FALSE)

## access GATES generic targets for best learner w/o plot
get_GATES(x, learner = "best", plot = FALSE)

## access GATES generic targets for ranger learner w/o plot
get_GATES(x, learner = "mlr3::lrn('ranger', num.trees = 10)", plot = FALSE)

## access CLAN generic targets for "V1" & best learner, w/o plot
get_CLAN(x, learner = "best", variable = "V1", plot = FALSE)

## access CLAN generic targets for "V1" & ranger learner, w/o plot
get_CLAN(x, learner = "mlr3::lrn('ranger', num.trees = 10)",
         variable = "V1", plot = FALSE)
}

}
\seealso{
\code{\link{GenericML}()},
\code{\link{get_BLP}()},
\code{\link{get_GATES}()},
\code{\link{get_CLAN}()}
}
