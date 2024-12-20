% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initializers.R
\name{setup_vcov}
\alias{setup_vcov}
\title{Setup function for \code{vcov_control} arguments}
\usage{
setup_vcov(estimator = "vcovHC", arguments = list(type = "const"))
}
\arguments{
\item{estimator}{Character specifying a covariance matrix estimator in the \href{https://CRAN.R-project.org/package=sandwich}{"sandwich"} package. Default is \code{"vcovHC"}. Supported estimators are \code{"vcovBS"}, \code{"vcovCL"}, \code{"vcovHAC"}, and \code{"vcovHC"}.}

\item{arguments}{A list of arguments that are to be passed to the function in the \code{"sandwich"} package that is specified in \code{estimator}. Default is \code{list(type = "const")}, which specifies the homoskedastic ordinary least squares covariance matrix estimator.}
}
\value{
An object of class \code{"setup_vcov"}, consisting of the following components:
\describe{
\item{\code{estimator}}{A character equal to covariance estimation function names in the \href{https://CRAN.R-project.org/package=sandwich}{"sandwich"} package.}
\item{\code{arguments}}{A list of arguments that shall be passed to the function specified in the \code{estimator} argument.}
}
See the description above for details.
}
\description{
Returns a list with two elements called \code{estimator} and \code{arguments}. The element \code{estimator} is a string specifying the covariance matrix estimator to be used in the linear regression regression of interest and needs to be a covariance estimator function in the \href{https://CRAN.R-project.org/package=sandwich}{"sandwich"} package. The second element, \code{arguments}, is a list of arguments that shall be passed to the function specified in the first element, \code{estimator}.
}
\details{
The output of this setup function is intended to be used as argument in the functions \code{\link{GenericML}()} and \code{\link{GenericML_single}()} (arguments \code{vcov_BLP}, \code{vcov_GATES}), as well as \code{\link{BLP}()} and \code{\link{GATES}()} (argument \code{vcov_control}).
}
\examples{
# use standard homoskedastic OLS covariance matrix estimate
setup_vcov(estimator = "vcovHC", arguments = list(type = "const"))

# use White's heteroskedasticity-robust estimator
setup_vcov(estimator = "vcovHC", arguments = list(type = "HC0"))

if (require("sandwich")){

# use HAC-robust estimator with prewhitening and Andrews' (Econometrica, 1991) weights
# since weightsAndrews() is a function in 'sandwich', require this package
setup_vcov(estimator = "vcovHAC", arguments = list(prewhite = TRUE, weights = weightsAndrews))

}

}
\references{
Zeileis A. (2004). \dQuote{Econometric Computing with HC and HAC Covariance Matrix Estimators.} \emph{Journal of Statistical Software}, \bold{11}(10), 1--17. \doi{10.18637/jss.v011.i10}

Zeileis A. (2006). \dQuote{Object-Oriented Computation of Sandwich Estimators.} \emph{Journal of Statistical Software}, \bold{16}(9), 1--16. \doi{10.18637/jss.v016.i09}
}
\seealso{
\code{\link{GenericML}()},
\code{\link{GenericML_single}()},
\code{\link{BLP}()},
\code{\link{GATES}()},
\code{\link{setup_X1}()},
\code{\link{setup_diff}()}
}
