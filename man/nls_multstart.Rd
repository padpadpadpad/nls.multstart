% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nls_multstart.R
\name{nls_multstart}
\alias{nls_multstart}
\title{Finds the best fit of non-linear model based on AIC score}
\arguments{
\item{formula}{a non-linear model formula, with the response on the left of a
~ operator and an expression involving parameters on the right.}

\item{data}{(optional) data.frame, list or environment in which to evaluate
the variables in \code{formula} and \code{modelweights}.}

\item{iter}{number of combinations of starting parameters which will be tried
. If a single value is provided, then a shotgun/random-search approach will
be used to sample starting parameters from a uniform distribution within
the starting parameter bounds. If a vector of the same length as the number
of parameters is provided, then a gridstart approach will be used to define
each combination of that number of equally spaced intervals across each of
the starting parameter bounds respectively. Thus, c(5,5,5) for three fitted
parameters yields 125 model fits.  Supplying a vector for \code{iter} will
override \code{convergence_count}.}

\item{start_lower}{lower boundaries for the start parameters. If missing,
this will default to -1e+10.}

\item{start_upper}{upper boundaries for the start parameters. If missing,
this will default to 1e+10.}

\item{supp_errors}{if \code{supp_errors = 'Y'}, then warning messages will be
suppressed and no error messages from \code{\link[minpack.lm]{nlsLM}} will
be shown, reducing the number of error messages printed while the model
attempts to converge using poor starting parameters. We advise to only use
\code{supp_errors = 'Y'} when confident in the bounds of your starting
parameters.}

\item{convergence_count}{The number of counts that the winning model should
be undefeated for before it is declared the winner. This argument defaults
to 100. If specified as \code{FALSE}, then all of the iterations will be
fitted, and the best model selected. Note that \code{convergence_count} can
only be used with a shotgun/random-search approach, and not with a
gridstart approach. This argument will be ignored if a gridstart approach
is specified by a vector input for \code{iter}.}

\item{control}{specific control can be specified using
\code{\link[minpack.lm]{nls.lm.control}}.}

\item{modelweights}{Optional model weights for the nls. If \code{data} is
specified, then this argument should be the name of the numeric weights
vector within the \code{data} object.}

\item{\dots}{Extra arguments to pass to \code{\link[minpack.lm]{nlsLM}} if
necessary.}
}
\value{
returns a nls object of the best estimated model fit.
}
\description{
Finds the best estimated model using non-linear least squares regression
using nlsLM(). The best fit is determined using AIC scores.
}
\note{
Useful additional arguments for \code{\link[minpack.lm]{nlsLM}}
  include: \code{na.action = na.omit}, \code{lower/upper = c()} where these
  represent upper and lower boundaries for parameter estimates.
}
\examples{
# load in data

data("Chlorella_TRC")
Chlorella_TRC_test <- Chlorella_TRC[Chlorella_TRC$curve_id == 1,]

# run nls_multstart()

# define the Sharpe-Schoolfield equation
schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
 Tc <- 273.15 + Tc
 k <- 8.62e-5
 boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp)))
 inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
 return(boltzmann.term + inactivation.term)
}

fits <- nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                data = Chlorella_TRC_test,
                iter = 500,
                start_lower = c(lnc=-10, E=0.1, Eh=0.5, Th=285),
                start_upper = c(lnc=10, E=2, Eh=5, Th=330),
                lower = c(lnc=-10, E=0, Eh=0, Th=0),
                supp_errors = 'Y')

}
\seealso{
\code{\link[minpack.lm]{nlsLM}} for details on additional arguments
  to pass to the nlsLM function.
}
\author{
Daniel Padfield

Granville Matheson
}
