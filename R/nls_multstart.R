#' Finds the best fit of non-linear model based on AIC score
#'
#' Finds the best estimated model using non-linear least squares regression using
#' nlsLM(). The best fit is determined using AIC scores.
#'
#' @param formula a non-linear model formula, with the response on the left of a
#'  ~ operator and an expression involving parameters on the right.
#' @param data (optional) data.frame, list or environment in which to evaluate
#'  the variables in \code{formula} and \code{modelweights}.
#' @param iter number of combinations of starting parameters which will be tried
#'  . If a single value is provided, then a shotgun/random-search approach will
#'  be used to sample starting parameters from a uniform distribution within the
#'  starting parameter bounds. If a vector of the same length as the number of
#'  parameters is provided, then a gridstart approach will be used to define
#'  each combination of that number of equally spaced intervals across each of
#'  the starting parameter bounds respectively. Thus, c(5,5,5) for three fitted
#'  parameters yields 125 model fits.  Supplying a vector for \code{iter}
#'  will override \code{convergence_count}.
#' @param start_lower lower boundaries for the start parameters. If missing, this
#'  will default to -1e+10.
#' @param start_upper upper boundaries for the start parameters. If missing, this
#'  will default to 1e+10.
#' @param supp_errors if \code{supp_errors = 'Y'}, then warning messages will be suppressed and no error messages from
#' \code{\link[minpack.lm]{nlsLM}} will be shown, reducing the number of error messages printed while the model attempts to converge using poor starting parameters.
#' We advise to only use \code{supp_errors = 'Y'} when confident in the bounds of
#'  your starting parameters.
#' @param convergence_count The number of counts that the winning model should be
#'  undefeated for before it is declared the winner. This argument defaults to
#'  100. If specified as \code{FALSE}, then all of the iterations will be fitted, and
#'  the best model selected. Note that \code{convergence_count} can only be used
#'  with a shotgun/random-search approach, and not with a gridstart approach.
#'  This argument will be ignored if a gridstart approach is specified by a
#'  vector input for \code{iter}.
#' @param control specific control can be specified using
#'  \code{\link[minpack.lm]{nls.lm.control}}.
#' @param modelweights Optional model weights for the nls. If \code{data} is specified, then this argument should be the name of the numeric weights vector within
#'  the \code{data} object.
#' @param \dots Extra arguments to pass to \code{\link[minpack.lm]{nlsLM}} if
#'  necessary.
#' @return returns a nls object of the best estimated model fit.
#' @note Useful additional arguments for \code{\link[minpack.lm]{nlsLM}} include:
#'  \code{na.action = na.omit}, \code{lower/upper = c()} where these represent
#'  upper and lower boundaries for parameter estimates.
#' @author Daniel Padfield
#' @author Granville Matheson
#' @seealso \code{\link[minpack.lm]{nlsLM}} for details on additional arguments
#' to pass to the nlsLM function.
#' @examples
#' # load in data
#'
#' data("Chlorella_TRC")
#' Chlorella_TRC_test <- Chlorella_TRC[Chlorella_TRC$curve_id == 1,]
#'
#' # run nls_multstart()
#'
#' # define the Sharpe-Schoolfield equation
#' schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
#'  Tc <- 273.15 + Tc
#'  k <- 8.62e-5
#'  boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp)))
#'  inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
#'  return(boltzmann.term + inactivation.term)
#' }
#'
#' fits <- nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
#'                 data = Chlorella_TRC_test,
#'                 iter = 500,
#'                 start_lower = c(lnc=-10, E=0.1, Eh=0.5, Th=285),
#'                 start_upper = c(lnc=10, E=2, Eh=5, Th=330),
#'                 lower = c(lnc=-10, E=0, Eh=0, Th=0),
#'                 supp_errors = 'Y')
#'
#' @name nls_multstart
#' @export

# getting rid of undefined variables note
if(base::getRversion() >= "2.15.1")  utils::globalVariables(c("AICval", 'iteration', 'startpars'))

nls_multstart <-
  # arguments needed for nls_multstart ####
  function(formula, data = parent.frame(), iter, start_lower, start_upper,
           supp_errors = c("Y", "N"), convergence_count = 100, control,
           modelweights, ...) {

    # set default values
    if (missing(supp_errors)) {
      supp_errors <- "N"
    }

    # create model ####
    formula <- stats::as.formula(formula)

    # define parameters to estimate and independent variable ####
    if (any(class(data) %in% c("data.frame", "list", "environment"))) {
      params_ind <- all.vars(formula[[3]])[all.vars(formula[[3]]) %in% names(data)]
      params_est <- all.vars(formula[[3]])[!all.vars(formula[[3]]) %in% names(data)]
      params_dep <- all.vars(formula[[2]])
    } else {
      stop("data should be a data.frame, list or an environment")
    }

    # set up parameter boundaries ####

    if (missing(start_lower) || missing(start_upper)) {
      cat(
        "No boundaries specified for the starting values of sought parameters. \n",
        "Default values of +/- 1e+10 will be used. This is likely \n",
        "to slow the process of finding the best model. \n"
      )
      r <- readline("Continue with default values [y/n]? ")

      if (tolower(r) == "n") {
        stop("Please enter upper and lower parameter boundaries as start_lower and start_upper in function argument.")
      }
    }

    if (missing(start_lower)) start_lower <- rep(-10 ^ 10, length(params_est))
    if (missing(start_upper)) start_upper <- rep(10 ^ 10, length(params_est))

    if (length(start_lower) != length(params_est) || length(start_upper) != length(params_est)) {
      stop("There must be as many parameter starting bounds as there are parameters")
    }

    params_bds <- data.frame(
      param = params_est,
      low.bds = unlist(start_lower),
      high.bds = unlist(start_upper),
      stringsAsFactors = FALSE
    )

    # nlsLM controls - this can stay the same, potential to be overridden
    if (missing(control)) {
      control <- minpack.lm::nls.lm.control(maxiter = 1000, ftol = .Machine$double.eps, ptol = .Machine$double.eps)
    }


    # transform input arguments
    silent <- ifelse(supp_errors == "Y", TRUE, FALSE)

    # if silent is TRUE, temporarily switch off warnings
    if(silent == TRUE){
      oo <- options(warn=-1)
      on.exit(options(oo))
    }

    if ("modelweights" %in% all.vars(formula)) {
      stop(paste0(
        "The variable name 'modelweights' is reserved for model weights. Please change the name\n",
        "of this variable"
      ))
    }

    if (missing(modelweights)) {
      data$modelweights <- rep(1, length(data[[params_dep]]))
    } else {
      data$modelweights <- eval(substitute(modelweights), data)
    }


    if (length(iter) == 1) {
      multistart_type <- "shotgun"
    } else if (length(iter) == nrow(params_bds)) {
      multistart_type <- "gridstart"
    } else {
      stop(paste0(
        "iter should be of length 1 for shotgun approach and of the same length as the\n",
        "number of parameters for the gridstart approach."
      ))
    }

    if (multistart_type == "gridstart" && convergence_count != FALSE) {
      warning("A gridstart approach cannot be applied with convergence_count. Convergence count will be set to FALSE")
      convergence_count <- FALSE
    }


    # set up start values ####

    ## SHOTGUN ##

    if (multistart_type == "shotgun") {
      strt <- purrr::map2(params_bds$low.bds, params_bds$high.bds, ~runif(iter, .x, .y))
      names(strt) <- params_bds$param
      strt <- dplyr::bind_rows(strt)
    }


    ## GRIDSTART ##

    if (multistart_type == "gridstart") {
      params_bds$iter <- iter
      strt <- purrr::pmap(
        as.list(params_bds[, -1]),
        function(low.bds, high.bds, iter) seq(from = low.bds, to = high.bds, length.out = iter)
      )
      names(strt) <- params_bds$param
      strt <- tibble::as_tibble(expand.grid(strt))
    }


    # Fit nls model using LM optimisation across multiple starting values


    #### With convergence ####

    if (multistart_type == "shotgun" && convergence_count != FALSE) {

      # set count to 0
      count <- 0

      # create AIC vector
      stored_AIC <- Inf

      # create empty fit
      fit <- NULL
      fit_best <- NULL

      for (j in 1:nrow(strt)) {
        start.vals <- as.list(strt[j, ])

        # try and fit the model for every set of searching parameters

        try(
          fit <- minpack.lm::nlsLM(
            formula,
            start = start.vals,
            control = control,
            data = data,
            weights = modelweights, ...
          ),
          silent = silent
        )


        # if the AIC score of the next fit model is < the AIC of the stored fit fit, replace the fit
        # the output to ensure the best model is selected

        if (is.null(fit) && is.null(fit_best)) {
          count <- 0
        }
        else {
          count <- ifelse(stored_AIC <= stats::AIC(fit), count + 1, 0)
        }
        if (count == convergence_count) break

        if (!is.null(fit) && stored_AIC > stats::AIC(fit)) {
          stored_AIC <- stats::AIC(fit)
          fit_best <- fit
        }
      }
    }


    #### Without convergence_count ####

    if (convergence_count == FALSE) {
      strt$iteration <- 1:nrow(strt)

      allfits <- tidyr::nest(strt, startpars = -iteration)


      # create empty fit
      fit <- NULL
      fit_best <- NULL


      # fit_function saves only AIC values so as not to fill memory with model fits

      fit_aic <- function(startpars) {
        start.vals <- as.list(startpars[[1]])

        try(
          fit <- minpack.lm::nlsLM(
            formula,
            start = start.vals,
            control = control,
            data = data,
            weights = modelweights, ...
          ),
          silent = silent
        )

        AICval <- ifelse(!is.null(fit), stats::AIC(fit), Inf)

        return(AICval)
      }

      allfits <- dplyr::group_by(allfits, iteration)
      allfits <- dplyr::mutate(allfits, AICval = purrr::map_dbl(
        startpars,
        ~fit_aic(startpars)
      ))
      allfits <- dplyr::arrange(allfits, AICval)

      fit_best <- minpack.lm::nlsLM(
        formula,
        start = allfits$startpars[[1]],
        control = control,
        data = data,
        weights = modelweights, ...
      )
    }

    return(fit_best)
  }
