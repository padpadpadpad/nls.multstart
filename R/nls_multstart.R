#' Finds the best fit of non-linear model based on AIC score
#'
#' Finds the best estimated model using non-linear least squares regression using nlsLM(). The best fit is determined using AIC scores.
#'
#' @param formula a non-linear model formula, with the response on the left of a ~ operator and an expression involving parameters on the right
#' @param data data.frame (optional) in which to evaluate the variables in \code{formula} and \code{weights}
#' @param iter number of combinations of starting parameters that
#' are tried on each curve.
#' @param param_bds lower and upper boundaries for the start parameters. If
#' missing these default to +/- 1e+09. Need to specified as a vector as :
#' c(lower bound param 1, upper bound param 1, lower bound param 2, upper bound
#' param 2 etc)
#' @param supp_errors if \code{supp_errors = 'Y'}, then no error messages will be shown
#' from the nlsLM function, reducing the number of error messages
#' printed while the model attempts to converge using poor starting parameters. Advised to only use \code{supp_errors = 'Y'} when you are confident in the bounds of your starting parameters.
#' @param AICc whether or not the small sample AIC should be used. Defaults to
#' \code{'Y'}. Override this using \code{AICc == 'N'}. AICc should be used instead of AIC
#' when sample size is small in comparison to the number of estimated
#' parameters (Burnham & Anderson 2002 recommend its use when n / n_param < 40).
#' @param control specific control can be specified using \code{\link[minpack.lm]{nls.lm.control}}.
#' @param \dots Extra arguments to pass to \code{\link[minpack.lm]{nlsLM}} if necessary.
#' @return returns a nls object of the best estimated model fit
#' @note Useful additional arguments for \code{\link[minpack.lm]{nlsLM}} include: \code{na.action = na.omit},
#' \code{lower/upper = c()} where these represent upper and lower boundaries for parameter estimates
#' @author Daniel Padfield
#' @seealso
#' \code{\link[minpack.lm]{nlsLM}} for details on additional arguments to pass to the nlsLM function.
#'
#' \code{\link[MuMIn]{AICc}} for application of AICc.
#' @examples
#' # load in data
#'
#' data("Chlorella_TRC")
#' Chlorella_TRC_test <- Chlorella_TRC[Chlorella_TRC$curve_id == 1,]
#'
#' # run nlsLoop()
#'
#'# define the Sharpe-Schoolfield equation
#' schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
#'  Tc <- 273.15 + Tc
#'  k <- 8.62e-5
#'  boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp)))
#'  inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
#'  return(boltzmann.term + inactivation.term)
#'}
#'
#' fits <- nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
#'                 data = Chlorella_TRC_test,
#'                 iter = 500,
#'                 param_bds = c(-10, 10, 0.1, 2, 0.5, 5, 285, 330),
#'                 lower = c(lnc=-10, E=0, Eh=0, Th=0),
#'                 supp_errors = 'Y')
#'
#' @export

nls_multstart <-
  # arguments needed for nlsLoop ####
function(formula, data = parent.frame(), iter, param_bds, supp_errors = c('Y', 'N'), AICc = c('Y', 'N'), control, ...){

  # set default values
  if(missing(supp_errors)){supp_errors <- 'N'}
  if(missing(AICc)){AICc <- 'Y'}

  # checking whether MuMIn is installed
  if (!requireNamespace("MuMIn", quietly = TRUE)){
    stop("The MuMIn package is needed for calculation of AICc. Please install. Can be bypassed by using classic AIC using AICc = 'N'",
         call. = FALSE)
  }

  # create model ####
  formula <- stats::as.formula(formula)

  # define parameters to estimate and independent variable ####
  if(is.data.frame(data)) {
    params_ind <- all.vars(formula[[3]])[all.vars(formula[[3]]) %in% colnames(data)]
    params_est <- all.vars(formula[[3]])[! all.vars(formula[[3]]) %in% colnames(data)]
  } else if(is.environment(data)) {
    params_ind <- all.vars(formula[[3]])[all.vars(formula[[3]]) %in% names(data)]
    params_est <- all.vars(formula[[3]])[! all.vars(formula[[3]]) %in% names(data)]
  } else {
    stop("data should be a data.frame or an environment")
  }

  # set up parameter boundaries ####
  params_bds <- data.frame(param = params_est, stringsAsFactors = FALSE)
  params_bds$low.bds <- NA
  params_bds$high.bds <- NA

  if(missing(param_bds)){
    cat('No boundaries specified for the sought parameters. \n',
        'Default values of +/- 1e+10 will be used. This is likely to slow the process of finding the best model. \n')
    r <- readline("Continue with default values [y/n]? ")
    if(tolower(r) == "y") {
      params_bds$low.bds <- 10^-10
      params_bds$high.bds <- 10^10
    }
    if(tolower(r) == "n"){
      stop('Please enter upper and lower parameter boundaries as param_bds in function argument.')
    }

  }

  else {

    for(i in 1:nrow(params_bds)){
      params_bds$low.bds[i] <- param_bds[(2*i)-1]
      params_bds$high.bds[i] <- param_bds[2*i]
    }

  }

  # nlsLM controls - this can stay the same, potential to be overridden
  if(missing(control)) {
    control <- minpack.lm::nls.lm.control(maxiter = 1000, ftol = .Machine$double.eps, ptol = .Machine$double.eps)
  }


  # set up start values ####
  make_strt_values <- function(x, iter){
    strt_values <- data.frame(param = rep(x, times = iter),
                              value = stats::runif(iter, min = params_bds$low.bds[params_bds$param == x], max = params_bds$high.bds[params_bds$param == x]), stringsAsFactors = FALSE)
    return(strt_values)
  }

  strt <- plyr::ldply(as.list(params_est), make_strt_values, iter)

  # fit nls model using LM optimisation and using shotgun approach to get starting values ####

  # set count to 0
  count <- 0

  # create AIC vector
  stored_AIC <- 0

  # create empty fit
  fit <- NULL
  fit_best <- NULL

  for (j in 1:iter){
    # if((j/10) %% 1 == 0){cat(j, ' ')}
    # create start list
    start.vals <- list()
    for(k in 1:length(params_est)){
      start.vals[[params_est[k]]] <- strt[strt$param == params_est[k],]$value[j]
    }
    # try and fit the model for every set of searching parameters
    if(supp_errors == 'Y'){
      try(fit <- minpack.lm::nlsLM(formula,
                                   start=start.vals,
                                   control = control,
                                   data = data, ...),
          silent = TRUE)}
    if(supp_errors != 'Y'){
      try(fit <- minpack.lm::nlsLM(formula,
                                   start=start.vals,
                                   control = control,
                                   data = data, ...))}

    # if the AIC score of the next fit model is < the AIC of the stored fit fit, replace the fit
    # the output to ensure the best model is selected
    if(AICc == 'N'){

      if(is.null(fit) && is.null(fit_best)){
        count <-  0
      }
      else{
        count <- ifelse(stored_AIC <= stats::AIC(fit), count + 1, 0)
      }
      if(count == 100) break

      if(!is.null(fit) && stored_AIC == 0 | !is.null(fit) && stored_AIC > stats::AIC(fit)){
        stored_AIC <- stats::AIC(fit)
        fit_best <- fit
      }

    }
    else{

      if(is.null(fit) && is.null(fit_best)){
        count <-  0
      }
      else{
        count <- ifelse(stored_AIC <= MuMIn::AICc(fit), count + 1, 0)
      }

      if(count == 100) break


      if(!is.null(fit) && stored_AIC == 0 | !is.null(fit) && stored_AIC > MuMIn::AICc(fit)){

        stored_AIC <- MuMIn::AICc(fit)
        fit_best <- fit
      }

    }
  }

  return(fit_best)

}
