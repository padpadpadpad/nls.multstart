context("test-nls_multstart.R")

data("Chlorella_TRC")
Chlorella_TRC_test <- Chlorella_TRC[Chlorella_TRC$curve_id == 1,]
Chlorella_TRC_test$mweights <- c(rep(0, 6), rep(1, 6))

# run nls_multstart()

# define the Sharpe-Schoolfield equation
schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
  Tc <- 273.15 + Tc
  k <- 8.62e-5
  boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp)))
  inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
  return(boltzmann.term + inactivation.term)
}



test_that("fit with convergence works", {
  expect_equal(
    round( coef(nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                              data = Chlorella_TRC_test,
                              iter = 500,
                              start_lower = c(-10, 0.1, 0.5, 285),
                              start_upper = c(10, 2, 5, 330),
                              lower = c(lnc=-10, E=0, Eh=0, Th=0),
                              supp_errors = 'Y') )),
    c(lnc=-1, E=1, Eh=4, Th=312) )
})

test_that("fit without convergence works", {
  expect_equal(
    round( coef(nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                              data = Chlorella_TRC_test,
                              iter = 250,
                              start_lower = c(-10, 0.1, 0.5, 285),
                              start_upper = c(10, 2, 5, 330),
                              lower = c(lnc=-10, E=0, Eh=0, Th=0),
                              supp_errors = 'Y', convergence_count = F) )),
    c(lnc=-1, E=1, Eh=4, Th=312) )
})


test_that("weighted fit with convergence works", {
  expect_equal(
    round( coef(nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                              data = Chlorella_TRC_test,
                              iter = 500,
                              start_lower = c(-10, 0.1, 0.5, 285),
                              start_upper = c(10, 2, 5, 330),
                              lower = c(lnc=-10, E=0, Eh=0, Th=0),
                              supp_errors = 'Y',
                              modelweights = mweights) )),
    c(lnc=-5, E=3, Eh=6, Th=311) )
})

test_that("weighted fit without convergence works", {
  expect_equal(
    round( coef(nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                              data = Chlorella_TRC_test,
                              iter = 250,
                              start_lower = c(-10, 0.1, 0.5, 285),
                              start_upper = c(10, 2, 5, 330),
                              lower = c(lnc=-10, E=0, Eh=0, Th=0),
                              supp_errors = 'Y', convergence_count = F,
                              modelweights = mweights) )),
    c(lnc=-5, E=3, Eh=6, Th=311) )
})
