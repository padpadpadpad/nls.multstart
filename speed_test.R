# nls_multstart test
library(nls.multstart)
library(dplyr)
library(broom)
library(tidyr)
library(purrr)
library(microbenchmark)

data("Chlorella_TRC")

schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
  Tc <- 273.15 + Tc
  k <- 8.62e-5
  boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp)))
  inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
  return(boltzmann.term + inactivation.term)
}

d_1 <- subset(Chlorella_TRC, curve_id == 1)

# run nls_multstart
fit <- nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                     data = d_1,
                     iter = 500,
                     param_bds = c(-10, 10, 0.1, 2, 0.5, 5, 285, 330),
                     supp_errors = 'Y',
                     AICc = 'Y',
                     na.action = na.omit,
                     lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))

# run speed test
res <- microbenchmark(fit1 = Chlorella_TRC %>%
  group_by(., flux, growth.temp, process, curve_id) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                                                data = .x,
                                                iter = 1000,
                                                param_bds = c(-10, 10, 0.1, 2, 0.5, 10, 285, 330),
                                                supp_errors = 'Y',
                                                AICc = 'Y',
                                                na.action = na.omit,
                                                lower = c(lnc = -10, E = 0, Eh = 0, Th = 0)))),
  fit2 = Chlorella_TRC %>%
  group_by(., flux, growth.temp, process, curve_id) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                                                data = .x,
                                                iter = 2000,
                                                param_bds = c(-10, 10, 0.1, 2, 0.5, 10, 285, 330),
                                                supp_errors = 'Y',
                                                AICc = 'Y',
                                                na.action = na.omit,
                                                lower = c(lnc = -10, E = 0, Eh = 0, Th = 0)))),
  times = 5
)
