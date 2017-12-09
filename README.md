<!-- README.md is generated from README.Rmd. Please edit that file -->
nls.multstart
-------------

Robust and reproducible non-linear regression in R

### Issues and suggestions

Please report any issues/suggestions for improvement in the [issues link](https://github.com/padpadpadpad/nls.multstart/issues) for the repository. Or please email <d.padfield@exeter.ac.uk>.

[![Travis-CI Build Status](https://travis-ci.org/padpadpadpad/nls.multstart.svg?branch=master)](https://travis-ci.org/padpadpadpad/nls.multstart)

### Licensing

This package is licensed under GPL-3.

### Overview

**nls.multstart** is an R package that allows more robust and reproducible non-linear regression compared to **nls** or **nlsLM**. These functions allow only a single starting value, meaning that it can be hard to get the best estimated model. This is especially true if the same model is fitted over the levels of a factor, which may have the same shape of curve, but be much different in terms of parameter estimates.

**nls\_multstart()** is the main (currently only) function of **nls.multstart**. Similar to the R package **nls2**, it allows multiple starting values for each parameter and then iterates through multiple starting values, attempting a fit with each set of start parameters. The best model is then picked on AIC score. This results in a more reproducible and reliable method of fitting non-linear least squares regression in R.

This package is designed to work with the **tidyverse**, harnessing the functions within **broom**, **tidyr**, **dplyr** and **purrr** to extract estimates and plot things easily with **ggplot2**. A slightly less tidy-friendly implementation is [**nlsLoop**](https://github.com/padpadpadpad/nlsLoop).

### Installation and examples

#### 1. Installation

R packages in GitHub can be installed using **devtools**.

``` r
# install package
devtools::install_github("padpadpadpad/nls.multstart")
```

#### 2. Run nls\_multstart()

**nls\_multstart** can be used to do non-linear regression on a single curve

``` r

# load in nlsLoop and other packages
library(nls.multstart)
library(ggplot2)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)

# load in example data set
data("Chlorella_TRC")

# define the Sharpe-Schoolfield equation
schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
  Tc <- 273.15 + Tc
  k <- 8.62e-5
  boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp)))
  inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
  return(boltzmann.term + inactivation.term)
  }
```

``` r
 # subset dataset
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

fit
#> Nonlinear regression model
#>   model: ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20)
#>    data: data
#>      lnc        E       Eh       Th 
#>  -1.3462   0.9877   4.3326 312.1887 
#>  residual sum-of-squares: 7.257
#> 
#> Number of iterations to convergence: 13 
#> Achieved convergence tolerance: 1.49e-08
```

#### 3. Clean up fit

This fit can then be tidied up in various ways using the R package **broom**. Each different function in **broom** returns a different set of information. **tidy** returns the estimated parameters, **augment** returns the predictions and **glance** returns information about the model such as AIC score.

``` r
# get info
info <- glance(fit)

# get params
params <- tidy(fit)

# get predictions
preds <- augment(fit)

info
#>       sigma isConv       finTol    logLik      AIC      BIC deviance
#> 1 0.9524198   TRUE 1.490116e-08 -14.00948 38.01896 40.44349 7.256827
#>   df.residual
#> 1           8
params
#>   term    estimate std.error statistic      p.value
#> 1  lnc  -1.3462105 0.4656398 -2.891098 2.016515e-02
#> 2    E   0.9877307 0.4521481  2.184529 6.043425e-02
#> 3   Eh   4.3326452 1.4877826  2.912149 1.952469e-02
#> 4   Th 312.1887459 3.8781636 80.499117 6.323537e-13
preds
#>        ln.rate      K     .fitted      .resid
#> 1  -2.06257833 289.15 -1.88694034 -0.17563799
#> 2  -1.32437939 292.15 -1.48002016  0.15564077
#> 3  -0.95416807 295.15 -1.08143501  0.12726694
#> 4  -0.79443675 298.15 -0.69121465 -0.10322210
#> 5  -0.18203642 301.15 -0.31058072  0.12854430
#> 6   0.17424007 304.15  0.05336433  0.12087574
#> 7  -0.04462754 307.15  0.36657463 -0.41120217
#> 8   0.48050690 310.15  0.49837148 -0.01786458
#> 9   0.38794188 313.15  0.17973799  0.20820389
#> 10  0.39365516 316.15 -0.64473314  1.03838830
#> 11 -3.86319577 319.15 -1.70300699 -2.16018878
#> 12 -1.72352435 322.15 -2.81272004  1.08919569
```

#### 4. Plot fit

The predictions can then easily be plotted alongside the actual data.

``` r

ggplot() +
  geom_point(aes(K, ln.rate), d_1) +
  geom_line(aes(K, .fitted), preds)
```

![](README-plot_one_fit-1.png)

#### 5. Fitting over levels of a factor with nls\_multstart

**nls\_multstart** is unlikely to speed you up very much if only one curve is fitted. However, if you have 10, 60 or 100s of curves to fit, it makes sense that at least some of them may not fit with the same starting parameters, no matter how many iterations it is run for.

This is where **nls\_multstart** can help. Multiple models can be fitted using **purrr**, **dplyr** and **tidyr**. These fits can then be tidied using **broom**, an approach Hadley Wickham has previously [written about](https://blog.rstudio.com/2016/02/02/tidyr-0-4-0/).

``` r
# fit over each set of groupings
fits <- Chlorella_TRC %>%
  group_by(., flux, growth.temp, process, curve_id) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                                   data = .x,
                                   iter = 1000,
                                   param_bds = c(-10, 10, 0.1, 2, 0.5, 10, 285, 330),
                                   supp_errors = 'Y',
                                   AICc = 'Y',
                                   na.action = na.omit,
                                   lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))))
```

A single fit can check to make sure it looks ok.

``` r
# look at output object
fits
#> # A tibble: 60 x 6
#>           flux growth.temp     process curve_id              data
#>          <chr>       <dbl>       <chr>    <dbl>            <list>
#>  1 respiration          20 acclimation        1 <tibble [12 x 3]>
#>  2 respiration          20 acclimation        2 <tibble [12 x 3]>
#>  3 respiration          23 acclimation        3 <tibble [12 x 3]>
#>  4 respiration          27 acclimation        4  <tibble [9 x 3]>
#>  5 respiration          27 acclimation        5 <tibble [12 x 3]>
#>  6 respiration          30 acclimation        6 <tibble [12 x 3]>
#>  7 respiration          30 acclimation        7 <tibble [12 x 3]>
#>  8 respiration          33 acclimation        8 <tibble [10 x 3]>
#>  9 respiration          33 acclimation        9  <tibble [8 x 3]>
#> 10 respiration          20 acclimation       10 <tibble [10 x 3]>
#> # ... with 50 more rows, and 1 more variables: fit <list>

# look at a single fit
summary(fits$fit[[1]])
#> 
#> Formula: ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20)
#> 
#> Parameters:
#>     Estimate Std. Error t value Pr(>|t|)    
#> lnc  -1.3462     0.4656  -2.891   0.0202 *  
#> E     0.9877     0.4521   2.185   0.0604 .  
#> Eh    4.3326     1.4878   2.912   0.0195 *  
#> Th  312.1887     3.8782  80.499 6.32e-13 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.9524 on 8 degrees of freedom
#> 
#> Number of iterations to convergence: 13 
#> Achieved convergence tolerance: 1.49e-08
```

#### 6. Clean up multiple fits

These fits can be cleaned up in a similar way to the single fit, but this time purrr::map() iterates the **broom** function over the grouping variables.

``` r
# get summary
info <- fits %>%
  unnest(fit %>% map(glance))

# get params
params <- fits %>%
  unnest(fit %>% map(tidy))

# get predictions
preds <- fits %>%
  unnest(fit %>% map(augment))
```

Looking at **info** allows us to see if all the models converged.

``` r
info
#> # A tibble: 60 x 14
#>           flux growth.temp     process curve_id              data
#>          <chr>       <dbl>       <chr>    <dbl>            <list>
#>  1 respiration          20 acclimation        1 <tibble [12 x 3]>
#>  2 respiration          20 acclimation        2 <tibble [12 x 3]>
#>  3 respiration          23 acclimation        3 <tibble [12 x 3]>
#>  4 respiration          27 acclimation        4  <tibble [9 x 3]>
#>  5 respiration          27 acclimation        5 <tibble [12 x 3]>
#>  6 respiration          30 acclimation        6 <tibble [12 x 3]>
#>  7 respiration          30 acclimation        7 <tibble [12 x 3]>
#>  8 respiration          33 acclimation        8 <tibble [10 x 3]>
#>  9 respiration          33 acclimation        9  <tibble [8 x 3]>
#> 10 respiration          20 acclimation       10 <tibble [10 x 3]>
#> # ... with 50 more rows, and 9 more variables: fit <list>, sigma <dbl>,
#> #   isConv <lgl>, finTol <dbl>, logLik <dbl>, AIC <dbl>, BIC <dbl>,
#> #   deviance <dbl>, df.residual <int>
```

#### 7. Plotting predictions

When plotting non-linear fits, it often looks better to have a smooth curve, even if there are not many points underlying the fit. This can be achieved by including `newdata` in the **augment** function and creating a higher resolution set of predictor values.

However, when predicting for many different fits, it is not certain that each curve has the same range of predictor variables. Consequently, we need to filter each new prediction by the **min** and **max** of the predictor variables.

``` r
# new data frame of predictions
new_preds <- Chlorella_TRC %>%
  do(., data.frame(K = seq(min(.$K), max(.$K), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(Chlorella_TRC, curve_id) %>%
  summarise(., min_K = min(K), max_K = max(K)) %>%
  ungroup()

# create new predictions
preds2 <- fits %>%
  unnest(fit %>% map(augment, newdata = new_preds)) %>%
  merge(., max_min, by = 'curve_id') %>%
  group_by(., curve_id) %>%
  filter(., K > unique(min_K) & K < unique(max_K)) %>%
  rename(., ln.rate = .fitted) %>%
  ungroup()
```

These can then be plotted using **ggplot2**.

``` r
# plot
ggplot() +
  geom_point(aes(K - 273.15, ln.rate, col = flux), size = 2, Chlorella_TRC) +
  geom_line(aes(K - 273.15, ln.rate, col = flux, group = curve_id), alpha = 0.5, preds2) +
  facet_wrap(~ curve_id, labeller = labeller(.multi_line = F)) +
  scale_colour_manual(values = c('green4', 'black')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab('log Metabolic rate') +
  xlab('Assay temperature (ºC)') +
  theme(legend.position = c(0.9, 0.15))
```

![](README-plot_many_fits-1.png)
