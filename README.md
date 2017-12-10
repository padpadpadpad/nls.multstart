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

**nls.multstart** is an R package that allows more robust and reproducible non-linear regression compared to **nls()** or **nlsLM()**. These functions allow only a single starting value, meaning that it can be hard to get the best estimated model. This is especially true if the same model is fitted over the levels of a factor, which may have the same shape of curve, but be much different in terms of parameter estimates.

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

**nls\_multstart()** can be used to do non-linear regression on a single curve

``` r

# load in nlsLoop and other packages
library(nls.multstart)
library(ggplot2)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)

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
#> Number of iterations to convergence: 18 
#> Achieved convergence tolerance: 1.49e-08
```

#### 3. Clean up fit

This fit can then be tidied up in various ways using the R package **broom**. Each different function in **broom** returns a different set of information. **tidy()** returns the estimated parameters, **augment()** returns the predictions and **glance()** returns information about the model such as AIC score. Confidence intervals of non-linear regression can also be estimated using **nlstools::confint2()**

``` r
# get info
info <- glance(fit)
info
#>       sigma isConv       finTol    logLik      AIC      BIC deviance
#> 1 0.9524198   TRUE 1.490116e-08 -14.00948 38.01896 40.44349 7.256827
#>   df.residual
#> 1           8

# get params
params <- tidy(fit)

# get confidence intervals using nlstools
CI <- confint2(fit) %>%
  data.frame() %>%
  rename(., conf.low = X2.5.., conf.high = X97.5..)

# bind params and confidence intervals
params <- bind_cols(params, CI)
select(params, -c(statistic, p.value))
#>   term    estimate std.error     conf.low   conf.high
#> 1  lnc  -1.3462105 0.4656398  -2.41997789  -0.2724432
#> 2    E   0.9877307 0.4521481  -0.05492466   2.0303860
#> 3   Eh   4.3326452 1.4877826   0.90181231   7.7634782
#> 4   Th 312.1887459 3.8781636 303.24568449 321.1318073

# get predictions
preds <- augment(fit)
preds
#>        ln.rate      K     .fitted      .resid
#> 1  -2.06257833 289.15 -1.88694035 -0.17563798
#> 2  -1.32437939 292.15 -1.48002017  0.15564078
#> 3  -0.95416807 295.15 -1.08143502  0.12726694
#> 4  -0.79443675 298.15 -0.69121466 -0.10322209
#> 5  -0.18203642 301.15 -0.31058073  0.12854430
#> 6   0.17424007 304.15  0.05336433  0.12087574
#> 7  -0.04462754 307.15  0.36657462 -0.41120216
#> 8   0.48050690 310.15  0.49837148 -0.01786458
#> 9   0.38794188 313.15  0.17973800  0.20820389
#> 10  0.39365516 316.15 -0.64473314  1.03838829
#> 11 -3.86319577 319.15 -1.70300698 -2.16018879
#> 12 -1.72352435 322.15 -2.81272003  1.08919568
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

**nls\_multstart()** is unlikely to speed you up very much if only one curve is fitted. However, if you have 10, 60 or 100s of curves to fit, it makes sense that at least some of them may not fit with the same starting parameters, no matter how many iterations it is run for.

This is where **nls\_multstart()** can help. Multiple models can be fitted using **purrr**, **dplyr** and **tidyr**. These fits can then be tidied using **broom**, an approach Hadley Wickham has previously [written about](https://blog.rstudio.com/2016/02/02/tidyr-0-4-0/).

``` r
# fit over each set of groupings
fits <- Chlorella_TRC %>%
  group_by(., flux, growth.temp, process, curve_id) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                                   data = .x,
                                   iter = 1000,
                                   param_bds = c(-1000, 1000, 0.1, 2, 0.5, 10, 285, 330),
                                   supp_errors = 'Y',
                                   AICc = 'Y',
                                   na.action = na.omit,
                                   lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))))
```

A single fit can check to make sure it looks ok. Looking at `fits` demonstrates that there is now a `fit` list column containing each of the non-linear fits for each combination of our grouping variables.

``` r
# look at output object
select(fits, curve_id, data, fit)
#> # A tibble: 60 x 3
#>    curve_id              data       fit
#>       <dbl>            <list>    <list>
#>  1        1 <tibble [12 x 3]> <S3: nls>
#>  2        2 <tibble [12 x 3]> <S3: nls>
#>  3        3 <tibble [12 x 3]> <S3: nls>
#>  4        4  <tibble [9 x 3]> <S3: nls>
#>  5        5 <tibble [12 x 3]> <S3: nls>
#>  6        6 <tibble [12 x 3]> <S3: nls>
#>  7        7 <tibble [12 x 3]> <S3: nls>
#>  8        8 <tibble [10 x 3]> <S3: nls>
#>  9        9  <tibble [8 x 3]> <S3: nls>
#> 10       10 <tibble [10 x 3]> <S3: nls>
#> # ... with 50 more rows

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
#> Number of iterations to convergence: 14 
#> Achieved convergence tolerance: 1.49e-08
```

#### 6. Clean up multiple fits

These fits can be cleaned up in a similar way to the single fit, but this time **purrr::map()** iterates the **broom** function over the grouping variables.

``` r
# get summary
info <- fits %>%
  unnest(fit %>% map(glance))

# get params
params <- fits %>%
  unnest(fit %>% map(tidy))

# get confidence intervals
CI <- fits %>% 
  unnest(fit %>% map(~ confint2(.x) %>%
  data.frame() %>%
  rename(., conf.low = X2.5.., conf.high = X97.5..))) %>%
  group_by(., curve_id) %>%
  mutate(., term = c('lnc', 'E', 'Eh', 'Th')) %>%
  ungroup()

# merge parameters and CI estimates
params <- merge(params, CI, by = intersect(names(params), names(CI)))

# get predictions
preds <- fits %>%
  unnest(fit %>% map(augment))
```

Looking at **info** allows us to see if all the models converged.

``` r
select(info, curve_id, logLik, AIC, BIC, deviance, df.residual)
#> # A tibble: 60 x 6
#>    curve_id      logLik       AIC       BIC  deviance df.residual
#>       <dbl>       <dbl>     <dbl>     <dbl>     <dbl>       <int>
#>  1        1 -14.0094789 38.018958 40.443491 7.2568273           8
#>  2        2  -1.1969914 12.393983 14.818516 0.8577249           8
#>  3        3  -7.3855699 24.771140 27.195673 2.4059814           8
#>  4        4  -0.5234387 11.046877 12.033000 0.5919502           5
#>  5        5 -10.8498521 31.699704 34.124237 4.2859330           8
#>  6        6  -8.5177734 27.035547 29.460080 2.9056540           8
#>  7        7  -1.2867882 12.573576 14.998110 0.8706583           8
#>  8        8 -13.3622709 36.724542 38.237467 8.4753522           6
#>  9        9   1.8203869  6.359226  6.756434 0.2971458           4
#> 10       10  -1.2689999 12.538000 14.050925 0.7546570           6
#> # ... with 50 more rows
```

#### 7. Plotting predictions

When plotting non-linear fits, it often looks better to have a smooth curve, even if there are not many points underlying the fit. This can be achieved by including `newdata` in the **augment()** function and creating a higher resolution set of predictor values.

However, when predicting for many different fits, it is not certain that each curve has the same range of predictor variables. Consequently, we need to filter each new prediction by the **min()** and **max()** of the predictor variables.

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
  facet_wrap(~ growth.temp + process, labeller = labeller(.multi_line = FALSE)) +
  scale_colour_manual(values = c('green4', 'black')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab('log Metabolic rate') +
  xlab('Assay temperature (ºC)') +
  theme(legend.position = c(0.9, 0.15))
```

![](README-plot_many_fits-1.png)

#### 8. Plotting confidence intervals

The confidence intervals of each parameter for each curve fit can also be easily visualised.

``` r
# plot
ggplot(params, aes(col = flux)) +
  geom_point(aes(curve_id, estimate)) +
  facet_wrap(~ term, scale = 'free_x', ncol = 4) +
  geom_linerange(aes(curve_id, ymin = conf.low, ymax = conf.high)) +
  coord_flip() +
  scale_color_manual(values = c('green4', 'black')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  theme(legend.position = 'top') +
  xlab('curve') +
  ylab('parameter estimate')
```

![](README-confint_plot-1.png)
