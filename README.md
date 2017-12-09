<!-- README.md is generated from README.Rmd. Please edit that file -->
nls.multstart
-------------

[![Travis-CI Build Status](https://travis-ci.org/padpadpadpad/nls.multstart.svg?branch=master)](https://travis-ci.org/padpadpadpad/nls.multstart)

### Overview

**nls.multstart** is an R package that allows more robust and reproducible non-linear regression compared to **nls** or **nlsLM**. These functions allow only a single starting value, meaning that it can be hard to get the best estimated model. This is especially true if the same model is fitted over the levels of a factor, which may have the same shape of curve, but be much different in terms of parameter estimates.

**nls\_multstart()** is the main (currently only) function of **nls.multstart**. Similar to the R package **nls2**, it allows multiple starting values for each parameter and then iterates through multiple starting values, attempting a fit with each set of start parameters. The best model is then picked on AIC score. This results in a more reproducible and reliable method of fitting non-linear least squares regression in R.

This package is designed to work with the **tidyverse**, harnessing the functions within **broom**, **tidyr**, **dplyr** and **purrr** to extract estimates and plot things easily with **ggplot2**. A slightly less tidy-friendly implementation is [**nlsLoop**](https://github.com/padpadpadpad/nlsLoop).

### Issues and suggestions

Please report any issues/suggestions for improvement in the [issues link](https://github.com/padpadpadpad/nls.multstart/issues) for the repository. Or please email <d.padfield@exeter.ac.uk>.
