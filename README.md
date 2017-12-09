<!-- README.md is generated from README.Rmd. Please edit that file -->
nlsiter
-------

### Overview

**nlsiter** is an R package that allows more robust and reproducible fitting than vanilla nls or even nlsLM. These functions only allow a single starting value, meaning that it can be hard to get the best estimated model from a single set of starting values. This is especially true if the same model is wanting to be fitted over the levels of a factor.

**nls\_iter()** is the main (currently only) function of **nlsiter**. Similar to the R package **nls2**, it allows multiple starting values for each parameter and then tries multiple starting values on the data. The best model is then picked on AIC score. This results in a more reproducible and reliable method of fitting non-linear least squares regression in R.

This package is designed to work with the **tidyverse**, harnessing the functions within **broom**, **tidyr**, **dplyr** and **purrr** to extract estimates and plot things easily with **ggplot2**. A slightly less tidy-friendly implementation is [**nlsLoop**](https://github.com/padpadpadpad/nlsLoop)

### Issues and suggestions

Please report any issues/suggestions for improvement in the [issues link](https://github.com/padpadpadpad/nlsLoop/issues) for the repository. Or please email <d.padfield@exeter.ac.uk>.
