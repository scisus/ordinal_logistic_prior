# meta for ordinal_logistic_prior

This folder contains files where I try to figure out what kind of prior to use on cutpoints in an ordinal logistic model.

The induced dirichlet prior was developed by [Michael Betancourt]](https://betanalpha.github.io/assets/case_studies/ordinal_regression.html)

# folders contain Stan code
* `gamma` contains Stan code for models with gamma priors on the cutpoints
* `induced_dirichlet` contains Stan code for models with induced dirichlet priors on the cutpoints
* `simulate` contains Stan code that simulates data from ordinal logistic models

# files
* `prior_analysis.Rmd` Document with analysis of different priors and description of results
* `prior_analysis_functions.R` helper functions
* `prior_analysis_gamma.Rmd` analysis of gamma prior on cutpoints
* `prior_analysis_induced_dirichlet.Rmd` analysis of induced dirichlet prior on cutpoints
* `tests` tests for functions in `prior_analysis_functions.R`

