# functions supporting priorAnalysisGamma.Rmd


# fit model with induced dirichlet prior
fit_indir_model <- function(simdatlist, pars, groups) {
    #choose whether to use data simulated with a rapid or slow transition
    if (pars$transition == "slow") {
        simdat <- simdatlist$slow
    }
    if (pars$transition == "medium") {
        simdat <- simdatlist$medium
    }
    if (pars$transition == "fast") {
        simdat <- simdatlist$fast
    }

    #extract parameters for prior distribtuions
    simdat$beta_rate <- pars$beta_rate
    simdat$anchor <- pars$anchor

    #fit the model

    fit <- stan(file='dirichlet prior/dirichlet/dirichlet_covar.stan', data=simdat, chains=4)
    return(fit)
}



#bind true parameters (in list parlist) and model parameters (in list fit) even tho it will make a giant df.
bind_true_model_pars <- function(fits, parlist) {
    # extract params from model object
    params <- lapply(fits, function(x) {data.frame(rstan::extract(x) ) } )

    # label params as coming from the model or as true params used to
    params <- map(params, label_names, label="model")
    parlist <- map(parlist, label_names, label="true")

    # combine model and true params in a big list of dataframes - each list entry is a dataframe for a single model
    params <- map2(params, parlist, cbind)

    return(params)
}

## append a label (string) to all columnnames in a dataframe (x)
label_names <- function(x, label) {
    colnames(x) <- paste0(colnames(x), "_", label)
    return(x)
}



# function to calculate the difference between modeled parameters (in the model_params dataframe) and true parameters (h is globally declared). Where true params is a dataframe
posterior_differencer <- function(pars, h) {
    c1_diff <- pars$c.1_model - pars$c.1_true
    c2_diff <- pars$c.2_model - pars$c.1_true
    h1_diff <- pars$h1_model - pars$h1_true
    h2_diff <- pars$h2_model - pars$h2_true
    beta_diff <- pars$beta_model - pars$beta_true
    diffframe <- data.frame(c1_diff, c2_diff, h1_diff, h2_diff, beta_diff)
    return(diffframe)
}

# function to plot histograms of differences between true params and modeled params (model_params dataframe). MAKE IT WORK WITH BOTH GAM AND INDIR
diffplotter <- function(diffs, pars) {
    #diffs <- posterior_differencer(pars)
    cuts <- mcmc_intervals(diffs, regex_pars = "c") +
        ggtitle("", subtitle = "differences between modeled and true params")
    opars <- mcmc_intervals(diffs, pars=c("h1_diff", "h2_diff", "beta_diff"))
    cowplot::plot_grid(cuts, opars, labels=paste("model", unique(pars$modelid_true),
                                                 "rate=", unique(pars$beta_rate_true),
                                                 "anchor=", unique(pars$anchor_true)))
}







