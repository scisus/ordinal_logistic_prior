# functions for prior analysis

# generate samples from a truncated normal distribution.
# n = how many samples, mean = mean, sd = sd. min and max are the limits of the distribution/truncation points.
# rtnorm <- function(n, mean, sd, min, max) {
#     x <- rnorm(n, mean=mean, sd=sd)
#     x <- x[x >= min & x <= max]
#     while(length(x) < n) {
#         newx <- rnorm(1, mean=mean, sd=sd)
#         while(newx <= min | newx >= max) {
#             newx <- rnorm(1, mean=mean, sd=sd)
#         }
#         x <- c(x, newx)
#     }
#     length(x)==n
#     return(x)
# }

# Build a list of objects that are used for simulating data in stan. simu_pars is a dataframe of beta and cutpoint parameters, inputlist is what gets fed directly to stan, and h is the transition points and is useful for plotting later
set_simulation_parameters <- function(N=100) {
    N <- N
    K <- 3

    beta <- data.frame(transition = c("slow", "medium", "fast"), beta=c(0.5, 1, 2))
    cutpoints <- data.frame(c.1= c(5, 10, 20), c.2=c(7.5, 15, 30), transition=c("slow", "medium", "fast"))

    simu_pars <- merge(beta, cutpoints)

    # half transition points, engineered to be identical all transitions
    h1 <- unique(simu_pars$c.1/simu_pars$beta )
    h2 <- unique(simu_pars$c.2/simu_pars$beta )
    h <- c(h1, h2)

    # covariate over full range of heat accumulation (risto scale) Jan-Julyish
    #x <- rtnorm(n=N, mean=mean(h), sd=2, min=0, max=20) #covariate
    x <- runif(n=N, min=8, max=18)

    inputs_for_sim <- split(simu_pars, simu_pars$transition) %>%
        purrr::map(.f=function(y) {list("N" = N, "K" = K, "c" = c(y$c.1, y$c.2), "beta"=y$beta, "h" = h, "x" = x)})

    return(list(pars=simu_pars, inputlist=inputs_for_sim, h=h))
}

# simulate data from an ordinal logistic model and format it as input for stan. input is a list for stan, as is output. groups is TRUE or FALSE
simulate_data <- function(input, prior, groups) {
    # simulate data
    if (isTRUE(groups)) {
        simu <- rstan::stan(file='simulate/covar_group_sim.stan', iter=1, chains=1, algorithm="Fixed_param", data=input)
    } else {
        simu <- rstan::stan(file='simulate/covar_sim.stan', iter=1, chains=1, algorithm="Fixed_param", data=input)
    }
    # extract data from stan model
    simu_params <- rstan::extract(simu)
    # format data as input to another stan model
    input_data_for_model <- list("N" = input$N, "K" = input$K, "x" = input$x, "y" = array(simu_params$y[1,]))

    if (groups==TRUE) {
        append(input_data_for_model, "G"=input$G, "GID" = input$GID)
    }

    return(input_data_for_model)
}


# Plot simulated data to make sure it's reasonable
plot_simulated_data <- function(simdat, simulation_input) {
    # format for plotting
    simdf <- purrr::map(simdat, .f = function(x) {x[c("x", "y")]}) %>%
        purrr::map_dfr(.f = bind_rows, .id=".id")

    # scatterplot
    p1 <- ggplot(simdf, aes(x=x, y=y)) +
        geom_jitter(shape=1, height=0.1, alpha=0.5) +
        geom_vline(xintercept = simulation_input[[1]]$h) +
        ggtitle("Simulated data with cutpoints") +
        facet_grid(.id ~ .)

    # cumulative plot
    p2 <- ggplot(simdf, aes(x=x, colour=as.factor(y))) +
        stat_ecdf() +
        geom_vline(xintercept=simulation_input[[1]]$h) +
        theme(legend.position = "none") +
        ggtitle("Cumulative x for 3 states") +
        facet_grid(.id ~ .)

    cowplot::plot_grid(p1, p2, ncol=2)
}

# make shape and rate parameters for a gamma distribution centered on the mean of h (a vector) with spread scaled by a factor (scaler). Larger factors make the distribution skinnier and smaller ones make it fatter. "h" are the half transition points, so the mean is basically the midpoint of your state 2 in a 3 state system when you have data like mine.
make_gamma_pars <- function(factor, h) {
    center <- mean(h)
    shape <- mean(h) * factor
    cut_rate <- shape/center
    sr <- data.frame(shape, cut_rate)
    return(sr)
}

# turn a dataframe of parameter values into a list of parameters for each model run. Each row of the dataframe becomes a dataframe in a list
make_parframe_list <- function(parframe) {
    parlist <- split(parframe, seq(nrow(parframe)))
    names(parlist) <- parframe$modelid

    return(parlist)
}

# fit a model with a gamma prior in stan. simdatlist is a list of simulated data (1 simulated dataset per list entry), pars is a list of parameter values (1 set of parameter values per list entry) and groups is TRUE or FALSE indicating whether you're trying to fit groups.
fit_gamma_model <- function(simdatlist, pars, groups) {
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
    simdat$shape <- pars$shape
    simdat$cut_rate <- pars$cut_rate
    simdat$beta_rate <- pars$beta_rate
    
    #fit the model
    if (isTRUE(groups)) {
        fitgam <- stan(file='gamma/gamma_covar_group.stan', data=simdat, chains=4, iter = 3500, warmup=1000)
    } else {
        fitgam <- stan(file='gamma/gamma_covar.stan', data=simdat, chains=4, iter=20, warmup=10)
    }
    return(fitgam)
}

## append a label (string) to all columnnames in a dataframe (x)
label_names <- function(x, label) {
    colnames(x) <- paste0(colnames(x), "_", label)
    return(x)
}

# take a stan object and find out if there's anything egregiously wrong with it
check_model <- function(stanobj) {
    nuts <- nuts_params(stanobj)
    divergences <- dplyr::filter(nuts, Parameter=="divergent__" & Value==1) %>%
        nrow()
    rhats <- rhat(stanobj)
    bad_rhats <- sum(rhats > 1)
    nefrats <- neff_ratio(stanobj)
    bad_neff <- sum(nefrats < 0.1, is.nan(nefrats), na.rm=TRUE)
    diagnostics <- data.frame(divergences = divergences, bad_rhats=bad_rhats, bad_neff=bad_neff)
    return(diagnostics)
}

# run check_model on a list of models and add a column that names each row by the list name
check_list_of_models <- function(model_list) {
    map_dfr(model_list, check_model, .id=".id") %>%
        rename(modelid=.id)
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

HPDIlow <- function(x, prob) {
    HPDI <- rethinking::HPDI(x, prob=prob)
    return(HPDI[1])
}

HPDIhigh <- function(x, prob) {
    HPDI <- rethinking::HPDI(x, prob=prob)
    return(HPDI[2])
}

calc_HPDI <- function(params, prob) {
    low <- params %>% dplyr::summarise_at(vars(ends_with("model")), HPDIlow, prob=prob)
    high <- params %>% dplyr::summarise_at(vars(ends_with("model")), HPDIhigh, prob=prob)
    
    # awkward formatting
    hdpis <- dplyr::full_join(low, high) %>%
        select(-contains("lp"))
    colnames(hdpis) <- stringr::str_replace(colnames(hdpis), "_model", "")
    
    # true param
    true <- params %>% dplyr::summarise_at(vars(ends_with("true")), unique)
    colnames(true) <- stringr::str_replace(colnames(true), "_true", "")
    true <- select(true, colnames(hdpis))
    
    # more awkward formatting
    compframe <- dplyr::full_join(hdpis, true) %>%
        t(.) %>%
        data.frame()
    colnames(compframe) <- c("low", "high", "true")
    compframe$params <- rownames(compframe)
    
    # true param in interval?
    tf <- compframe %>% mutate(inint = true > low & true < high)
    return(tf)
    
}

# function to plot modeled parameters with label being a string to label the graph (usually prior and whether groups were included and pardf is the modeled parameter dataframe) trueparams is a one row dataframe of true parameters. h is a global parameter have fun!
parplot <- function(pars) {
    pars <- dplyr::select(pars, starts_with("c"), starts_with("beta"), starts_with("h"), starts_with("modelid")) # drop character cols
    cutplot <- mcmc_areas(pars, regex_pars="c.\\d_model") +
        geom_vline(xintercept=c(unique(pars$c.1_true), unique(pars$c.2_true)))
    
    betaplot <- mcmc_areas(pars, pars="beta_model") +
        geom_vline(xintercept=unique(pars$beta_true))
    
    h1plot <- mcmc_areas(pars, regex_pars = "h1_model") +
        geom_vline(xintercept=unique(pars$h1_true))
    h2plot <- mcmc_areas(pars, regex_pars = "h2_model") +
        geom_vline(xintercept=unique(pars$h2_true))
    allplot <- cowplot::plot_grid(cutplot, betaplot, h1plot, h2plot,
                                  nrow=1, ncol=4,
                                  labels=paste("modelid =", unique(pars$modelid_true)))
    print(allplot)
}
