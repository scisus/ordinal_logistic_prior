#functions for prior_analysis_groups.Rmd
#
# Set parameters for a simulation in stan. N is how many observations, hmod is a vector of length # groups describing how much each half transition is shifted. Number of groups is set to 7. The output is a list of three
# simu_pars: a dataframe of parameters to be used in the simulation
# inputs_for_sim: a Stan-readable list of information in simu_pars combined with a covariate
# h_base: the population level transition
set_simulation_parameters <- function(N=100, hmod, G=G, noeffect=TRUE) {
  K <- 3

  beta <- data.frame(transition = c("medium", "fast"), beta=c(1, 2))
  cutpoints <- data.frame(c.1= c(10, 20), c.2=c(15, 30), transition=c("medium", "fast"))

  h_base <- data.frame(h.1 = 10, h.2 = 15) # half transitions

  groups <- data.frame(group = 1:G, h_mod = hmod)

  simu_pars <- merge(beta, cutpoints) %>%
    merge(h_base) %>%
    merge(groups) %>%
    mutate(h.1_group = h.1 + h_mod, h.2_group = h.2 + h_mod) %>%
    mutate(alpha_group = (beta*(h.1_group) - c.1))


  # covariate around active period
  x <- runif(n=N, min=8, max=18)

  if (isFALSE(noeffect)) {
    inputs_for_sim <- split(simu_pars, simu_pars$transition) %>%
      purrr::map(.f=function(y) {list("N" = N, "K" = K, "G" = G, "c" = unique(c(y$c.1, y$c.2)), "beta"=unique(y$beta), "alpha_g" = unique(y$alpha_group), "group" = unique(y$group), "x" = x)})
  }

  if (isTRUE(noeffect)) {
    inputs_for_sim <- split(simu_pars, simu_pars$transition) %>%
      purrr::map(.f=function(y) {list("N" = N, "K" = K, "G" = G, "c" = unique(c(y$c.1, y$c.2)), "beta"=unique(y$beta), "alpha_g" = rep(0, G), "group" = unique(y$group), "x" = x)})
  }

  return(list(pars=simu_pars, inputlist=inputs_for_sim, h=h_base))
}


# simulate from models with identical transition points. Choose whether effects are on alpha (intercept) or beta (slope)

simulate_data <- function(input) {
  # simulate data
    simu <- rstan::stan(file='simulate/covar_group_alpha_sim.stan', iter=1, chains=1, algorithm="Fixed_param", data=input)

  # extract data from stan model
  simu_params <- data.frame(rstan::extract(simu)$y[1,,]) %>%
    mutate(x = input$x) %>%
    pivot_longer(cols = -x, names_to = "group", values_to = "y") %>%
    extract(group, into = "group", regex = "([:digit:]{1,2})") %>%
    mutate(group = as.integer(group))

  # format data as input to another stan model

  input_data_for_model <- list("N" = input$N * input$G, "K" = input$K, "G" = input$G, "x" = simu_params$x, "y" = simu_params$y, "group" = simu_params$group)

  return(input_data_for_model)
}

# Plot simulated data to make sure it's reasonable
plot_simulated_data <- function(simdat, simulation_input) {
  # format for plotting
  simdf <- purrr::map(simdat, .f = function(x) {x[c("x", "y", "group")]}) %>%
    purrr::map_dfr(.f = bind_rows, .id=".id")

  # scatterplot
  p1 <- ggplot(simdf, aes(x=x, y=y)) +
    geom_jitter(shape=1, height=0.1, alpha=0.5) +
    ggtitle("Simulated data") +
    facet_grid(group ~ .id)

  # cumulative plot
  p2 <- ggplot(simdf, aes(x=x, colour=as.factor(y))) +
    stat_ecdf() +
    theme(legend.position = "none") +
    ggtitle("Cumulative x for 3 states") +
    facet_grid(group ~ .id)

  cowplot::plot_grid(p1, p2, nrow=2)
}

# format parframe so it works with parLapply better
make_parframe_list <- function(parframe) {
  parlist <- split(parframe, seq(nrow(parframe)))
  names(parlist) <- parframe$modelid

  return(parlist)
}

# fit a model with an induced dirichlet prior on cutpoints in stan. simdatlist is a list of simulated data (1 simulated dataset per list entry), pars is a list of parameter values (1 set of parameter values per list entry) and groups is TRUE or FALSE indicating whether you're trying to fit groups.
fit_indir_model <- function(simdatlist, pars) {
  #choose whether to use data simulated with a rapid or slow transition
  if (pars$transition == "medium") {
    simdat <- simdatlist$medium
  }
  if (pars$transition == "fast") {
    simdat <- simdatlist$fast
  }
  #extract parameters for prior distributions
  simdat$anchor <- pars$anchor
  simdat$beta_rate <- pars$beta_rate

  #fit the model

  fitindir <- stan(file='induced_dirichlet/dirichlet_covar_alpha_group.stan', data=simdat, chains=4)

  return(fitindir)
}


## append a label (string) to all columnnames in a dataframe (x)
label_names <- function(x, label) {
  colnames(x) <- paste0(colnames(x), "_", label)
  return(x)
}

# make_long_param <- function(params) {
#   params <- params %>%
#     pivot_longer(starts_with("alpha"), names_to = "group", values_to = "alpha_group") %>%
#       extract(group, regex="([:digit:])", into = "group")
#   return(params)
# }
#
table_problems <- function(problemlist) {
  problems <- problemlist$problems %>%
    map_dfr(bind_rows, .id=".id") %>%
    rename(run=.id) %>%
    group_by(modelid) %>%
    summarize(bad_proportion = n()/reps, bad_count = n(), divergences_mean=mean(divergences), bad_rhats_mean=mean(bad_rhats), bad_neff_mean=mean(bad_neff))

  return(problems)
}



# calculate half transition points with group effects
calc_h <- function(pars) {

    newpars <- pars %>%
      mutate(h.1_group.1 = (c.1 + alpha_g.1)/beta,
             h.1_group.2 = (c.1 + alpha_g.2)/beta,
             h.1_group.3 = (c.1 + alpha_g.3)/beta,
             h.1_group.4 = (c.1 + alpha_g.4)/beta,
             h.1_group.5 = (c.1 + alpha_g.5)/beta,
             h.1_group.6 = (c.1 + alpha_g.6)/beta,
             h.1_group.7 = (c.1 + alpha_g.7)/beta,
             h.1_group.8 = (c.1 + alpha_g.8)/beta,
             h.1_group.9 = (c.1 + alpha_g.9)/beta,
             h.1_group.10 = (c.1 + alpha_g.10)/beta,


             h.2_group.1 = (c.2 + alpha_g.1)/beta,
             h.2_group.2 = (c.2 + alpha_g.2)/beta,
             h.2_group.3 = (c.2 + alpha_g.3)/beta,
             h.2_group.4 = (c.2 + alpha_g.4)/beta,
             h.2_group.5 = (c.2 + alpha_g.5)/beta,
             h.2_group.6 = (c.2 + alpha_g.6)/beta,
             h.2_group.7 = (c.2 + alpha_g.7)/beta,
             h.2_group.8 = (c.2 + alpha_g.8)/beta,
             h.2_group.9 = (c.2 + alpha_g.9)/beta,
             h.2_group.10 = (c.2 + alpha_g.10)/beta)

  return(newpars)
}

# calc_true_h <- function(truepars) {
#   newtruepars <- truepars %>%
#     mutate(h.1_group.1 = (c.1 + alpha_g.1)/beta,
#            h.1_group.2 = (c.1 + alpha_g.2)/beta,
#            h.1_group.3 = (c.1 + alpha_g.3)/beta,
#            h.1_group.4 = (c.1 + alpha_g.4)/beta,
#            h.1_group.5 = (c.1 + alpha_g.5)/beta,
#            h.1_group.6 = (c.1 + alpha_g.6)/beta,
#            h.1_group.7 = (c.1 + alpha_g.7)/beta,
#
#            h.2_group.1 = (c.2 + alpha_g.1)/beta,
#            h.2_group.2 = (c.2 + alpha_g.2)/beta,
#            h.2_group.3 = (c.2 + alpha_g.3)/beta,
#            h.2_group.4 = (c.2 + alpha_g.4)/beta,
#            h.2_group.5 = (c.2 + alpha_g.5)/beta,
#            h.2_group.6 = (c.2 + alpha_g.6)/beta,
#            h.2_group.7 = (c.2 + alpha_g.7)/beta)
#   return(newtruepars)
# }

#bind true parameters (in list parlist) and model parameters (in list fit) even tho it will make a giant df.
bind_true_model_pars <- function(fits, parlist) {
  # extract params from model object
  params <- lapply(fits, function(x) {data.frame(rstan::extract(x) ) } )

  # calculate transitions
  params <- purrr::map(params, calc_h)
  parlist <- purrr::map(parlist, calc_h)

  # label params as coming from the model or as true params used to
  params <- purrr::map(params, label_names, label="model")
  parlist <- purrr::map(parlist, label_names, label="true")

  # combine model and true params in a big list of dataframes - each list entry is a dataframe for a single model
  params <- purrr::map2(params, parlist, cbind)

  return(params)
}

# take a stan object and find out if there's anything egregiously wrong with it
check_model <- function(stanobj) {
  nuts <- nuts_params(stanobj)
  divergences <- dplyr::filter(nuts, Parameter=="divergent__" & Value==1) %>%
    nrow()
  rhats <- rhat(stanobj)
  bad_rhats <- sum(rhats > 1.01 | is.na(rhats))
  nefrats <- neff_ratio(stanobj)
  bad_neff <- sum(nefrats < 0.1, is.nan(nefrats), na.rm=TRUE)
  diagnostics <- data.frame(divergences = divergences, bad_rhats=bad_rhats, bad_neff=bad_neff)
  return(diagnostics)
}

# run check_model on a list of models and add a column that names each row by the list name
check_list_of_models <- function(model_list) {
  purrr::map_dfr(model_list, check_model, .id=".id") %>%
    rename(modelid=.id)
}

# Extract model configurations and the most obvious problems with fit/convergence from stanfit objects. Stanfit objects are in lists in .rds files saved in path. path is a string denoting the directory where the rds files are stored. Nothing other than .rds files with lists of stanfit objects should be stored in path.
extract_pars_and_problems <- function(path, parlist) {
  fits <- list.files(path=path)
  reps <- length(fits)
  pars <- list()
  bad_models <- list()

  for (i in 1:reps) {
    run <- readRDS(paste0(path, fits[i])) # read in first run
    run_pars <- bind_true_model_pars(run, parlist=parlist) # extract parameters
    run_fit <- check_list_of_models(run) %>% # id problems
      mutate(all_bads = divergences + bad_rhats + bad_neff) %>%
      filter(all_bads > 0) %>%
      select(-all_bads)
    pars[[i]] <- run_pars
    bad_models[[i]] <- run_fit
    rm(run, run_pars, run_fit)
    gc()
  }

  return(list(pars = pars, problems = bad_models))
}


# parplot <- function(modelpars) {
#   p1 <- ggplot(modelpars, aes(x=beta_model, color=run)) +
#     geom_density()+
#     theme(legend.position="none") +
#     xlab("beta") +
#     geom_vline(xintercept = unique(modelpars$beta_true))
#
#   p2 <- ggplot(modelpars, aes(x=c.1_model, color=run)) +
#     geom_density() +
#     geom_density(aes(x=c.2_model)) +
#     theme(legend.position="none") +
#     xlab("cutpoints") +
#     geom_vline(xintercept = c(unique(modelpars$c.1_true), unique(modelpars$c.2_true)))
#
#   p3 <- ggplot(modelpars, aes(x=h.1_group.1_model, color=run)) +
#     geom_density() +
#     geom_density(aes(x=h.2_group.1_model)) +
#     theme(legend.position="none") +
#     xlab("group 1 transitions") +
#     geom_vline(xintercept = c(unique(modelpars$h.1_group.1_true), unique(modelpars$h.2_group.1_true)))
#
#   cowplot::plot_grid(p1, p2, p3, labels=paste("model", modelpars$modelid_true))
#
# }

plot_densities_anchorxbeta <- function(data, modelpar, truepar) {
  plot <- ggplot(data, aes(x=get(modelpar), group=interaction(run, modelid_true), fill=transition_true, colour=transition_true)) +
    geom_density(alpha=0.1) +
    geom_vline(aes(xintercept=get(truepar), colour = transition_true), alpha=1) +
    facet_grid(anchor_true ~ beta_rate_true) +
    ggtitle(strsplit(modelpar, split = "_model"), subtitle = "facet rows are anchor, cols are beta prior param")
  print(plot)
}


# calculate whether true value is in HPDI
#
HPDIlow <- function(x, prob) {
  HPDI <- rethinking::HPDI(x, prob=prob)
  return(HPDI[1])
}

HPDIhigh <- function(x, prob) {
  HPDI <- rethinking::HPDI(x, prob=prob)
  return(HPDI[2])
}

calc_HPDI <- function(params, prob) {
  params <- params %>%
    tidyr::pivot_longer(ends_with("model"), names_to = "param", values_to = "param_value") %>%
    filter(param != "lp___model") %>%
    tidyr::extract(param, into="param", regex="(.*)_model") %>% # drop model ending
    rename(modelid=modelid_true)

  # calculate bottom and top of HPDI at prob for each model run and each parameter
  low <- params %>%
    group_by(modelid, run, param) %>%
    summarise(low= HPDIlow(param_value, prob=prob))
  high <- params %>%
    group_by(modelid, run, param) %>%
    summarise(high= HPDIhigh(param_value, prob=prob))


  hdpis <- dplyr::full_join(low, high)

  # true params
  true <- params %>% dplyr::summarise_at(vars(ends_with("true")), unique)
  colnames(true) <- stringr::str_replace(colnames(true), "_true", "")
  true <- select(true, c("c.1", "c.2", "beta", starts_with("alpha"), starts_with("h."))) %>%
    tidyr::pivot_longer(cols=c("c.1", "c.2", "beta", starts_with("alpha"), starts_with("h.")), names_to = "param", values_to="true")

  # combine true and hdpis for comparison
  compframe <- full_join(hdpis, true)


  # true param in interval?
  tf <- compframe %>% mutate(inint = true > low & true < high)
  return(tf)
}

# For each parameter (5) in each model (27) in each run (30), is a given parameter in the 50% or 90% HPDI? modelpars is a list of dataframes. Each dataframe in the list contains parameters from 30 runs of 1 model along with the true parameters used to simulate the datasets.
# output is a list of 2 dataframes 4050 rows each (5x27x30) - with each parameter estimated by the model and an inint column with TRUE if it falls in the HPDI interval and FALSE if not.
which_params_recaptured <- function(modelpars) {
  in50 <- purrr::map(modelpars, calc_HPDI, prob=0.5)

  #in75 <- purrr::map(params_indir, calc_HPDI, prob=0.75)

  in90 <- purrr::map(modelpars, calc_HPDI, prob=.90)

  # recaptured parameters (refactor this later so it's all in one df
  perform50 <- purrr::map_dfr(in50, bind_rows)
  perform90 <- purrr::map_dfr(in90, bind_rows)

  return(list(fifty = perform50, ninety=perform90))
}

calc_prop_recaptured_overall <- function(inint, truepars) {

  # trues <- inint$fifty %>%
  #   ungroup() %>%
  #   select(modelid, param, true) %>%
  #   distinct()
  #
  # nas <- which(is.na(inint$fifty$inint))
  # inint$fifty[nas,]

  # proportion of parameters recaptured
  prop_recaptured50 <- inint$fifty %>%
    filter(param != "sigma_group") %>%
    group_by(modelid, run)  %>%
    summarise(captured = sum(inint, na.rm = TRUE))  %>%
    summarise(mean_captured = mean(captured), sd_captured=sd(captured)) %>%
    full_join(truepars) %>%
    arrange(beta, desc(mean_captured))


  prop_recaptured90 <- inint$ninety %>%
    filter(param != "sigma_group") %>%
    group_by(modelid, run) %>%
    summarise(captured = sum(inint, na.rm=TRUE)) %>%
    summarise(mean_captured = mean(captured), sd_captured=sd(captured)) %>%
    full_join(truepars) %>%
    arrange(beta, desc(mean_captured))

  return(list(fifty = prop_recaptured50, ninety = prop_recaptured90))
}

calc_recaptured_by_param <- function(inint, truepars) {
  perform50_summary <- inint$fifty %>%
    filter(param != "sigma_group") %>%
    group_by(modelid, param) %>%
    summarise(prop_inint = mean(inint)) %>%
    left_join(truepars)

  perform90_summary <- inint$ninety %>%
    filter(param != "sigma_group") %>%
    group_by(modelid, param) %>%
    summarise(prop_inint = mean(inint)) %>%
    left_join(truepars)

  return(list(fifty=perform50_summary, ninety=perform90_summary))
}
