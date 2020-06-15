#functions for prior_analysis_groups.Rmd
#
# Set parameters for a simulation in stan. N is how many observations, hmod is a vector of length # groups describing how much each half transition is shifted. Number of groups is set to 7. The output is a list of three
# simu_pars: a dataframe of parameters to be used in the simulation
# inputs_for_sim: a Stan-readable list of information in simu_pars combined with a covariate
# h_base: the population level transition
set_simulation_parameters <- function(N=100, hmod=c( 0.35, -0.52, -1.32, -0.72, -0.78,  0.73, -0.45)) {
K <- 3
G <- 7 # number of groups (sites)
N <- N

beta <- data.frame(transition = c("medium", "fast"), beta=c(1, 2))
cutpoints <- data.frame(c.1= c(10, 20), c.2=c(15, 30), transition=c("medium", "fast"))

h_base <- data.frame(h.1 = 10, h.2 = 15) # half transitions

groups <- data.frame(group = 1:G, h_mod = h_mod)

simu_pars <- merge(beta, cutpoints) %>%
  merge(h_base) %>%
  merge(groups) %>%
  mutate(h.1_group = h.1 + h_mod, h.2_group = h.2 + h_mod) %>%
  mutate(beta_group = (c.1/(h.1_group) - beta), alpha_group = (beta*(h.1_group) - c.1))


# covariate around active period
x <- runif(n=N, min=8, max=18)

inputs_for_sim <- split(simu_pars, simu_pars$transition) %>%
  purrr::map(.f=function(y) {list("N" = N, "K" = K, "G" = G, "c" = unique(c(y$c.1, y$c.2)), "beta"=unique(y$beta), "beta_g" = unique(y$beta_group), "alpha_g" = unique(y$alpha_group), "group" = unique(y$group), "x" = x)})

return(list(pars=simu_pars, inputlist=inputs_for_sim, h=h_base))
}
