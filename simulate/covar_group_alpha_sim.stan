// Simulate ordinal logistic data with a covariate and group effect

data {
  int<lower=1> N;
  int<lower=2> K;
  int<lower=1> G;

  positive_ordered[K-1] c; //cutpoints
  real beta; //slope

  vector[N] x; //covariate

  int group[G]; // group ids
  vector[G] alpha_g; //group effects

}

generated quantities {
  int<lower=1, upper=K> y[N,G];                     // Simulated ordinals
  matrix[N,G] gamma; //simulated linear model

  for (n in 1:N) {
    for (g in 1:G) {
      gamma[n,g] = (x[n] * beta) + alpha_g[group[g]]; // Latent effect
      y[n,g] = ordered_logistic_rng(gamma[n,g], c); // ordinals
    }
  }
}

