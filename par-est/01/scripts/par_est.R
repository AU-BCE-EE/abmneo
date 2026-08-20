
# Parallel setup
n_cores <- min(maxcores, parallel::detectCores() - 2)  # Leave at least 2 cores free
if (!foreach::getDoParRegistered()) {
  registerDoParallel(cores = n_cores)
}

# Reproducible RNG
seed_main <- 123
registerDoRNG(seed_main)

res <- optim(
  par = c(h_rate_ref.PS = -1.3, qhat_opt.m1 = -1, qhat_opt.m2 = -0.1),
  fn = residuals,
  stors = stors,
  times = times, 
  infls = infls, 
  grp_pars = grp_pars, 
  subs = subs, 
  temps = temps
)
