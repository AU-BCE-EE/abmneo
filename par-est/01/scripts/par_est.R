
# Parallel setup
n_cores <- min(maxcores, parallel::detectCores() - 2)  # Leave at least 2 cores free
if (!foreach::getDoParRegistered()) {
  registerDoParallel(cores = n_cores)
}

## Test
#residuals(add_pars = c(h_rate_ref.PS = -1.3), stors, times, infls, grp_pars, subs, temps)
#residuals(add_pars = c(h_rate_ref.PS = -1.3, qhat_opt.scale = 0), stors, times, infls, grp_pars, subs, temps)
#residuals(add_pars = c(h_rate_ref.PS = -1.3, qhat_opt.scale = 0.4), stors, times, infls, grp_pars, subs, temps)
#residuals(add_pars = c(h_rate_ref.PS = -1.3, qhat_opt.scale = -1), stors, times, infls, grp_pars, subs, temps)

# Reproducible RNG
seed_main <- 123
registerDoRNG(seed_main)

res <- optim(
  par = c(h_rate_ref.PS = -1.3, h_rate_q10.PS = log10(2), qhat_opt.scale = 0),
  fn = residuals,
  stors = stors,
  times = times, 
  infls = infls, 
  grp_pars = grp_pars, 
  subs = subs, 
  temps = temps
)

best_pars <- as.list(10^res$par)
