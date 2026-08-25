
# Parallel setup
if (!foreach::getDoParRegistered()) {
  registerDoParallel(cores = min(maxcores, parallel::detectCores() - 2))
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
  par = c(h_rate_ref.PS = log10(0.003), resid_enrich = log10(0.9)),
  fn = residuals,
  stors = stors,
  times = times, 
  infls = infls, 
  grp_pars = grp_pars, 
  subs = subs, 
  temps = temps,
  ids = c('AD1', 'AD2', 'AD3'),
  obj = 'ss'
)

best_pars <- as.list(10^res$par)
