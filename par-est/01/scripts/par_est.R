# Parameter estimation

# Parallel setup
if (!foreach::getDoParRegistered()) {
  registerDoParallel(cores = min(maxcores, parallel::detectCores() - 2))
}

## Test
#residuals(add_pars = c(h_rate_ref.PS = -1.3), stors, times, infls, grp_pars, subs, temps, meas = emis, ids = c('AD1', 'AD2', 'AD3'), weights = 'weight', obj = 'sae')
#residuals(add_pars = c(h_rate_ref.PS = -2), stors, times, infls, grp_pars, subs, temps, meas = emis, ids = c('AD1', 'AD2', 'AD3'), weights = 'weight', obj = 'sae')

# Weight by inverse of mean emission (so smaller tanks don't count less)
emis[, weight := 1 / mean(CH4_emis_rate), by = tank]

# Reproducible RNG
seed_main <- 123
registerDoRNG(seed_main)

res <- optim(
  par = c(h_rate_ref.PS = log10(0.03), resid_enrich = log10(0.9)),
  fn = residuals,
  stors = stors,
  times = times, 
  infls = infls, 
  grp_pars = grp_pars, 
  subs = subs, 
  temps = temps,
  meas = emis,
  ids = c('AD1', 'AD2', 'AD3'),
  weights = 'weight',
  obj = 'sae'
)

best_pars <- as.list(10^res$par)
