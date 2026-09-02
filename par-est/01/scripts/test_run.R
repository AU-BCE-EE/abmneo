# Test ABM runs
if (!foreach::getDoParRegistered()) {
  registerDoParallel(cores = min(maxcores, parallel::detectCores() - 2))
}

abm_out <- foreach (i = ids, .combine = rbind, .multicombine = TRUE, .packages = 'data.table') %dorng% {
  out <- abmneo(
    storage = stors[[i]],
    days = 365,
    times = times[[i]],
    inf_pars = infls[[i]],
    grp_pars = grp_pars,
    sub_pars = subs[[i]],
    var_pars = list(var = temps[[i]]),
    startup = 1
  )
  setDT(out)
  out[, tank := i]
  out
}

# For clarity
abm_out[, CH4_emis_rate_mod := CH4_emis_rate_ave / 1000]
abm_out[, doy_end_abm := time]
emis[, CH4_emis_rate_meas := CH4_emis_rate]
emis[, doy_end_meas := doy_end]

# Merge
emis_comp <- merge(
  abm_out[, .(tank, time, doy_end_abm, CH4_emis_rate_mod)],
  emis[, .(tank, doy_mid, doy_end_meas, CH4_emis_rate_meas)],
  by.x = c('tank', 'doy_end_abm'),
  by.y = c('tank', 'doy_end_meas')
)

# Residuals
resids <- emis_comp$CH4_emis_rate_mod - emis_comp$CH4_emis_rate_meas

ggplot(emis_comp, aes(doy_end_abm, CH4_emis_rate_mod, colour = tank)) +
  geom_step() +
  geom_point(aes(x = doy_mid, y = CH4_emis_rate_meas)) +
  facet_wrap(~ tank) +
  theme_bw()
ggsave('../plots/test_run.png')

ggplot(abm_out, aes(time, slurry_depth, colour = tank)) +
  geom_step() +
  facet_wrap(~ tank) +
  theme_bw()
ggsave('../plots/slurry_level.png')


