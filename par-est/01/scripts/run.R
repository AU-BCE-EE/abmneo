# Run ABM with best-fit pars

abm_out <- data.table()
for (i in ids) {
  cat(i, '\n')
  out <- abmneo(
    storage = stors[[i]],
    days = 365,
    times = 0:365,
    inf_pars = infls[[i]],
    grp_pars = grp_pars,
    sub_pars = subs[[i]],
    var_pars = list(var = temps[[i]]),
    add_pars = best_pars,
    startup = 1
  )
  setDT(out)
  out[, tank := i]
  abm_out <- rbind(abm_out, out)
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
  by.y = c('tank', 'doy_end_meas'),
  all = TRUE
)

ggplot(emis_comp, aes(doy_end_abm, CH4_emis_rate_mod, colour = tank)) +
  geom_step() +
  geom_point(aes(x = doy_mid, y = CH4_emis_rate_meas)) +
  facet_wrap(~ tank, scale = 'free') +
  theme_bw()
ggsave('../plots/emis_comp.png', height = 4, width = 6)

