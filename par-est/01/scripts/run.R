# Run ABM with best-fit pars

abm_out <- foreach (i = ids, .combine = rbind, .multicombine = TRUE, .packages = 'data.table') %dorng% {
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
  out
}

abm_out_ave <- foreach (i = ids, .combine = rbind, .multicombine = TRUE, .packages = 'data.table') %dorng% {
  out <- abmneo(
    storage = stors[[i]],
    days = 365,
    times = times[[i]],
    inf_pars = infls[[i]],
    grp_pars = grp_pars,
    sub_pars = subs[[i]],
    var_pars = list(var = temps[[i]]),
    add_pars = best_pars,
    startup = 1
  )
  setDT(out)
  out[, tank := i]
  out
}

# For clarity
abm_out_ave[, CH4_emis_rate_mod := CH4_emis_rate_ave / 1000]
abm_out_ave[, doy_end_abm := time]
emis[, CH4_emis_rate_meas := CH4_emis_rate]
emis[, doy_end_meas := doy_end]

# Merge
emis_comp <- merge(
  abm_out_ave[, .(tank, time, doy_end_abm, CH4_emis_rate_mod)],
  emis[, .(tank, doy_mid, doy_end_meas, CH4_emis_rate_meas)],
  by.x = c('tank', 'doy_end_abm'),
  by.y = c('tank', 'doy_end_meas'),
  all = TRUE
)

# plot high-res mode output first
ggplot(abm_out, aes(time, CH4_emis_rate / 1000, colour = tank)) +
  geom_step() +
  geom_point(data = emis, aes(x = doy_mid, y = CH4_emis_rate)) +
  facet_wrap(~ tank, scale = 'free') +
  theme_bw() +
  labs(x = 'Time (d)', y = 'CH4 emis. rate (g/d)')
ggsave('../plots/emis_comp.png', height = 4, width = 6)

# Then average, which was used for par est
ggplot(emis_comp, aes(doy_end_abm, CH4_emis_rate_mod, colour = tank)) +
  geom_step() +
  geom_point(aes(x = doy_mid, y = CH4_emis_rate_meas)) +
  facet_wrap(~ tank) +
  theme_bw() +
  labs(x = 'Time (d)', y = 'Meas. period ave. CH4 emis. (g/d)')
ggsave('../plots/emis_comp_ave.png', height = 4, width = 6)

ggplot(abm_out, aes(time, 1000 * VFA_conc, colour = tank)) +
  geom_step() +
  facet_wrap(~ tank, scale = 'fixed') +
  theme_bw() +
  labs(x = 'Time (d)', y = 'Total VFA (mg/kg)')
ggsave('../plots/VFA_conc.png', height = 4, width = 6)

ggplot(abm_out, aes(time, PS_conc, colour = tank)) +
  geom_step() +
  facet_wrap(~ tank, scale = 'free') +
  theme_bw() +
  labs(x = 'Time (d)', y = 'Substrate conc. (g/kg)')
ggsave('../plots/PS_conc.png', height = 4, width = 6)

grp_out <- melt(abm_out, id.vars = c('tank', 'time'), measure.vars = c('m1_conc', 'm2_conc', 'm3_conc', 'm4_conc', 'm5_conc'))
ggplot(grp_out, aes(time, value, colour = variable)) +
  geom_step() +
  facet_wrap(~ tank, scale = 'fixed') +
  theme_bw() +
  labs(x = 'Time (d)', y = 'Methanogen conc. (g/kg)')
ggsave('../plots/mic_conc.png', height = 4, width = 6)

