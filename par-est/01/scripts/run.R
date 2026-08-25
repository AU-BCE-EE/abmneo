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

names(emis)
ggplot(abm_out, aes(time, CH4_emis_rate / 1000, colour = tank)) +
  geom_step() +
  geom_point(data = emis, aes(x = doy_mid, y = CH4_emis_rate)) +
  facet_wrap(~ tank, scale = 'free') +
  theme_bw()
ggsave('../plots/emis_comp.png', height = 4, width = 6)


ggplot(abm_out, aes(time, VFA_conc, colour = tank)) +
  geom_step() +
  facet_wrap(~ tank, scale = 'free') +
  theme_bw()
ggsave('../plots/VFA_conc.png', height = 4, width = 6)

ggplot(abm_out, aes(time, PS_conc, colour = tank)) +
  geom_step() +
  facet_wrap(~ tank, scale = 'free') +
  theme_bw()
ggsave('../plots/PS_conc.png', height = 4, width = 6)

grp_out <- melt(abm_out, id.vars = c('tank', 'time'), measure.vars = c('m1_conc', 'm2_conc', 'm3_conc', 'm4_conc', 'm5_conc'))
ggplot(grp_out, aes(time, value, colour = variable)) +
  geom_step() +
  facet_wrap(~ tank, scale = 'fixed') +
  theme_bw()
ggsave('../plots/mic_conc.png', height = 4, width = 6)

