
residuals <- function(add_pars, stors, times, infls, grp_pars, subs, temps) {

  if ('qhat_opt.scale' %in% names(add_pars)) {
    qscale <- add_pars['qhat_opt.scale']
    add_pars[c('qhat_opt.m1', 'qhat_opt.m2', 'qhat_opt.m3', 'qhat_opt.m4', 'qhat_opt.m5')] <- log10(grp_pars[['qhat_opt']][c('m1', 'm2', 'm3', 'm4', 'm5')]) + qscale
  }
  add_pars <- 10^add_pars
  add_pars <- as.list(add_pars)

  abm_out <- foreach (i = ids, .combine = rbind, .multicombine = TRUE, .packages = 'data.table') %dorng% {
    out <- abmneo(
      storage = stors[[i]],
      days = 365,
      times = times[[i]],
      inf_pars = infls[[i]],
      grp_pars = grp_pars,
      sub_pars = subs[[i]],
      var_pars = list(var = temps[[i]]),
      add_pars = add_pars,
      startup = startup,
      quiet = TRUE
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

  print(sum(abs(resids)))
  print(sum(resids^2))

  #return(sum(abs(resids)))
  return(sum(resids^2))

}
