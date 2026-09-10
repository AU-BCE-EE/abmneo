inf_pars_pig_01 <- list(
    VFA_fresh = 1.7,
    VFA_init = 1.7,
    comps = c('SO4'),
    comp_fresh = c(SO4 = 0.01),
    comp_init =  c(SO4 = 0.01),
    pH = 7,
    dens = 1000
)

grp_pars_pig_01 <- list(
  grps = c('m0', 'm1', 'm2', 'm3', 'm4', 'm5', 'sr1'),
  yield = c(default = 0.05, sr1 = 0.065),
  xa_fresh = c(default = 0.06),
  xa_init = c(default = 0.06),
  d_max = c(default = 0.02),
  qhat_opt =  c(m0 = 0.45, m1 = 0.71, m2 = 1.0, m3 = 2.7, m4 = 4.4, m5 = 7.1, sr1 = 8.0),
  T_opt = c(m0 = 18, m1 = 18,   m2 = 28,   m3 = 36, m4 = 44, m5 = 55, sr1 = 44),
  T_min = c(m0 =  0, m1 =  9.5, m2 =  9.5, m3 = 15, m4 = 26, m5 = 30, sr1 =  0),
  T_max = c(m0 = 25, m1 = 25,   m2 = 38,   m3 = 45, m4 = 51, m5 = 60, sr1 = 51),
  ksmat = matrix(
    c(1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 0.46,
       NA,  NA,  NA,  NA,  NA,  NA, 0.0069),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c('VFA', 'SO4'),
      c('m0', 'm1', 'm2', 'm3', 'm4', 'm5', 'sr1')
    )
  ),
  mstoich = matrix(
    c(
      -1, -1, -1, -1, -1, -1, -1,
       0,  0,  0,  0,  0,  0, -0.50156,
       1,  1,  1,  1,  1,  1,  0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c('VFA', 'SO4', 'CH4'),
      c('m0', 'm1', 'm2', 'm3', 'm4', 'm5', 'sr1')
    )
  )
)

sub_pars_pig_01 <- list(
  subs = 'PS',
  sub_enrich = 'PS',
  xd = 'xd',
  h_rate_ref = c(PS = 0.04, xd = 0.032),
  h_rate_q10 = c(PS = 2.1, xd = 2.1),
  sub_fresh = c(PS = 63, xd = 0),
  sub_init  = c(PS = 63, xd = 0)
)

save(inf_pars_pig_01, file = '../data/inf_pars_pig_01.rda')
save(grp_pars_pig_01, file = '../data/grp_pars_pig_01.rda')
save(sub_pars_pig_01, file = '../data/sub_pars_pig_01.rda')

