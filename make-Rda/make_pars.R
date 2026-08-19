# Creates parameter objects

grp_VS_dig3 <- list(
  grps = c('m0', 'm1', 'm2','m3','m4','m5'),
  yield = c(default = 0.05),
  xa_fresh = c(default = 0.06),
  xa_init = c(default = 0.06),
  d_max = c(default = 0.02),
  qhat_opt =  c(m0 = 1.025, m1 = 1.642, m2 = 2.839, m3 = 6.154, m4 = 10.26, m5 = 16.4),
  T_opt = c(m0 = 18, m1 = 18,    m2 = 28,    m3 = 36, m4 = 43.75, m5 = 55),
  T_min = c(m0 =  0, m1 =  2.44, m2 =  2.44, m3 = 15, m4 = 26.25, m5 = 30),
  T_max = c(m0 = 25, m1 = 25,    m2 = 38,    m3 = 45, m4 = 51.25, m5 = 60),
  ksmat = matrix(
    c(1.153, 1.153, 1.153, 1.153, 1.153, 1.153),
    nrow = 1,
    byrow = TRUE,
    dimnames = list(
      c('VFA'),
      c('m0', 'm1', 'm2', 'm3', 'm4', 'm5')
    )
  ),
  mstoich = matrix(
    c(
      -1, -1, -1, -1, -1, -1,
       1,  1,  1,  1,  1,  1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c('VFA', 'CH4'),
      c('m0', 'm1', 'm2', 'm3', 'm4', 'm5')
    )
  )
)

sub_VS_dig3 <- list(
  subs = c('PS'),
  sub_enrich = c('PS'),
  xd = 'xd',
  arrA = c(default = 2.250e7),
  arrE = c(default = 55609),
  sub_fresh = c(PS = 18.22, xd = 0),
  sub_init  = c(PS = 18.22, xd = 0)
)

inf_VS_dig3 <- list(
  VFA_fresh = 0.5,
  VFA_init = 0.t,
  pH = 8,
  dens = 1000
)


