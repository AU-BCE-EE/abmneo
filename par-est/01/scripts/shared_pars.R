# Shared pars

10^seq(log10(.5), log10(8), length.out = 5) 
grp_pars <- list(
  grps = c('m1', 'm2','m3','m4','m5'),
  yield = c(default = 0.05),
  xa_fresh = c(default = 0.06),
  xa_init = c(default = 0.06),
  d_max = c(default = 0.02),
  qhat_opt =  c(m1 = 0.5, m2 = 1, m3 = 2, m4 = 4, m5 = 8),
  T_opt = c(m1 = 18, m2 = 28, m3 = 36, m4 = 44, m5 = 55),
  T_min = c(m1 =  2, m2 =  2, m3 = 15, m4 = 26, m5 = 30),
  T_max = c(m1 = 25, m2 = 38, m3 = 45, m4 = 51, m5 = 60),
  ksmat = matrix(
    c(1.1, 1.1, 1.1, 1.1, 1.1),
    nrow = 1,
    byrow = TRUE,
    dimnames = list(
      c('VFA'),
      c('m1', 'm2', 'm3', 'm4', 'm5')
    )
  ),
  mstoich = matrix(
    c(
      -1, -1, -1, -1, -1,
       1,  1,  1,  1,  1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c('VFA', 'CH4'),
      c('m1', 'm2', 'm3', 'm4', 'm5')
    )
  )
)


