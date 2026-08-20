

res <- optim(
  par = c(h_rate_ref.PS = -1, qhat_opt.m1 = -1, qhat_opt.m2 = -1),
  fn = residuals,
  stors = stors,
  times = times, 
  infls = infls, 
  grp_pars = grp_pars, 
  subs = subs, 
  temps = temps
)

args(optim)


residuals(stors, times, infls, grp_pars, subs, temps, add_pars = NULL)
residuals(stors, times, infls, grp_pars, subs, temps, add_pars = list(h_rate_ref.PS = 0.1))
residuals(stors, times, infls, grp_pars, subs, temps, add_pars = list(h_rate_ref.PS = 0.0))
residuals(stors, times, infls, grp_pars, subs, temps, add_pars = list(xa_fresh.m1 = 0.5))
Q

residuals(stors, times, infls, grp_pars, subs, temps, add_pars = list(h_rate_ref = c(PS = 0.5, xd = 0)))

out <- abmneo(
  storage = stors[[1]],
  days = 365,
  times = times[[1]],
  inf_pars = infls[[1]],
  grp_pars = grp_pars,
  sub_pars = subs[[1]],
  var_pars = list(var = temps[[1]]),
  add_pars = NULL,
  startup = 1
)

# Does not work
out <- abmneo(
  storage = stors[[1]],
  days = 365,
  times = times[[1]],
  inf_pars = infls[[1]],
  grp_pars = grp_pars,
  sub_pars = subs[[1]],
  var_pars = list(var = temps[[1]]),
  add_pars = list(h_rate_ref.PS = 0.1),
  startup = 1
)

# Internally, I see:
# Browse[2]> pars$h_rate_ref
# $PS
# [1] 0.5
# 
# $<NA>
# NULL

# Works as expected
out <- abmneo(
  storage = stors[[1]],
  days = 365,
  times = times[[1]],
  inf_pars = infls[[1]],
  grp_pars = grp_pars,
  sub_pars = subs[[1]],
  var_pars = list(var = temps[[1]]),
  add_pars = list(h_rate_ref = c(PS = 0.1, xd = 0.01)),
  startup = 1
)

# Internally, I see:
# Browse[1]> pars$h_rate_ref
#  PS  xd 
# 0.1 0.1 
  
p  
pars$xa_fresh
pars$h_rate_ref
Q
c
pars
optim()
