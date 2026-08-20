# Set fixed inputs

stors <- list()

ids <- c('AD1', 'AD2', 'AD3', 'AD4')

for (i in ids) {
  stors[[i]] <- list(
    type = 'series',
    dat = as.data.frame(mass[tank == i, -3]),
    storage_depth = tanks[tank == i, depth],     
    area = tanks[tank == i, area],              
    temp_C = 0,
    resid_enrich = 0.2 
  )
}

infls <- list()

for (i in ids) {
  infls[[i]] <- list(
    VFA_fresh = tanks[tank == i, VFA],
    VFA_init = tanks[tank == i, VFA],
    pH = tanks[tank == i, pH],
    dens = 1000
  )
}

subs <- list()

for (i in ids) {
  subs[[i]] <- list(
    subs = 'PS',
    sub_enrich = 'PS',
    xd = 'xd',
    h_rate_ref = c(PS = 0.1, xd = 0.1),
    h_rate_q10 = c(PS = 2, xd = 2),
    sub_fresh = c(PS = tanks[tank == i, substrate], xd = 0),
    sub_init  = c(PS = tanks[tank == i, substrate], xd = 0)
  )
}

temps <- list()
for (i in ids) {
  dd <- temp[tank == i, .(doy, slurry_temp_interp)]
  names(dd) <- c('time', 'temp_C')
  temps[[i]] <- as.data.frame(dd)
}

times <- list()
for (i in ids) {
  times[[i]] <- sort(unique(c(0, unname(unlist(emis[tank == i, .(doy_start, doy_end)])))))
}


