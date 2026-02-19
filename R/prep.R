# Creates schedule for both regular and time-variable scenarios


# Create initial state variable vector
get_init_state <- function(
  structure, 
  pars
) {
 
  slurry_mass_init <- structure$dat[1, 'slurry_mass']

  if (!is.null(pars$kl)) {
    emis <- pars$kl * 0
    names(emis) <- paste0(names(emis), '_emis_cum')
  } else {
    emis <- NULL
  }
  
  y <- c(
    pars$xa_init * slurry_mass_init,                       # Multiple microbial groups
    pars$sub_init[pars$subs] * slurry_mass_init,           # Multiple particulate substrates 
    pars$conc_init * slurry_mass_init,                     # VFA and conservative solutes
    slurry_mass = slurry_mass_init, 
    CH4_emis_cum = 0, 
    emis,
    slurry_load = 0,
    COD_load = 0
  )
  
  return(y)

}


