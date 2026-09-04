# Creates schedule for both regular and time-variable scenarios


# Create initial state variable vector
get_init_state <- function(
  pars,
  series
) {
 
  slurry_mass_init <- series[1, 'slurry_mass']

  gas_init <- rep(0, length(pars$gases))
  names(gas_init) <- pars$gases

  y <- c(
    pars$conc_init * slurry_mass_init,                     # Microbial groups, (multiple) particulate substrates, VFA and conservative solutes
    slurry_mass = slurry_mass_init, 
    gas_init, 
    slurry_load = 0,
    COD_load = 0
  )

  return(y)

}


