# Creates schedule for both regular and time-variable scenarios


# Create initial state variable vector
get_init_state <- function(
  storage, 
  pars
) {
 
  if (storage$type == 'series') {
    slurry_mass_init <- storage$dat[1, 'slurry_mass']
  } else {
    slurry_mass_init <- storage$slurry_mass
  }

  y <- c(
    pars$conc_init * slurry_mass_init,                     # Microbial groups, (multiple) particulate substrates, VFA and conservative solutes
    slurry_mass = slurry_mass_init, 
    CH4 = 0, 
    slurry_load = 0,
    COD_load = 0
  )

  return(y)

}


