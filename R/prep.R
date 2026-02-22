# Creates schedule for both regular and time-variable scenarios


# Create initial state variable vector
get_init_state <- function(
  structure, 
  pars
) {
 
  if (structure$type == 'series') {
    slurry_mass_init <- structure$dat[1, 'slurry_mass']
  } else {
    slurry_mass_init <- structure$slurry_mass
  }

  if (!is.null(pars$kl)) {
    emis <- pars$kl * 0
    names(emis) <- paste0(names(emis), '_emis_cum')
  } else {
    emis <- NULL
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


