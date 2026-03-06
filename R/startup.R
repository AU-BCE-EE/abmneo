# Internal function for repeated abm() calls
# Argument list should exactly match abm()

abm_startup <- function(
  storage,
  days,
  delta_t,
  times,
  pars,
  startup,
  starting,
  warn
) {

  for (i in 1:(startup + 1)) {
    if (i > startup) {
      cat('and final run')
      cat('\n')
    } else {
      cat(paste0(i, 'x -> '))
    }

    if (i > 1) {
      starting <- out
    } 
    
    # Call abm() with arguments given in outside call except for startup
    out <- abmneo(
      storage = storage,
      days = days,
      delta_t = delta_t,
      times = times,
      pars = pars,
      startup = 0,
      starting = starting,
      warn = warn
    )

  }

  return(out)

}



# Move pars from starting argument into place
starting_pars <- function(pars, starting) {
  
  message('Using starting conditions from `starting` argument')
  pars$conc_init[pars$grps] <- as.numeric(starting[nrow(starting), paste0(pars$grps, '_conc')])
  # Remember xa_init and other stuff is grouped into conc_init, so next line not needed
  # Potential source of confusion!
  #pars$xa_init[pars$grps] <- as.numeric(starting[nrow(starting), paste0(pars$grps, '_conc')])
  pars$conc_init['VFA'] <- as.numeric(starting[nrow(starting), 'VFA_conc'])
  pars$sub_init[pars$subs] <- as.numeric(starting[nrow(starting), paste0(pars$subs, '_conc')])
  # Set slurry_mass as well?
  # NTS: Set comp solutes also? Yes!

  return(pars)
  
}


