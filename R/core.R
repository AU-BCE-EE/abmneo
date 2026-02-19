# Heart of abm(), with lsoda call

abm_core <- function(
  days,
  schedule,
  y, 
  series,
  pars, 
  warn
) {
    
  n_int <- length(schedule)
  
  # Empty data frame for holding results
  dat <- NULL

  # Time trackers
  # Time remaining to run
  t_rem <- days
  # Time that has already run
  t_run <- 0

  # Start the time (emptying) loop
  for (i in 2:n_int) {

    # Sort out call duration
    t_call <- min(max(schedule[[i]]), t_rem)

    # Fill in current pars from var
    # Also adds 2 temperature-dependent derivative vectors
    pars <- update_var_pars(series, pars, y, i - 1)

    # Calculate log Ka for speciation
    pars <- calcKa(pars)
    
    # Create default y.eff vector with zeros because washing could occur, and dat needs columns
    y.eff <- 0 * emptyStore(y)$eff

    # If there is a removal event, remove slurry before calling up ODE solver
    if (pars$removal) {
      y <- emptyStore(y, resid_mass = pars$resid_mass, resid_enrich = pars$resid_enrich)
      y.eff <- y$eff
      y <- y$store
    }

    # Get times for lsoda() call
    # Need some care with times to make sure t_call is last one in case it is not multiple of delta_t
    tt <- schedule[[i]]

    # Call up ODE solver
    out <- deSolve::lsoda(y = y, 
                          times = tt, 
                          rates, 
                          parms = pars)
     
    # Change format of output and drop first (time 0) row (duplicated in last row of previous)
    if (i == 2) {
      out <- data.frame(out)
    } else {
      out <- data.frame(out[-1, , drop = FALSE])
    }

    # Extract new state variable vector from last row of lsoda output
    y <- getLastState(out, y)

    # Add effluent results
    out <- add_eff(out, y.eff)
 
    # Clean up and stack output with earlier results
    dat <- add_out(dat, out)
    
    # Update time remaining and total time run so far
    t_rem <- t_rem - t_call
    t_run <- t_run + t_call
    
  }

  return(dat)
}
