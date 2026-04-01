# Exported abm() function
abmneo <- function(
  storage,
  days = 365,
  delta_t = 1,
  times = NULL,
  inf_pars = NULL,
  grp_pars = NULL,
  sub_pars = NULL,
  chem_pars = list(COD_conv = c(CH4 = 5.32), gases = c('CH4', 'CO2', 'H2S')),
  ctrl_pars = list(
    approx_method = 'early',
    fill_method = 'interp',
    par_key = '\\.',
    arr_max_temp_K = 313
  ),
  var_pars = list(var = NULL),
  add_pars = NULL,
  pars = NULL,
  startup = 0,
  starting = NULL,
  warn = TRUE
) {

  # Sort out parameters, package all parameters into a single pars list, add some others
  # Includes sorting out var_pars
  # All these steps are skipped if pars is provided
  if (is.null(pars)) {
    pars <- pack_pars(
      storage = storage,
      inf_pars = inf_pars,
      grp_pars = grp_pars,
      sub_pars = sub_pars,
      chem_pars = chem_pars,
      ctrl_pars = ctrl_pars,
      var_pars = var_pars,
      add_pars = add_pars,
      days = days
    )
  }
  
  if (!is.null(starting) & inherits(starting, 'data.frame')) {
    pars$conc_init[pars$supercomps] <- as.numeric(starting[nrow(starting), paste0(pars$supercomps, '_conc')])
  }

  # If startup repetitions are requested, repeat some number of times before returning results
  # Uses pars, already packed and with starting values added
  if (startup > 0) {
    cat('\nStartup run ')
    out <- abm_startup(
      storage = storage,
      days = days,
      delta_t = delta_t,
      times = times,
      pars = pars,
      startup = startup,
      starting = starting,
      warn = warn
    )
    return(out)
  } 

  # Sort out storage, extending time as needed, checking for required components, and adding in any var_pars
  series <- extract_series(storage, pars, days)
  series <- clean_series(series, pars, days)

  # Create initial state variable vector
  y <- get_init_state(storage, pars) 

  # Get timing of intervals (list of times)
  schedule <- get_schedule(
    series, 
    times, 
    days, 
    delta_t
  )

  dat <- abm_core(
    days = days,
    schedule = schedule, 
    y = y, 
    series = series,
    pars = pars, 
    warn = warn
  )

  # Clean up and extend output
  dat <- clean_output(
    dat, 
    days,
    times,
    pars, 
    addcols = TRUE, 
    addconcs = TRUE, 
    cumeff = TRUE
  )

  # Check COD balance
  codbal <- check_COD(
    dat = dat, 
    pars = pars,
    rtol = 0.01
  )

  # Return results
  return(dat)

}

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

    # Create default y.eff vector with zeros because washing could occur, and dat needs columns
    y.eff <- 0 * empty_store(y, skip = TRUE, ignore_names = unique(c(pars$gases, '_emis', '_load', '_cum_', '_conv_', 'cum$')), warn = FALSE)$eff

    # If there is a removal event, remove slurry before calling up ODE solver
    # Gases are not removed (their state is cumulative emission)
    if (pars$removal) {
      y <- empty_store(
        y, 
	    resid_mass = pars$resid_mass, 
	    resid_enrich = pars$resid_enrich,
	    enrich_names = c(pars$grps, pars$sub_enrich),
	    ignore_names =  unique(c(pars$gases, '_emis', '_load', '_cum_', '_conv_', 'cum$'))
      )
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
    y <- get_last_state(out, y)

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
