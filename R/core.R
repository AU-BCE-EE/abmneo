
# Exported abm() function
abm <- function(
  structure,
  days = 365,
  delta_t = 1,
  times = NULL,
  mng_pars = NULL,
  man_pars = NULL,
  init_pars = list(conc_init =  c(man_pars$comp_fresh, man_pars$VFA_fresh)),
  grp_pars = NULL,
  sub_pars = NULL,
  chem_pars = NULL,
  inhib_pars = NULL,
  mt_pars = NULL,
  ctrl_pars = list(
    respir = TRUE,
    pH_inhib = FALSE, 
    approx_method = 'early',
    fill_method = 'interp',
    par_key = '\\.',
    rates_calc = 'instant'
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
      mng_pars = mng_pars,
      man_pars = man_pars,
      init_pars = init_pars,
      grp_pars = grp_pars,
      sub_pars = sub_pars,
      chem_pars = chem_pars,
      inhib_pars = inhib_pars,
      ctrl_pars = ctrl_pars,
      var_pars = var_pars,
      add_pars = add_pars,
      days = days
    )
  }
  
  if (!is.null(starting) & inherits(starting, 'data.frame')) {
    pars <- starting_pars(pars, starting)
  }

  # If startup repetitions are requested, repeat some number of times before returning results
  # Uses pars, already packed and with starting values added
  if (startup > 0) {
    cat('\nStartup run ')
    out <- abm_startup(
      structure = structure,
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

  # Sort out structure, extending time as needed, checking for required components, and adding in any var_pars
  series <- extract_series(structure, pars, days)
  series <- clean_series(series, pars, days)

  # Create initial state variable vector
  y <- get_init_state(structure, pars) 

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
    times,
    pars, 
    addcols = TRUE, 
    addconcs = TRUE, 
    cumeff = TRUE
  )

  # Check COD balance
  codbal <- check_COD(
    dat = dat, 
    grps = pars$grps, 
    subs = pars$subs, 
    COD_conv = pars$COD_conv, 
    stoich = pars$stoich, 
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

    # Calculate log Ka for speciation
    pars <- calc_ka(pars)
    
    # Create default y.eff vector with zeros because washing could occur, and dat needs columns
    y.eff <- 0 * empty_store(y)$eff

    # If there is a removal event, remove slurry before calling up ODE solver
    if (pars$removal) {
      y <- empty_store(y, resid_mass = pars$resid_mass, resid_enrich = pars$resid_enrich)
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
