# Main abm function

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
      structure = structure,
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
  series <- make_series(structure, pars, days)

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
  codbal <- checkCOD(
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

