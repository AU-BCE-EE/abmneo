# Creates schedule for both regular and time-variable scenarios

get_schedule <- function(
  series, 
  times, 
  days, 
  delta_t
) {

  # Sort out times returned by ODE solver
  if (is.null(times)) {
    times <- seq(0, days, by = delta_t)
  }

  # Include days argument in times vector
  times <- sort(unique(c(times, series$time)))
    
  # Notes about time: 1) All simulations start at 0, 2) days must be at least as long as var seriesa
  # Note that this works even with t_end = NULL (is ignored)
  # Note the "dummy" placeholder in position 1 (and extra + 1 in n_int)
  n_int <- nrow(series)
  st <- series[, 'time']
  schedule <- cumtime <- as.list(rep(0, n_int))
  for (i in 2:n_int) {
    tt <- times[times > st[i - 1] & times <= st[i]] # Earlier note suggested this line was slow
    # Simulation intervals should end exactly at emptying time, so it is added here through st[i] for cases where there is not alignment
    tt <- unique(c(tt, st[i]))
    if (length(tt) == 0) { 
      # Not clear when this might happen, but it would be bad
      stop('No times for interval (row) ', i, ' in time list for some reason. xyn917')
    } else {
      cumtime[[i]] <- max(tt)
      tt <- tt - cumtime[[i - 1]]
      tt <- unique(c(0, tt))
    }
    schedule[[i]] <- tt
  }

  return(schedule)

}


getRegTimes <- function(pars, days) {
  
  # If empty interval is set to 0 or NA the storage is never emptied. 
  empty_int <- pars$empty_int
  if(empty_int == 0 || is.na(empty_int)) {
    empty_int <- days + 1
  }
  
  # Figure out time intervals for loop
  if (!is.na(pars$wash_int) && pars$wash_water > 0) {  
    wash_int <- pars$wash_int
    rest_d <- pars$rest_d
  } else {
    wash_int <- Inf
    rest_d <- 0
  }
  wash_rest_int <- wash_int + rest_d

  # Continue sorting out intervals
  i <- 0
  t_int <- 0
  t_nowash <- 0
  wash <- FALSE

  # Continute . . .
  # Each interval is either 1) the fixed empty_int or if time between washings would be exceeded, 
  # 2) time to get to a washing event, or 3) time until end of simulation
  while (sum(t_int, wash * rest_d) < days) {
    i <- i + 1
    t_int[i] <- min(wash_int - t_nowash, empty_int, days - sum(t_int, wash * rest_d))
    if (t_int[i] == wash_int - t_nowash) {
      wash[i] <- TRUE
      t_nowash <- 0
    } else {
      wash[i] <- FALSE
      t_nowash <- t_nowash + t_int[i]
    }
  }

  # Number of empty or wash intervals
  n_int <- length(t_int)

  return(list(n_int = n_int,
              empty_int = empty_int,
              wash_int = wash_int,
              rest_d = rest_d,
              t_int = t_int,
              wash = wash))

}


# Checks and prepares slurry mass series
# Applies approx_method to slurry_mass
make_series <- function(
  structure,
  pars,
  days
) {

  # If inputs are for regular schedule, create var data frame
  if (structure$type == 'regular') {
    stop('Yo! Ya gotta add code for making var out of regular!')
  }

  # Extract dat
  dat <- structure$dat

  # Add missing time 0
  if (dat[1, 'time'] > 0) {
    dat <- rbind(c(0, dat$slurry_mass[1]), dat)
  }

  # Check for the right columns
  if (ncol(dat) != 2 || !identical(names(dat), c('time', 'slurry_mass'))) {
    stop('The structure dat element must have two columns: time and slurry_mass.')
  }

  # Cannot have no slurry present because is used in all concentration calculations
  dat[dat[, 'slurry_mass'] == 0, 'slurry_mass'] <- 1E-10

  # Add in var_pars data frame (if NULL, then no effect)
  series <- merge(dat, pars$var, by = 'time', all = TRUE)

  # Fill in missing values
  if (pars$fill_method == 'interp') {
    series <- interpm(series, 'time', names(series)[-1], rule = 2)
  } else {
    stop('ctrl_pars element fill_method--only available option is \"interp\".')
  }

  series <- calc_prod_rem(series, pars)

  # Trim unused times
  series <- series[series$time <= days, ]

  # Make sure series is sorted by time
  series <- series[order(series$time), ]

  # Check for duplicated times or other problems
  if (any(duplicated(series$time))) {
    stop('Duplicated times in series object. Check var_pars and structure inputs.')
  }
  
  # If simulation continues past pars seriesa frame time, extend last row all the way
  if (series[nrow(series), 'time'] < days) {
    t_end <- days
    series <- rbind(series, series[nrow(series), ])
    series[nrow(series), 'time'] <- days
    # But make sure washing is not repeated!
    if (ncol(series) > 2) {
      series[nrow(series), 3:ncol(series)] <- 0
    }
  }

  # For 'mid' option, other variables are copied from previous time
  if (pars$approx_method == 'mid') {
    # Get midpoint time
    ir <- which(- c(0, diff(series[, 'slurry_mass'])) > 0)
    tt <- (series[ir, 'time']  + series[ir - 1, 'time']) / 2
    nr <- series[ir, ]
    nr$time <- tt
    series <- rbind(series, nr)
    series <- series[order(series$time), ]
  } 

  return(series)

}

# Sorts out slurry production values and removal timing
calc_prod_rem <- function(
  series,
  pars
) {
  
  # Removals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Note final 0--alignment is a bit tricky
  if (pars$approx_method %in% c('late', 'mid')) {
    removals <- - c(0, diff(series[-nrow(series), 'slurry_mass']), 0) > 0
  } else if (pars$approx_method == 'early') {
    removals <- - c(diff(series[, 'slurry_mass']), 0) > 0
  } 
  series$removal <- removals

  # Slurry production rate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  slurry_prod_rate_t <- c(diff(series[, 'slurry_mass']) / diff(series[, 'time']), 0) 
  slurry_prod_rate_t[slurry_prod_rate_t < 0] <- 0
  slurry_prod_rate_t[!is.finite(slurry_prod_rate_t)] <- 0
  series$slurry_prod_rate <- slurry_prod_rate_t

  # Residual slurry for emptying ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (pars$approx_method == 'late') {
    series$resid_mass <- series$slurry_mass
  } else {
    series$resid_mass <- c(series$slurry_mass[-1], 0)
  }

  return(series)
}



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


getRegTimes <- function(pars, days) {
  
  # If empty interval is set to 0 or NA the storage is never emptied. 
  empty_int <- pars$empty_int
  if(empty_int == 0 || is.na(empty_int)) {
    empty_int <- days + 1
  }
  
  # Figure out time intervals for loop
  if (!is.na(pars$wash_int) && pars$wash_water > 0) {  
    wash_int <- pars$wash_int
    rest_d <- pars$rest_d
  } else {
    wash_int <- Inf
    rest_d <- 0
  }
  wash_rest_int <- wash_int + rest_d

  # Continue sorting out intervals
  i <- 0
  t_int <- 0
  t_nowash <- 0
  wash <- FALSE

  # Continute . . .
  # Each interval is either 1) the fixed empty_int or if time between washings would be exceeded, 
  # 2) time to get to a washing event, or 3) time until end of simulation
  while (sum(t_int, wash * rest_d) < days) {
    i <- i + 1
    t_int[i] <- min(wash_int - t_nowash, empty_int, days - sum(t_int, wash * rest_d))
    if (t_int[i] == wash_int - t_nowash) {
      wash[i] <- TRUE
      t_nowash <- 0
    } else {
      wash[i] <- FALSE
      t_nowash <- t_nowash + t_int[i]
    }
  }

  # Number of empty or wash intervals
  n_int <- length(t_int)

  return(
    list(
      n_int = n_int,
      empty_int = empty_int,
      wash_int = wash_int,
      rest_d = rest_d,
      t_int = t_int,
      wash = wash
    )
  )

}
