# Schedule and related timing functions

get_schedule <- function(
  series, 
  times, 
  days, 
  delta_t
) {

  # Trim 
  series <- series[series$time <= days,]

  # Sort out times returned by ODE solver
  if (is.null(times)) {
    times <- seq(0, days, by = delta_t)
  }

  # Include days argument in times vector
  times <- sort(unique(c(times, series$time)))
    
  # Notes about time: 1) All simulations start at 0, 2) days must be at least as long as var series
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


# Checks and prepares slurry mass series
# Applies approx_method to slurry_mass
extract_series <- function(
  storage,
  pars,
  days
) {

  # If inputs are for regular schedule, create var data frame
  if (storage$type == 'regular') {

    # If empty interval is set to 0 or NA the storage is never emptied. 
    # Set to later than complete simulation
    empty_int <- storage$empty_int
    if(empty_int == 0 || is.na(empty_int)) {
      empty_int <- days + 1
    }
    
    # Figure out time intervals for loop
    if (!is.na(storage$wash_int) && storage$wash_water > 0) {  
      wash_int <- storage$wash_int
      rest_d <- storage$rest_d
    } else {
      wash_int <- Inf
      rest_d <- 0
    }

    # Continue sorting out intervals
    i <- 0
    t_int <- 0
    t_nowash <- 0
    wash <- FALSE

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

    # Create dat data frame
    # First resid mass, where first row is end of first interval
    resid_mass <- c(storage$slurry_mass, rep(storage$resid_depth * pars$area * pars$dens, length(t_int) - 1))
    # dat has an additional (time 0) row
    dat <- data.frame(
      time = cumsum(c(0, t_int)), 
      slurry_mass = c(
        storage$slurry_mass,
	storage$slurry_prod_rate * t_int + resid_mass
      ),
      resid_mass = storage$resid_depth * pars$area * pars$dens,
      removal = TRUE,
      wash_water = storage$wash_water,
      rest_days = storage$rest_d,
      slurry_prod_rate = storage$slurry_prod_rate
    )

    # No intial (time 0) and no final removal
    dat[c(1, nrow(dat)), 'removal'] <- FALSE

    return(dat)

  } else {

    # Extract dat
    dat <- storage$dat

    # Add any missing time 0
    if (dat[1, 'time'] > 0) {
      dat <- rbind(c(0, dat$slurry_mass[1]), dat)
    }

    # For 'mid' option, other variables are copied from previous time
    if (pars$approx_method == 'mid') {
      # Get midpoint time
      ir <- which(- c(0, diff(series[, 'slurry_mass'])) > 0)
      tt <- (series[ir, 'time']  + series[ir - 1, 'time']) / 2
      nr <- series[ir, ]
      nr$time <- tt
      dat <- rbind(dat, nr)
      dat <- dat[order(dat$time), ]
    } 

    # Trim unused times
    dat <- dat[dat$time <= days, ]

    # Calculate slurry production rate and identify removal events
    dat <- calc_prod_rem(dat, pars)

    return(dat)

  }
}

# Additional series processing
clean_series <- function(
  series,
  pars,
  days
){

  ## Check for the right columns
  #if (ncol(series) != 2 || !identical(names(series), c('time', 'slurry_mass'))) {
  #  stop('The storage series element must have two columns: time and slurry_mass.')
  #}

  # Cannot have no slurry present because is used in all concentration calculations
  series[series[, 'slurry_mass'] == 0, 'slurry_mass'] <- 1E-10

  # Add in var_pars seriesa frame if present
  if (inherits(pars$var, 'data.frame')) {
    series <- merge(series, pars$var, by = 'time', all = TRUE)

    # Fill in missing values after merge (if series and var_pars have different times)
    # But added removal should = FALSE
    series[is.na(series$removal), 'removal'] <- FALSE

    # And slurry_prod_rate should always use forward
    series$slurry_prod_rate <- fill_down(series$slurry_prod_rate)

    if (pars$fill_method == 'interp') {
      if (any(sapply(series, class) == 'list')) {
        warning('ctrl_pars element fill_method \"interp\" cannot be used with list elements, so reverting to \"forward\"')
        pars$fill_method <- 'forward'
      } else {
        series <- interpm(series, 'time', names(series)[-1], rule = 2)
      }
    } 

    if (pars$fill_method == 'forward') {
       series <- fill_down_df(series)
    }
  } 

  # Make sure series is sorted by time
  series <- series[order(series$time), ]

  # Check for duplicated times or other problems
  if (any(duplicated(series$time))) {
    stop('Duplicated times in series object. Check var_pars and storage inputs.')
  }
  
  # If simulation continues past pars series time, extend last row all the way
  if (series[nrow(series), 'time'] < days) {
    t_end <- days
    series <- rbind(series, series[nrow(series), ])
    series[nrow(series), 'time'] <- days
    # But may need to make sure washing is not repeated somehow!
  }

  return(series)

}

# Sorts out slurry production values and removal timing for mass series input
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
  series$resid_mass[!series$removal] <- NA

  return(series)
}


