# Creates schedule for both regular and time-variable scenarios

get_schedule <- function(
  structure, 
  times, 
  days, 
  delta_t
) {

  if (structure$type != 'ready') {
    stop('structure type is not yet \"ready\"--something is wrong.')
  }

  # Sort out times returned by ODE solver
  if (is.null(times)) {
    times <- seq(0, days, by = delta_t)
  }

  dat <- structure$dat

  # Include days argument in times vector
  times <- sort(unique(c(times, dat$time)))
    
  # Notes about time: 1) All simulations start at 0, 2) days must be at least as long as var data
  # Note that this works even with t_end = NULL (is ignored)
  # Note the "dummy" placeholder in position 1 (and extra + 1 in n_int)
  n_int <- nrow(dat)
  st <- dat[, 'time']
  timelist <- cumtime <- as.list(rep(0, n_int))
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
    timelist[[i]] <- tt
  }

  return(timelist)

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
fix_structure <- function(
  structure,
  pars,
  days
) {

  if (structure$type == 'ready') {
    return(structure)
  }

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

  # Trim unused times
  dat <- dat[dat$time <= days, ]

  # Check for sorted time
  if (is.unsorted(dat$time)) {
    stop('Column `time` must be sorted when time-variable parameters are used (pars), but it is not: ',
         head(dat$time))
  }
  
  # If simulation continues past pars data frame time, extend last row all the way
  if (dat[nrow(dat), 'time'] < days) {
    t_end <- days
    dat <- rbind(dat, dat[nrow(dat), ])
    dat[nrow(dat), 'time'] <- days
    # But make sure washing is not repeated!
    if (ncol(dat) > 2) {
      dat[nrow(dat), 3:ncol(dat)] <- 0
    }
  }

  # For 'mid' option, other variables are copied from previous time
  if (pars$approx_method == 'mid') {
    # Get midpoint time
    ir <- which(- c(0, diff(dat[, 'slurry_mass'])) > 0)
    tt <- (dat[ir, 'time']  + dat[ir - 1, 'time']) / 2
    nr <- dat[ir, ]
    nr$time <- tt
    dat <- rbind(dat, nr)
    dat <- dat[order(dat$time), ]
  } 

  structure$dat <- dat
  
  return(structure)

}

# Sorts out slurry production values and removal timing
calc_prod_rem <- function(
  structure,
  pars
) {
  
  if (structure$type == 'ready') {
    return(structure)
  }
  if (structure$type == 'regular') {
    stop('Yo! Ya gotta add code for making var out of regular!')
  }

  dat <- structure$dat

  # Removals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Note final 0--alignment is a bit tricky
  if (pars$approx_method %in% c('late', 'mid')) {
    removals <- - c(0, diff(dat[-nrow(dat), 'slurry_mass']), 0) > 0
  } else if (pars$approx_method == 'early') {
    removals <- - c(diff(dat[, 'slurry_mass']), 0) > 0
  } 
  dat$removal <- removals

  # Slurry production rate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  slurry_prod_rate_t <- c(diff(dat[, 'slurry_mass']) / diff(dat[, 'time']), 0) 
  slurry_prod_rate_t[slurry_prod_rate_t < 0] <- 0
  slurry_prod_rate_t[!is.finite(slurry_prod_rate_t)] <- 0
  dat$slurry_prod_rate <- slurry_prod_rate_t

  # Residual slurry for emptying ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (pars$approx_method == 'late') {
    dat$resid_mass <- dat$slurry_mass
  } else {
    dat$resid_mass <- c(dat$slurry_mass[-1], 0)
  }

  return(structure)
}



