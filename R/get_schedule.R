# Creates schedule for both regular and time-variable scenarios

get_schedule <- function(
  pars, 
  times, 
  days, 
  delta_t
) {

  # Figure out if this is regular or variable
  if (is.null(pars$var)) {
    isvar <- FALSE
  } else if (inherits(pars$var, 'data.frame')) {
    isvar <- TRUE
  } else {
    stop('pars_var must be NULL or a data frame')
  }
  
  if (isvar) {

    # Sort out times returned by ODE solver
    if (is.null(times)) {
      times <- seq(0, days, by = delta_t)
    }

    # Include days argument in times vector
    times <- sort(unique(c(times, pars$var$time)))
      
    # Notes about time: 1) All simulations start at 0, 2) days must be at least as long as var data
    # Note that this works even with t_end = NULL (is ignored)
    # Note the "dummy" placeholder in position 1 (and extra + 1 in n_int)
    n_int <- nrow(pars$var)
    st <- pars$var[, 'time']
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

  } else {

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
