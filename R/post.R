# Clean up abm_*() output before returning it

clean_output <- function(
  dat, 
  times, 
  pars, 
  addcols, 
  addconcs, 
  cumeff
) {

  # Drop times that were not requested in output
  if (!is.null(times)) {
    dat <- dat[dat$time %in% c(times, days), ]
  }

  # Replace . in names with _
  names(dat) <- gsub('\\.', '_', names(dat))
  
  if (addcols) {
    # Add slurry depth
    dat$slurry_depth <- dat$slurry_mass / pars$area / pars$dens
  }

  # Make effluent results cumulative
  if (cumeff) {
    dat[, names(dat)[grepl('_eff', names(dat))]] <- lapply(dat[, names(dat)[grepl('_eff', names(dat))]], cumsum)
  }

  # Add concentrations (after cumulative effluent)
  if (addconcs) {
    dat <- add_concs(dat, pars)
  }
  
  return(dat)

}

# Add concetrations to abm() output
add_concs <- function(dat, pars) {
  
    dat[, paste0(c(pars$grps, pars$subs, pars$sols), '_conc')] <-         dat[,        c(pars$grps, pars$subs, pars$sols)] / dat$slurry_mass
    dat[, paste0(c(pars$grps, pars$subs, pars$sols), '_eff', '_conc')] <- dat[, paste0(c(pars$grps, pars$subs, pars$sols), '_eff')] / dat$slurry_mass_eff

    return(dat)

}

# Add effluent results
add_eff <- function(dat, y.eff) {

  dat <- data.frame(dat)

  # Create new effluent columns only if they do not exist
  # Should never be a case when some but not all exist
  if (!all(names(y.eff) %in% names(dat))) {
    dat[, names(y.eff)] <- 0
  }
  
  # Add in effluent results
  dat[nrow(dat), names(y.eff)] <- y.eff

  return(dat)

}

# Add interval output to earlier results
# new can be lsoda output or state variable vector (after emptying)
add_out <- function(main, new = NULL) {
  
  # Change output from matrix to data frame
  # Do not drop first (time 0) row
  if (inherits(new, 'matrix')) {
    new <- data.frame(new)
  } 
  
  if (inherits(main, 'matrix')) {
    main <- data.frame(main)
  } 

  if (inherits(new, 'data.frame')) {
    # Get previous time from main, if it exists
    if ('time' %in% names(main)) {
      t_add <- max(main$time)
    } else {
      t_add <- 0
    }
    # Change time in output to cumulative time for complete simulation
    new$time <- new$time + t_add
    # Add results to earlier ones
    main <- rbind(main, new)
  } else if (inherits(new, 'numeric')){
    # Add row to main, duplicate last row
    main <- main[c(1:nrow(main), nrow(main)), ]
    # Replace values that are present in new
    main[nrow(main), names(new)] <-new 
  } else {
    stop('Class of new is not data frame or vector but is ', class(new))
  }
  
  return(main)

}
