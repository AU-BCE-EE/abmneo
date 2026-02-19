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

# Function for emptying storage struction
# Note: enrich_names could also be length 1 vector: '^xa|^RFd$|^iNDF$' etc.
empty_store <- function(
  y, 
  resid_mass = 0, 
  resid_enrich = 0, 
  skip = FALSE,
  enrich_names = c('^xa', '^m[0-9]', '^sr[0-9]', '^RFd$', '^CPs$', '^iNDF$', '^VSd$', '^Cfat$', '^starch$', '^ash$'),
  ignore_names = c('_emis', '_load', '_cum_', '_conv_', 'cum$'),
  warn = TRUE
) {

  y <- unlist(y)
  slurry_mass <- y['slurry_mass']

  which.ignore <- grepl(paste(ignore_names, collapse = '|'), names(y))

  if (slurry_mass > resid_mass & !skip) {
    # Masses before emptying
    y.before <- y
    # Calculate mass of each variable remaining in storage
    resid_frac <- resid_mass / slurry_mass
    resid_par <- logistic(logit(resid_frac) + resid_enrich)
    which.enrich <- grepl(paste(enrich_names, collapse = '|'), names(y))
    y[which.enrich] <- y[which.enrich] * resid_par
    y[!which.enrich & !which.ignore] <- y[!which.enrich & !which.ignore] * resid_frac
    # Effluent
    y.eff <- y.before - y
    y.eff <- y.eff[!which.ignore]
    names(y.eff) <- paste0(names(y.eff), '_eff')
  } else {
    if (warn) {
      warning('Emptying skipped.')
    }
    y.eff <- 0 * y
    y.eff <- y.eff[!which.ignore]
    names(y.eff) <- paste0(names(y.eff), '_eff')
  }

  return(list(store = y, eff = y.eff))

}


get_last_state <- function(out, y) {
  y <- unlist(out[nrow(out), names(y)])
}
