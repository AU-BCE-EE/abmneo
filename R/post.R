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
