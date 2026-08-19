
# Sorts out parameters and packages them all together in the output list
# This is a central function that does a lot and is (unfortunately) complicated
pack_pars <- function(
  storage,
  inf_pars,
  grp_pars,
  sub_pars,
  chem_pars,
  ctrl_pars,
  var_pars,
  add_pars,
  days
) {

  # Time-variable pars initial processing
  if(!is.null(var_pars) && !is.null(var_pars$var)) {
    # Move extra var_pars into first (var) element, possibly as lists within each data frame element~
    var_pars <- combine_var_pars(var_pars)
    # Shift to mid time
    tt <- var_pars[['var']]$time
    var_pars[['var']]$time <- c(tt[1], tt[-length(tt)] + diff(tt) / 2)
  }

  # Combine pars to make extraction and pass to rates() easier ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # First extract out those pars that come through storage argument
  stor_pars <- storage[! names(storage) %in% c('type', 'dat')]
  pars <- c(stor_pars, inf_pars, grp_pars, sub_pars, chem_pars, ctrl_pars, var_pars)

  # Sort out add_pars and similar parameter inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pars <- fix_add_pars(pars, add_pars)

  # Check pars - NTS needs to be done!
  # * series dat must have only time and slurry_mass
  # * order of all grp pars must match
  # * check that all reactants in mstoich have a ksmat value, otherwise we get zero metabolism for that grp 
  # 

  # Check for identical dimensions in inhibition pars
  if (!is.null(pars$ic0) && !isTRUE(all.equal(dimnames(pars$ic0), dimnames(pars$ic100)))) {
    stop('ic0 and ic100 row and column names must be identical but are not.')
  }

  # Add xd to subs if missing
  if (! pars$xd %in% pars$subs) {
    pars$subs <- c(pars$subs, pars$xd)
  }

  # Fill in default values for grp_pars if keyword name `default` is used
  # Note: Microbial groups are defined by grps element
  # Note: `default` does *not* work with add_pars argument because grps are already defined in defaults
  # expandPars() will also sort out element order and drop excluded elements
  # NTS: Could these vectors of names be set in some kind of defaults?
  # Note that mstoich is *not* expanded! Too complicated. So it needs all the elements (rows/columns)
  grp_par_nms <- c("yield", "xa_fresh", "xa_init", "d_max", "qhat_opt", "T_opt", "T_min", "T_max")
  sub_par_nms <- c("T_opt_hyd", "T_min_hyd", "T_max_hyd", "hydrol_opt", "sub_fresh", "sub_init", "h_rate_ref", "h_rate_q10")
  pars <- expand_pars(pars = pars, elnms = pars$grps, parnms = grp_par_nms)
  pars <- expand_pars(pars = pars, elnms = pars$subs, parnms = sub_par_nms)

  # Do we need to check grp arguments, including order of element names in some pars?
  # After above block, this should be redundant
  # Otherwise see old check_grp_names() code

  # Size-variable elements ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # List of names of *all* components, supercomps
  pars$supercomps <- c(pars$grps, pars$subs, 'VFA', pars$comps)
  # And supercomps plus gases, at least to check stoich matrix and trim missing products
  pars$supercomps_gases <- c(pars$supercomps, pars$gases)

  # Sort out stoichiometry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Fermentation stoichiometry
  # Two possibilities: 1) NULL -> assume 1.0, 2) given
  # If missing, assume only VFA is produced
  if (is.null(pars$fstoich)) {
    # Fill in missing stochiometry
    pars$fstoich <- matrix(
      rep(1, length(pars$subs)),
      nrow = 1,
      dimnames = list(c('VFA'), c(pars$subs))
    )
  } 
  # Else given

  # Add substrates (consumption) to fstoich matrix as rows
  pars$fstoich <- comb_fstoich(pars$fstoich, pars$subs)

  # Add microbial yield to mstoich matrix as rows and subtract from product formation)
  pars$mstoich <- comb_mstoich(pars$mstoich, pars$yield)

  # Drop any gases not actually produced
  # comps treated differently because they could be conservative, or even accumulate if user wants
  pars$gases <- pars$gases[pars$gases %in% rownames(pars$mstoich)]

  # Check for reactants or products not included in simulation
  if (any(! rownames(pars$mstoich) %in% pars$supercomps_gases)) {
    missing <- rownames(pars$mstoich)[!rownames(pars$mstoich) %in% pars$supercomps_gases]
    stop('Reactant(s)/product(s) present in mstoich matrix but not a solute or gas: ', missing)
  }

  # Check for components or gases not included in mstoich matrix, add if missing
  if (any(!c(pars$comps, pars$gases) %in% rownames(pars$mstoich))) {
    warning('Expanding mstoich matrix with new rows to include missing comps or gases.')
    missing <- c(pars$comps, pars$gases)[!c(pars$comps, pars$gases) %in% rownames(pars$mstoich)]
    pars$mstoich <- rbind(pars$mstoich, matrix(0, nrow = length(missing), ncol = ncol(pars$mstoich), dimnames = list(missing, colnames(pars$mstoich))))
  }

  # Check for missing COD conversion factors
  if (any(!c(pars$comps, pars$gases) %in% names(pars$COD_conv))) {
    missing <- c(pars$comps, pars$gases)[!c(pars$comps, pars$gases) %in% names(pars$COD_conv)]
    stop('COD_conv missing some elements: ', paste(missing, collapse = ', '))
  }

  # Create biomass death stoichiometry matrix
  pars$dstoich <- make_dstoich(pars)

  # Extend ks matrix if not full (typically not)
  pars$ksmat <- fix_ksmat(pars$ksmat, pars$mstoich)

  # Solutes and particle groups in slurry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # VFA intermediate is always present; force name regardless of user input (expect no name)
  names(pars$VFA_fresh) <- names(pars$VFA_init) <- 'VFA'
  # Solutes (exclude gases from components) (and does not include particulate substrates)
  pars$sols <- unique(c(pars$comps[!names(pars$comps) %in% pars$gases], 'VFA'))
  # Grouping for easier rates() calculations
  pars$conc_fresh <- c(pars$xa_fresh, pars$sub_fresh, pars$VFA_fresh, pars$comp_fresh)
  pars$conc_init <- c(pars$xa_init, pars$sub_init, pars$VFA_init, pars$comp_init)

  # All comps (sols is part of comps) plus subs plus grps should have a COD conversion value for COD balance!
  # Units: g COD' per g whatever comp units are 
  # VFA is defined as 1
  pars$COD_conv['VFA'] <- 1
  # subs typically (always?) 1 but user could put in something else in fstoich matrix
  pars$COD_conv[pars$subs] <- pars$fstoich[1, ][pars$subs]
  # grps intended to be 1
  pars$COD_conv[pars$grps] <- 1
  # Values should be entered in chem_pars for all comps

  # Check for missing
  if (any(! pars$supercomps %in% names(pars$COD_conv))) {
    stop(paste0('Missing COD_conv for ', pars$supercomps[!pars$supercomps %in% names(pars$COD_conv)], '.'))
  }

  # Substrates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pars$n_subs <- length(pars$subs)
  pars$i_subs <- length(pars$grps) + 1:length(pars$subs)
  
  # Other constants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pars$g_NH4 <- 0.7
  pars$temp_standard <- 298
  pars$temp_zero <- 273
  pars$pH_floor <- 7

  # Conversions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Convert temperature constants to K if needed
  pars <- tempsC2K(pars, cutoff = 200)

  # Convert some supplied parameters
  # Maximum slurry mass in kg
  pars$max_slurry_mass <- pars$storage_depth * pars$area * pars$dens
  pars$resid_mass <- pars$resid_depth / pars$storage_depth * pars$max_slurry_mass

  # Respiration precomputation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pars$has_resp <- !is.null(pars$O2_flux) && pars$O2_flux > 0
  if (pars$has_resp && is.null(pars$rstoich)) {
    pars$rstoich <- c(VFA = -1)
  }

  # Inhibition precomputation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (!is.null(pars$ic0)) {
    # Slopes
    pars$ics <- 1 / (pars$ic100 - pars$ic0)

    # Figure out where inhibitors comes from
    rn <- rownames(pars$ic0)
    is_conc <- grepl('_conc$', rn)
    pars$inhib_y_rows <- which(is_conc)
    pars$inhib_y_vars <- gsub('_conc$', '', rn[is_conc])
    pars$inhib_p_rows <- which(!is_conc)
    pars$inhib_p_vars <- rn[!is_conc]
  }

  return(pars)
 
}

# Apply all keyword stuff
expand_pars <- function(
  pars, 
  elnms, 
  parnms
) {
  
  # Trim to only those elements in pars
  parnms <- parnms[parnms %in% names(pars)]

  for (i in parnms) {
    ppo <- pars[[i]]
    p_nms <- names(pars[[i]])

    # 'default' keyword fills in missing values
    if (any(p_nms == 'default')) {
      pars[[i]][elnms] <- pars[[i]]['default']
      if (any(p_nms != 'default')) {
        pars[[i]][p_nms[p_nms != 'default']] <- ppo[p_nms[p_nms != 'default']]
      }
    }

   # Fix order, drop default element if present, drop unused names
    pars[[i]] <- pars[[i]][elnms]
    # Check for missing values
    if (any(is.na(pars[[i]]))) {
      warning('Size-variable parameter problem: Missing element(s) in ', i, '.')
    }
  }

  return(pars)
  
}
  

combine_var_pars <- function(var_pars) {
  
  # If there is more than one data frame in var_pars, combine into one
  if (length(var_pars) > 1) {
    for (i in 2:length(var_pars)) {
      dv <- var_pars[[i]]
      nv <- names(var_pars)[[i]]
      ll <- list()
      # Collapse any multiple columns into one column
      for (j in 1:nrow(dv)) {
        # It is challenging to get list into each element using indexing
        x <- as.numeric(dv[j, -1, drop = FALSE])
        names(x) <- names(dv)[-1]
        ll[[j]] <- x
      }
      # Put combined columns into a new data frame with time
      dd <- dv[, 1, drop = FALSE]
      dd[[nv]] <- ll
      # And merge (any NAs will be fixed in clean_series())
      var_pars$var <- merge(var_pars$var, dd, by = 'time', all = TRUE)
    }
    # Take only the var element (but note that var_pars is still a list--see [] not [[]]) because all other info from other elements are now in it
    # And var_pars must remain a list (see single brackets) to avoid duplicate par elements from combining into par
    var_pars <- var_pars['var']
  }

  # Check for time and at least one other column name
  if (ncol(var_pars$var) == 1 || names(var_pars$var)[1] != 'time') {
    stop('Optional var_pars element var must be a data frame with first column time and at least one additional column.')
  }

  return(var_pars)

}



fix_add_pars <- function(pars, add_pars) {

  # Note: pe.pars = add_pars that use par.element approach, these are converted to normal (simple) add_par elements here
  # Note: sa.pars = normal (simple) add_pars that do not need to be converted
  # Note: Use of [] vs. [[]] affect how code works--needs to work for both lists and vector pars
  # Note: par.element approach is only designed to work for vector elements
  if (!is.null(add_pars) && length(add_pars) > 0 && any(ii <- grepl(pars$par_key, names(add_pars)))) {
    pe.pars <- add_pars[ii]
    sa.pars <- add_pars[!ii]
    apnames <- names(pe.pars)
    pe.pars[!grepl('[A-Za-z]', pe.pars)] <- as.numeric(pe.pars[!grepl('[A-Za-z]', pe.pars)])
    split.pars <- strsplit(apnames, pars$par_key)
    pnames <- sapply(split.pars, '[[', 1)
    enames <- sapply(split.pars, '[[', 2)
    names(pe.pars) <- enames
    pe.pars <- split(pe.pars, pnames)
    add_pars <- c(sa.pars, pe.pars)
  }

  # If any additional parameters were added (or modified) using add_pars, update them in pars list here
  # But grp_pars and sub_pars work differently than the others because of the default = keyword
  # Needs to work in a case where default is all but e.g., m1 is given in add_pars
  grp_par_nms <- c('yield', 'xa_fresh', 'xa_init', 'd_max', 'qhat_opt', 'T_opt', 'T_min', 'T_max')
  grp_par_nms <- grp_par_nms[grp_par_nms %in% names(pars)]
  sub_par_nms <- c('T_opt_hyd', 'T_min_hyd', 'T_max_hyd', 'hydrol_opt', 'sub_fresh', 'sub_init')
  sub_par_nms <- sub_par_nms[sub_par_nms %in% names(pars)]
  if (!is.null(add_pars) && length(add_pars) > 0) {
    if (any(bad.names <- !names(add_pars) %in% names(pars))) {
      stop ('Some `add_pars` names not recognized as valid parameters: ', names(add_pars)[bad.names]) 
    }
    # Add in pars (or replace existing elements unless it is time series data added)
    for (i in names(add_pars)) {
      if (inherits(add_pars[[i]], 'matrix')) {
	# For stoichiometry matrices
        pars[[i]] <- add_pars[[i]]
      } else if (!is.data.frame(add_pars[[i]]) && length(pars[[i]]) > 1) {
        pars[[i]][names(add_pars[[i]])] <- unlist(add_pars[[i]])
      } else {
        def <- pars[[i]]['all']
        pars[[i]] <- add_pars[[i]]
        if (i %in% c(grp_par_nms, sub_par_nms)) {
          pars[[i]]['default'] <- def
        }
      }
    }
  }
  
  # Unlike others, grps and subs in add_pars *will* override default vector (i.e., can be used to remove groups)
  # This must be clear in the help file!
  if ('grps' %in% names(add_pars)) {
    pars$grps <- add_pars$grps
  }
  if ('subs' %in% names(add_pars)) {
    pars$subs <- add_pars$subs
  }

  return(pars)

}

# Converts temperature parameters from C to K if values look like K (> cutoff, e.g., 270)
tempsC2K <- function(pars, cutoff) {

  for (i in grep('T_', names(pars))) {
    pars[[i]][pars[[i]] < cutoff] <- pars[[i]][pars[[i]] < cutoff] + 273.15
  }

  return(pars)
}


# Calculate all temperature-dependent parameters
calc_temp_pars <- function(pars, y) {
  
  pars$temp_K <- pars$temp_C + 273.15
  
  # Hydrolysis rate (vectorized)
  pars$h_rate[pars$subs] <- q10(pars$temp_K, pars$h_rate_ref, pars$h_rate_q10, pars$h_rate_max_temp)

  # Microbial substrate utilization rate (vectorized calculation)
  pars$qhat[pars$grps] <- CTM(pars$temp_K, pars$T_opt, pars$T_min, pars$T_max, pars$qhat_opt)

  # Temperature effect on ks values
  pars$ksmat_t <- pars$ksmat * (0.8157 * exp(-0.063 * pars$temp_C))

  # Death rate
  if (pars$temp_C > 40) {
    pars$d_rate <- pars$d_max
  } else {
    pars$d_rate <- CTM(pars$temp_K, 313, 273, 325, pars$d_max)
  }

  return(pars)

}

# CTM function
CTM <- function(temp_K, t_opt, t_min, t_max, y_opt) {

  # Check lengths
  if (length(temp_K) != 1) {
    stop('Length of temp_K must be 1')
  }
  if (!all.equal(length(t_opt), length(t_min), length(t_max), length(y_opt))) {
    stop('Length of t_opt, t_min, t_max, and y_opt must be identical')
  }

  # When temp_K outside bounds y = 0
  y <- numeric(length(t_opt))
  in_range <- temp_K[1] > t_min & temp_K[1] < t_max
  
  if (any(in_range)) {

    # Select which formula to use
    use_first <- (t_opt - t_min) < (t_max - t_min) / 2
    
    # First formula
    y1 <- y_opt * ((temp_K[1] - t_min) * (temp_K[1] - t_max)^2) /
      ((t_opt - t_max) * ((t_opt - t_max) * (temp_K[1] - t_opt) - (t_opt - t_min) * (t_opt + t_max - 2*temp_K[1])))
    
    # Second formula
    y2 <- y_opt * ((temp_K[1] - t_max) * (temp_K[1] - t_min)^2) /
      ((t_opt - t_min) * ((t_opt - t_min) * (temp_K[1] - t_opt) - (t_opt - t_max) * (t_opt + t_min - 2*temp_K[1])))
    
    # Combine using ifelse and enforce in-range
    y[in_range] <- pmax(0, ifelse(use_first, y1, y2)[in_range])
  }
  
  return(y)
}

# Arrhenius function
arrhenius <- function(temp_K, A, E, arr_max_temp_K = 313, R = 8.314){
  
  if (temp_K >= arr_max_temp_K) temp_K <- arr_max_temp_K

  y <- A * exp(-E/(R * temp_K))
  
  return(y)
}

# Q10 function for hydrolysis rate, possibly others
q10 <- function(temp_K, yref, qpar, temp_ref = 273.15 + 20, max_temp = 313) {

  temp_K[temp_K >= max_temp] <- max_temp

  y <- yref * qpar^((temp_K - temp_ref) / 10)

  return(y)

}

# Inhibition
update_inhib <- function(pars, y) {

  if (is.null(pars$ic0) || is.null(pars$ic100)) {
    return(pars$qhat)
  }

  # Build inhibitor value vector using precomputed indices
  x <- numeric(nrow(pars$ic0))
  x[pars$inhib_y_rows] <- y[pars$inhib_y_vars] / y['slurry_mass']
  x[pars$inhib_p_rows] <- unlist(pars[pars$inhib_p_vars])

  # Inhibition matrix (x recycled row-wise across columns of ic0 so this works without turning it into a matrix)
  im <- pars$ics * (x - pars$ic0)
  im[im < 0] <- 0
  im[im > 1] <- 1
  im <- 1 - im

  # Apply to qhat
  qhat <- pars$qhat
  # Use product of columns in im matrix 
  qhat[colnames(im)] <- qhat[colnames(im)] * apply(im, 2, prod)

  # And return
  return(qhat)
 
}
