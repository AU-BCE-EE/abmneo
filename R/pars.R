
# Sorts out parameters and packages them all together in the output list
# This is a central function that does a lot and is (unfortunately) complicated
pack_pars <- function(
  structure,
  mng_pars,
  man_pars,
  init_pars,
  grp_pars,
  sub_pars,
  chem_pars,
  inhib_pars,
  ctrl_pars,
  var_pars,
  add_pars,
  days
) {

  # Move extra var_pars into first (var) element, possibly as lists within each data frame element~
  if(!is.null(var_pars) && !is.null(var_pars$var)) {
    var_pars <- combine_var_pars(var_pars)
  }

  # Combine pars to make extraction and pass to rates() easier ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pars <- c(mng_pars, man_pars, init_pars, grp_pars, sub_pars, chem_pars, inhib_pars, ctrl_pars, var_pars)

  # Sort out add_pars and similar parameter inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pars <- fix_add_pars(pars, add_pars)

  # Multiple microbial groups ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # NTS: need to sort out how this works with above mess for add_pars with grps
  if (is.null(pars$meths)) {
    pars$meths <- pars$grps[! pars$grps %in% pars$srs]
  } else {
    pars$grps <- unique(c(pars$grps, pars$meths, pars$srs))
  }

  # Fill in default values for grp_pars if keyword name `default` or `all` is used
  # Note: Microbial groups are defined by grps element
  # Note: `default` does *not* work with add_pars argument because grps are already defined in defaults
  # Note: But `all` *does* work
  # expandPars() will also sort out element order and drop excluded elements
  # NTS: Could these vectors of names be set in some kind of defaults?
  grp_par_nms <- c("yield", "xa_fresh", "xa_init", "dd_rate", "ksv", "kss", "qhat_opt", "T_opt", "T_min", "T_max")
  sub_par_nms <- c("T_opt_hyd", "T_min_hyd", "T_max_hyd", "hydrol_opt", "sub_fresh", "sub_init")
  pars <- expand_pars(pars = pars, elnms = pars$grps, parnms = grp_par_nms)
  pars <- expand_pars(pars = pars, elnms = pars$subs, parnms = sub_par_nms)

  # Check grp arguments, including order of element names in some pars
  # After above block, this should be redundant
  check_grp_names(pars)

  # O2 kl ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # O2 kl is handled differently from other kl values, so here it becomes a separate pars element
  if (!is.null(pars$kl) && 'O2' %in% names(pars$kl)) {
    pars$O2kl <- pars$kl['O2']
    pars$kl <- pars$kl[names(pars$kl) != 'O2']
    if (length(pars$kl) == 0) {
      pars$kl <- NULL
    }
  }

  # Size-variable elements ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # For size-variable parameters, get number of elements and indices 
  # NTS: I expect to change to a different approach, only using block below this one with names
  pars$n_mic <- length(pars$grps)
  pars$i_mic <- grep('^sr|^p|^m', pars$grps)
  pars$i_meth <- grep('^[mp]', pars$grps)
  pars$i_sr <- grep('^sr', pars$grps)
  pars$i_aer <- grep('^aer', pars$grps)
  pars$i_hyd <- grep('^hyd', pars$grps)

  # Get names of variable elements 
  # Remember pars$grps/pars$mics and pars$subs already exist (set in pars input)
  pars$meths <- pars$grps[pars$i_meth]
  pars$srs <- pars$grps[pars$i_sr]
  pars$aers <- pars$grps[pars$i_aer]
  pars$hyds <- pars$grps[pars$i_hyd]
  
  # Drop sulfate reducers if SO4-2 is not available ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (length(pars$srs) > 0 & 'SO4m2' %in% pars$comps) {
    pars$sromit <- FALSE 
  } else {
    pars$sromit <- TRUE
  }

  # Solutes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # CH3COOH is always present
  pars$sols <- unique(c(pars$comps, 'CH3COOH'))
  pars$conc_fresh <- c(pars$comp_fresh, pars$VFA_fresh)

  # Master species, fill in masters = masters
  mmspec <- pars$sols 
  names(mmspec) <- mmspec
  pars$mspec <- c(pars$mspec, mmspec)
  pars$mspec <- pars$mspec[!duplicated(names(pars$mspec))]

  # Mass conversion factors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Need to be above stoich matrix determination (below)
  pars$mcf <- unlist(lapply(c(pars$form, pars$mspec), get_mass_conv))

  # Sort out stoichiometry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Three possibilities: 1) NULL, 2) "calc", 3) given
  # If missing, assume only VFA is produced
  if (is.null(pars$stoich)) {
    # Fill in missing stochiometry
    pars$stoich <- matrix(
      rep(1, length(pars$subs)),
      nrow = 1,
      dimnames = list(c('CH3COOH'), c(pars$subs))
    )
  } else if (all(tolower(pars$stoich) == 'calc')) {
    # Or calculated from substrate chemical formulas
    pars$stoich <- get_stoich(pars)
  }
  # Else given

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

    # 'all' keyword is no longer encouraged, but functionality is still here
    if (any(p_nms == 'all')) {
      pars[[i]][elnms] <- pars[[i]]['all']
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
      if (nrow(dv) != nrow(var_pars$var)) {
        stop('var_pars has multiple elements with different sizes.')
      }
      if (!all(dv$time == var_pars$var$time)) {
        stop('Multiple elements in var_pars have different time values in first column.')
      }
      ll <- list()
      for (j in 1:nrow(var_pars$var)) {
        # It is challenging to get list into each element using indexing
        x <- as.numeric(dv[j, -1, drop = FALSE])
        names(x) <- names(dv)[-1]
        ll[[j]] <- x
      }
      var_pars$var[[nv]] <- ll
    }
    # Take only the var element (but note that var_pars is still a list--see [] not [[]])
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
  grp_par_nms <- c("yield", "xa_fresh", "xa_init", "dd_rate", "ksv", "kss", "qhat_opt", "T_opt", "T_min", "T_max")
  grp_par_nms <- grp_par_nms[grp_par_nms %in% names(pars)]
  sub_par_nms <- c("T_opt_hyd", "T_min_hyd", "T_max_hyd", "hydrol_opt", "sub_fresh", "sub_init")
  sub_par_nms <- sub_par_nms[sub_par_nms %in% names(pars)]
  if (!is.null(add_pars) && length(add_pars) > 0) {
    if (any(bad.names <- !names(add_pars) %in% names(pars))) {
      stop ('Some `add_pars` names not recognized as valid parameters: ', names(add_pars)[bad.names]) 
    }
    # Add in pars (or replace existing elements unless it is time series data added)
    for (i in names(add_pars)) {
      if (!is.data.frame(add_pars[[i]]) && length(pars[[i]]) > 1) {
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


# NTS: Outdated! Not all these elements are present.
check_grp_names <- function(pars) {

  if (
    !all.equal(
      names(pars$yield), 
      names(pars$xa_fresh), 
      names(pars$xa_init), 
      names(pars$decay_rate),
      names(pars$ks_coefficient), 
      names(pars$resid_enrich), 
      names(pars$qhat_opt), 
      names(pars$T_opt), 
      names(pars$T_min), 
      names(pars$T_max), 
      names(pars$ki_NH3_min), 
      names(pars$ki_NH3_max), 
      names(pars$ki_NH4_min), 
      names(pars$ki_NH4_max), 
      names(pars$pH_lwr), 
      names(pars$upr)
    )
  ) {
    stop('Microbial groups, i.e., element names in `grp_pars`, must match but do not.')
  }

}

# Converts temperature parameters from C to K if values look like K (> cutoff, e.g., 270)
tempsC2K <- function(pars, cutoff) {

  for (i in grep('T_', names(pars))) {
    pars[[i]][pars[[i]] < cutoff] <- pars[[i]][pars[[i]] < cutoff] + 273.15
  }

  return(pars)
}

