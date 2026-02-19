# Sorts out parameters and packages them all together in the output list
# This is a central function that does a lot and is (unfortunately) complicated

packPars <- function(
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
    var_pars <- combineVarPars(var_pars)
  }

  # Combine pars to make extraction and pass to rates() easier ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pars <- c(mng_pars, man_pars, init_pars, grp_pars, sub_pars, chem_pars, inhib_pars, ctrl_pars, var_pars)

  # Add indicator for already variable inputs
  if(!is.null(var_pars) && !is.null(var_pars$var)) {
    pars$regular <- FALSE  
  } else {
    pars$regular <- TRUE
  }

  # Sort out add_pars and similar parameter inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pars <- add_pars(pars, add_pars)

  # Finish working with var_pars ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This must come after add_par block because approx_method (and other relevant pars?) could be set with add_pars
  pars <- fixVarDat(pars, days)
  pars <- calcProdRem(pars)

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
  grp_par_nms <- c("yield", "xa_fresh", "xa_init", "dd_rate", "ksv", "kss", "qhat_opt", "T_opt", "T_min", "T_max")
  grp_par_nms <- grp_par_nms[grp_par_nms %in% names(pars)]
  sub_par_nms <- c("T_opt_hyd", "T_min_hyd", "T_max_hyd", "hydrol_opt", "sub_fresh", "sub_init")
  sub_par_nms <- sub_par_nms[sub_par_nms %in% names(pars)]
  pars <- expandPars(pars = pars, elnms = pars$grps, parnms = grp_par_nms)
  pars <- expandPars(pars = pars, elnms = pars$subs, parnms = sub_par_nms)

  # Check grp arguments, including order of element names in some pars
  # After above block, this should be redundant
  checkGrpNames(pars)

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
  pars$mcf <- unlist(lapply(c(pars$form, pars$mspec), getMassConv))

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
    pars$stoich <- getStoich(pars)
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
