
add_pars <- function(pars, add_pars) {

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
