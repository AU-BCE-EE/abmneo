
To work on:
* does var_pars need first element named var? How about the multi-element pars?
* Don't need this in pars.R do we?

  # Multiple microbial groups ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # NTS: need to sort out how this works with above mess for add_pars with grps
  if (is.null(pars$meths)) {
    pars$meths <- pars$grps[! pars$grps %in% pars$srs]
  } else {
    pars$grps <- unique(c(pars$grps, pars$meths, pars$srs))
  }


