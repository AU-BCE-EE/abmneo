# Stoichometry functions

mol_mass <- function(form, elements = NULL) {

  ## Check argument
  #checkArgClassValue(form, 'character')

  # Loop through all elements in form
  mmass <- NULL
  for(f in form) {
    # If and only if first letter of form is lowercase, entire string is capitalized
    if(grepl('^[a-z]', f)) f <- toupper(f) 

    # Get coefficients of formula
    fc <- read_formula(f)

    if (!is.null(elements)) {
      fc <- fc[intersect(names(fc), elements)]
    }

    # Check for unidentified element
    if(any(!names(fc) %in% names(atom.weights))) stop('One or more elements in \"form\" is not in the database. You can add it to the \"atom.weights\" vector if you want to modify the function code. Otherwise send a request to sasha.hafner@bce.au.dk.')

    # Calculate molar mass, using names of fc for indexing
    mmass <- c(mmass, sum(atom.weights[names(fc)]*fc))
  }

  return(mmass)
}

# Returns COD per mol substrate
calc_COD <- function(form) {

  # If and only if first letter of form is lowercase, entire string is capitalized
  if(grepl('^[a-z]', form)) form <- toupper(form)
  # Read formula (function not vectorized)
  fc <- read_formula(form, elements = c('C', 'H', 'O', 'N'))
  # Calculate COD based on Rittmann and McCarty
  COD <- as.vector((2*fc['C'] + 0.5*fc['H'] - 1.5*fc['N'] - fc['O']) * mol_mass('O'))

  return(COD)
}


# Get mass conversion factor to go from moles of component to mass COD, N, C, S, or total, in that order
get_mass_conv <- function(form) {

  # Remove p and m (+/-)
  form <- gsub('p$|m$', '', form)
  
  cod <- calc_COD(form)
  fn <- read_formula(form)
  
  if (cod > 0) {
    cf <- cod
  } else if ('N' %in% names(fn)) {
    cf <- mol_mass(form, elements = 'N')
  } else if ('C' %in% names(fn)) {
    cf <- mol_mass(form, elements = 'C')
  } else if ('S' %in% names(fn)) {
    cf <- mol_mass(form, elements = 'S')
  } else {
    cf <- mol_mass(form)
  }

  return(cf)

}

# Check COD balance
check_COD <- function(dat, 
                     grps,
                     subs,
                     COD_conv,
                     fstoich,
                     rtol = 0.001
                    ) {

  first <- unlist(dat[1, ])
  last <- unlist(dat[nrow(dat), ])

  CODin <- sum(last[['COD_load']], first[grps], first[subs] * fstoich['CH3COOH', subs], first[['CH3COOH']])
  CODeff <- sum(last[paste0(subs, '_eff')] * fstoich['CH3COOH', subs]) + last[['CH3COOH_eff']]
  CODemis <- last[['CH4']] * COD_conv[['CH4']]
  CODrem <- sum(last[grps], last[subs] * fstoich['CH3COOH', subs], last[['CH3COOH']])
  bal <- CODin - CODeff - CODemis - CODrem
  rbal <- bal / CODin

  if (abs(rbal) > rtol) {
    warning('COD balance is off by ', signif(100 * rbal, 2), '%')
    return(invisible(rbal))
  } 

  return(invisible(rbal))

}

comb_fstoich <- function(fstoich, subs) {

  # Short name
  fsc <- fstoich

  # Substrate matrix
  sm <- diag(-1, ncol(fsc), ncol(fsc))
  rownames(sm) <- colnames(fsc)
  combstoich <- rbind(fsc, sm)

  return(combstoich)

}

comb_mstoich <- function(mstoich, yield) {

  # Subtract microbial yields from product formation
  # Note that this is not simply subtraction because some products may not be in COD units!
  msc <- mstoich
  msc[msc > 0] <- (msc * matrix(1 - yield, nrow = nrow(msc), ncol = length(yield), byrow = TRUE))[msc > 0]

  # Yield matrix
  ym <- diag(yield, ncol(msc), ncol(msc))
  rownames(ym) <- names(yield)
  combstoich <- rbind(mstoich, ym)

  return(combstoich)

}

fix_ksmat <- function(ksmat, mstoich) {

  ks <- 0 * mstoich + 1e10
  ks[rownames(ksmat), colnames(ksmat)] <- ksmat

  return(ks)

}
