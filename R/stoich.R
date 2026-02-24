# Stoichometry functions

# Figure out stoichiometry matrix of fermentation substrates from elemental formula
get_stoich <- function(pars) {

  # Get molar stoichiometric coefficients
  # Vectorize, return matrix without substrate (for 1 mole substrate)
  res <- lapply(as.list(pars$forms), predFerm)
  # Align names and sort before combining in matrix
  nn <- unique(unlist(lapply(res, names)))
  for (i in 1:length(res)) {
    res[[i]][nn[!nn %in% names(res[[i]])]] <- 0
    res[[i]] <- res[[i]][nn]
  }
  st <- matrix(unlist(res), ncol = length(pars$forms), byrow = FALSE)
  rownames(st) <- names(res[[1]])
  colnames(st) <- pars$forms

  # Drop 0
  st <- st[rowSums(st) != 0, , drop = FALSE]
  
  # Drop water (ignored, treated as conservative in system)
  st <- st[rownames(st) != 'H2O', , drop = FALSE] 

  # Switch from chemical formulas of columns to substrate names
  colnames(st) <- names(pars$forms)

  # Switch to master species names for products (some match)
  rownames(st) <- pars$mspec[rownames(st)]

  # Adjust coefficients to COD mass, N mass, C mass, S mass, or total mass
  for (i in 1:nrow(st)) {
    ff <- rownames(st)[i]
    st[i, ] <- st[i, ] * pars$mcf[ff]
  }
  
  for (i in 1:ncol(st)) {
    ff <- colnames(st)[i]
    st[, i] <- st[, i] * 1 / pars$mcf[ff]
  }

  return(st)

}

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

# Fermentation stoichiometry
# Example calls:
# source('read_formula.R')
# predFerm('C6H10O5', acefrac = 0, fs = 0.1)
# predFerm('C6H10O5', acefrac = 0.5, fs = 0.1)
# predFerm('C6H10O5', acefrac = 1)
# predFerm('C6H10O5', acefrac = 1)

# Function to get stoichiometry for custom organic reaction (O-19)
customOrgStoich <- function(
  form, 
  elements =  c('C', 'H', 'O', 'N'),
  dover = FALSE
  ) {
  
  fc <- read_formula(form, elements)

  # Use symbols from O-19 in R&M
  n <- as.numeric(fc['C'])
  a <- as.numeric(fc['H'])
  b <- as.numeric(fc['O'])
  cc <- as.numeric(fc['N'])
  
  # dover for d override deals with case of no available electrons
  if (dover) {
    d <- 1
    H.c <- 0
  } else {
    d <- 4 * n + a - 2 * b - 3 * cc
    H.c <- -1
  }

  if (d == 0) {
    return(form)
  }
  
  # Put together
  #rr <- c(CO2 = - (n - cc) / d, NH4. = - cc / d, HCO3. = - cc / d, H. = -1, H2O = (2*n - b + cc) /d)
  rr <- c(CO2 = - n / d, NH3 = - cc / d, H. = H.c, H2O = (2*n - b) / d)
  rr[form] <-  1/d

  return(rr)
  
}

# General Rittmann and McCarty stype stoichiometry calculations
RMStoich <- function(subform, rd, ra, rc, fs, dropzero, dropsub, order, tol) {

  ii <- unique(names(c(rd, rc, ra)))

  # Blanks
  rd[ii[!ii %in% names(rd)]] <- 0
  rc[ii[!ii %in% names(rc)]] <- 0
  ra[ii[!ii %in% names(ra)]] <- 0

  # Order
  rd <- rd[ii]
  rc <- rc[ii]
  ra <- ra[ii]

  fe <- 1 - fs
  
  # Combine
  rtot <- fe * ra + fs * rc  - rd

  # Drop substrate if requested
  if (isTRUE(dropsub)) {
    # Adjust coefficients to 1 mol substrate
    rtot <- - rtot / rtot[subform]
    rtot <- rtot[names(rtot) != subform] 
  }

  rtot[abs(rtot) < tol] <- 0 
  
  # Drop empty elements
  if (dropzero) {
    rtot <- rtot[rtot != 0]
  }

  if (!is.na(order[1]) && tolower(order[1]) == 'sort') {
    rtot <- rtot[order(rtot < 0, abs(rtot), decreasing = TRUE)]
  } else if (!is.na(order[1]) && all(sort(order) == sort(names(rtot)))) {
    rtot <- rtot[order]
  } else if (inherits(order, 'logical') && isFALSE(order)) {
    # Skips order
  } else if (!is.na(order[1])) {
    warning('order argument ignored')
  }

  return(rtot)

}

predFerm <- function(
  subform = NULL,           # Character chemical formula of substrate
  biomassform = 'C5H7O2N',  # Biomass empirical formula
  acefrac = 1,              # Acetate (vs. H2) fraction
  fs = 0,                    # Fraction substrate going to cell synthesis, fs in Rittmann and McCarty
  elements = c('C', 'H', 'O', 'N'),
  order = 'sort',
  dropzero = TRUE,
  dropsub = TRUE,
  tol = 1E-10
  ) {

  if (biomassform %in% subform) {
    stop('The biomassform formula is the/a substrate! Change one (change in just order is OK).')
  }

  # Donor half reaction
  rd <- customOrgStoich(subform, elements = elements)

  # If donor has no available electrons
  if (length(rd) == 1 && rd == subform) {
    rtot <- - customOrgStoich(subform, elements = elements, dover = TRUE)
    rtot <- c(rtot, CH3COOH = 0, H2 = 0)
    # Drop substrate
    rtot <- rtot[names(rtot) != subform]
    return(rtot)
  }
  
  ## Donor needs H2 and CH3COOH
  #for (sp in c('H2', 'CH3COOH')) {
  #  if (! sp %in% names(rd)) {
  #    rd[sp] <- 0
  #  }
  #}
  #rd <- rd[sort(names(rd))]

  # Synthesis half reaction
  rc <- customOrgStoich(biomassform, elements = elements)

  # Acceptor reactions
  # Acetate production
  #raa <- c(CO2 = - 1/8, HCO3. = - 1/8, H. = -1, CH3COO. = 1/8, H2O = 3/8)
  raa <- c(CO2 = - 1/4, H. = -1,  CH3COOH = 1/8, H2O = 1/4, H2 = 0)
  # Hydrogen|production |            |            |          |
  rah <- c(CO2 = 0,     H. = - 1, CH3COOH = 0,   H2O = 0,   H2 = 1/2)
  
  # Acceptor reaction
  ra <- acefrac * raa + (1 - acefrac) * rah 

  rtot <- RMStoich(subform = subform, rd = rd, ra = ra, rc = rc, fs = fs, 
                   dropzero = dropzero, dropsub = dropsub, order = order, 
                   tol = tol)

  return(rtot)

}

# For methanogenesis, only difference from ferm is acceptor reaction
# These should be combined--too much code copied!
predMethan <- function(
  subform = NULL,           # Character chemical formula of substrate
  biomassform = 'C5H7O2N',  # Biomass empirical formula
  fs = 0,                   # Fraction substrate going to cell synthesis, fs in Rittmann and McCarty
  elements = c('C', 'H', 'O', 'N'),
  order = 'sort',
  dropzero = TRUE,
  dropsub = FALSE,
  tol = 1E-10
  ) {

  if (biomassform %in% subform) {
    stop('The biomassform formula is the/a substrate! Change one (change in just order is OK).')
  }

  # Vectorize, return matrix without substrate (for 1 mole substrate)
  if (length(subform) > 1) {
    res <- lapply(as.list(subform), predMethan, biomassform = biomassform, 
                  fs = fs, elements = elements, order = FALSE, dropzero = FALSE, dropsub = TRUE, 
                  tol = tol)
    resmat <- matrix(unlist(res), nrow = length(subform), byrow = TRUE)
    colnames(resmat) <- names(res[[1]])
    rownames(resmat) <- subform

    if (dropzero) {
      resmat <- resmat[, colSums(resmat) != 0]
    }
    return(resmat)
  }

  # Donor half reaction
  rd <- customOrgStoich(subform, elements = elements)

  ## Donor needs H2 and CH3COOH
  #for (sp in c('H2', 'CH3COOH')) {
  #  if (! sp %in% names(rd)) {
  #    rd[sp] <- 0
  #  }
  #}
  #rd <- rd[sort(names(rd))]

  # Synthesis half reaction
  rc <- customOrgStoich(biomassform, elements = elements)

  # Acceptor reaction
  ra <- c(CO2 = - 1/8, H. = -1, CH4 = 1/8, H2O = 1/4)

  rtot <- RMStoich(subform = subform, rd = rd, ra = ra, rc = rc, fs = fs, 
                   dropzero = dropzero, dropsub = dropsub, order = order, 
                   tol = tol)

  return(rtot)

}


# Modified: 4 April 2016 SDH
# NTS: apparently *not* vectorized! Revisit. Had to modify calc_COD 10 Mar 2017 to fix it.

read_formula <- function(
  form,
  elements = NULL,        # Set of elements returned, all others ignored, e.g., c('C', 'H', 'N', 'O')
  min.elements = NULL,    # Minimum set of elements, will return error if these at least are not included
  cdigits = 6,
  value = 'numeric'       # Type of output, 'numeric' for named vector, 'shortform' for shortened formula
) {

  form.orig <- form

  # Remove spaces
  form <- gsub(' ', '', form)

  # Add implied coefficients of 1 (also after ")")
  form <- gsub('([a-zA-Z\\)])([A-Z\\)\\(])', '\\11\\2', form)
  form <- gsub('([a-zA-Z\\)])([A-Z\\)\\(])', '\\11\\2', form) # Repeated for e.g., COOH
  form <- gsub('([a-zA-Z\\)])$', '\\11', form)

  # Find parentheses and remove them, multipying coefficients inside by coefficient at end
  # So (CH2)2 ---> C2H4
  # First add ( after N), e.g., (CH2)2CH3 ---> (CH2)2(CH3 for separation below
  form <- gsub('(\\)[0-9\\.]+)', '\\1(', form)
  # Drop extra (
  form <- gsub('^\\(', '', form)
  form <- gsub('\\($', '', form)
  form <- gsub('\\(\\(', '(', form)
  s1  <- strsplit(form, '\\(')[[1]]

  # Build up elementwise formula piecewise
  formpw <- NULL
  for(i in 1:length(s1)) {
    xx <- s1[i]
    if(grepl('\\)', xx)) {
      nn <- as.numeric(gsub('.+\\)','',xx))
      ff <- gsub('\\).+', '', xx)
      cc <- nn*as.numeric(strsplit(ff, '[A-Za-z]+')[[1]][-1])
      ee <- strsplit(ff, '[0-9.]+')[[1]]
      formpw <- paste0(formpw, paste0(ee, cc, collapse = ''))
    } else {
      formpw <- paste0(formpw, xx)
    }
  }

  form <- formpw

  # Extract integer coefficients 
  cc <- as.numeric(strsplit(form, '[A-Za-z]+')[[1]][-1])
  names(cc) <- strsplit(form, '[0-9.]+')[[1]]

  # Sort out elements to return
  if(is.null(elements)) elements <- unique(names(cc))
  fc <- numeric(length(elements))
  names(fc) <- elements

  # Fill in fc, summing elements of cc if required (if elements are repeated)
  for(i in elements) {
    for(j in 1:length(cc)) {
      if(names(cc)[j]==i) fc[i] <- fc[i] + cc[j]
    }
  }

  # Simplify form based on fc (for output only)
  # format() is for the rare case with something like C0.00001, to avoid C1e-5 which will result in an error
  form <- paste0(names(fc), format(signif(fc/min(fc), cdigits), scientific = FALSE), collapse = '')
  # Drop spaces that come in with format(...scientific = FALSE) (scipen fix)
  form <- gsub(' ', '', form)
  # And drop coefficients of 1
  form <- gsub('([a-zA-Z])1([a-zA-Z])', '\\1\\2', form)
  form <- gsub('([a-zA-Z])1$', '\\1\\2', form)
  
  # Check for minimum set of elements
  if(!is.null(min.elements)) if(any(!min.elements %in% names(fc)) | any(fc[min.elements] == 0)) stop('Minimum elements required are ', min.elements, ' (from min.elements argument), but form is ', form.orig, ', interpreted as ', form)

  if(value == 'numeric') return(fc)
  if(value == 'shortform') as.vector(form)

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

comb_stoich <- function(mstoich, yield) {

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
