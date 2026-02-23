# rates() and related functions

rates <- function(t, y, parms) {
    
  # Short name for parms to make indexing in code below simpler 
  p <- parms

  # Calculation speciation, used for emission and inhibition
  p <- calc_spec(p, y)

  # Determine inhibition reductions
  p <- calc_inhib(p, y)

  # Initialize vectors with derivative components, all with same order of y elements
  inflow <- meth <- hyferm <- 0 * y

  # Inflow from slurry addition
  inflow[names(p$conc_fresh)] <- p$conc_fresh * p$slurry_prod_rate
  inflow[c('slurry_mass', 'slurry_load')] <- 1 * p$slurry_prod_rate
  inflow['COD_load'] <- sum(inflow[names(p$conc_fresh)])

  # VFA consumption rates, CH4 production, and growth
  # First utilization rate
  browser()
  p$mstoich
  i = 1
  j = 1
  rut <- 0 * p$qhat

  # Option 1
  rut <- p$qhat * y['CH3COOH'] / (y['CH3COOH'] + y['slurry_mass'] * p$kss)

  # Option 2
  for (i in seq_len(ncol(p$mstoich))) {
    for (j in seq_len(nrow(p$mstoich[, i, drop = FALSE]))) {
      ss <- p$mstoich[j, i, drop = FALSE]
      if (ss < 0) {
        grp <- colnames(ss)
        rct <- rownames(ss)
        rut[grp] <- (y[rct] / (y[rct] + y['slurry_mass'] * p$ksmat[j, i]) * rut[grp]
      }
    }
    rut[grp] <- qhat[grp] * rut[grp]
  }

  # Option 1
  rut <- p$qhat * y['CH3COOH'] / (y['CH3COOH'] + y['slurry_mass'] * p$kss)
  # And consumption, growth, production all in one matrix operation
  meth[rownames(p$mstoich)] <- p$mstoich %*% rut * y['slurry_mass']
  # Convert CH4 from g COD / d to g C / d
  meth['CH4'] <- meth['CH4'] / p$COD_conv['CH4']

  # Hydrolysis of particulate substrates and fermentation
  hyferm[p$subs] <- - p$alpha[p$subs] * y[p$subs]
  # Production of arbitrary products based on specified fermentation stoichiometry (can omit components)
  hyferm[rownames(p$fstoich)] <- - p$fstoich %*% hyferm[colnames(p$fstoich)]
  
  # Add vectors to get derivatives
  # All elements in g/d as COD except 
  #   * slurry_mass (kg/d as fresh slurry mass)
  #   * CH4 (g/d as CH4)
  #   * user-defined solutes other than VFA (as C, N, or S as described in documentation)
  ders <- inflow + meth + hyferm

  return(list(ders, c(CH4_emis_rate = meth[['CH4']], temp_C = p$temp_C, pH = p$pH)))

}

# Copy time-variable parameters into their normal par position for an interval
update_var_pars <- function(
  series,
  pars, 
  y, 
  i
) {

  vdat <- series[, names(series) != 'time', drop = FALSE]

  for (j in 1:ncol(vdat)) {
    pn <- names(vdat)[j]
    newval <- vdat[i, j] 
    if (is.list(newval)) {
      newval <- newval[[1]]
    }
    pars[[pn]] <- newval
  }

  # Temperature and temperature-dependent pars could change
  # Convert temperature to K in case any C values were given in var
  pars <- tempsC2K(pars, cutoff = 200)

  # Calculate temperature-dependent par values
  pars <- calc_temp_pars(pars, y)

  return(pars)

}

# Calculate inhibition from pars values
calc_inhib <- function(pars, y) {
  
  # Skip all this if there are no inhibition parameters
  if (!is.null(pars$ilwr) & !is.null(pars$iupr)) {
    # Get concentrations of inhibitor species
    concs <- pars$conc_sp[colnames(pars$ilwr)]
    pars$ired <- inhib(pars$ilwr, pars$iupr, concs)
  } else {
    pars$ired <- rep(1, length(pars$grps))
    names(pars$ired) <- pars$grps
  }

  return(pars)

}

# Calculate chemical speciation
# Always includes VFA
# May include optional parameter-defined species 
calc_spec <- function(pars, y) {

  # Get concentrations of all solutes
  concs <- y[pars$sols] / y['slurry_mass']

  # Save totals 
  pars$conc_tot <- concs
  
  # Calculate and save species concentrations
  pars$conc_sp <- as.numeric(concs[pars$mspec[pars$specs]]) * 1 / (1 + 10^(-pars$lkan[pars$specs] - pars$pH))
  pars$conc_sp[pars$mspec[pars$spec]] <- pars$conc_tot[pars$mspec[pars$spec]] - pars$conc_sp[pars$spec]
  
  return(pars)
}

# Calculate log Ka for speciation
calc_ka <- function(pars) {

  # Caclulate numeric pKa from expressions and temp_K
  temp_K <- pars$temp_K
  pars$lkan <- as.numeric(lapply(pars$lka, function(x) eval(parse(text = x))))
  names(pars$lkan) <- names(pars$lka)

  return(pars)
}


calc_volat <- function(p, vi) {

  volat <- 0 * vi

  # Emission of other species                                           
  if (!is.null(p$kl)) {
    volat[paste0(names(p$kl), '_emis_cum')] <- p$kl * p$conc_sp[names(p$kl)] * p$area
    # Remove emitted amount from component pool
    volat[p$mspec[names(p$kl)]] <- - volat[paste0(names(p$kl), '_emis_cum')]
  } 

  return(volat)
  
}


# Calculate inhibition
# Matrix approach
inhib <- function(im0, 
                  im1, 
                  iconc, 
                  combine = prod) {

  # Extract groups and inhibitors for dimensions but mainly for debugging
  inhibs <- colnames(im0)
  grps <- rownames(im0)

  # Get just those species that have inhibition parameters
  iconc <- iconc[inhibs]

  # Convert concentrations to matrix
  cmat <- matrix(rep(iconc, length(grps)), 
                 nrow = length(grps), 
                 byrow = TRUE,
                 dimnames = list(grps, 
                                 inhibs))
  # Slope matrix
  smat <- (1 - 0) / (im1 - im0)

  # Reduction matrix
  rmat <- 1 - (cmat - im0) * smat
  rmat[rmat > 1] <- 1
  rmat[rmat < 0] <- 0

  # Get product by microbial group
  red <- apply(rmat, 1, combine)

}
