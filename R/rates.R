# rates() and related functions

rates <- function(t, y, parms) {
    
  # Short name for parms to make indexing in code below simpler 
  p <- parms

  # Calculation speciation, used for emission and inhibition
  p <- calc_spec(p, y)

  # Determine inhibition reductions
  p <- calc_inhib(p, y)

  # Get indices for microbial groups
  # Remember order in y and in all pars is identical (it is forced and checked)
  igrp <- which(p$grps %in% names(y))
  # And substrate
  isub <- which(p$subs %in% names(y))

  # Initialize vectors with derivative components, all with same order of y elements
  inflow <- rxn <- hydrol <- 0 * y

  # Other (temperature-dependent) derivative vectors, brought in with pars
  alpha <- p$alpha
  qhat <- p$qhat

  # Inflow from slurry addition
  # First only concentrations are set, and multiplied by inflow in last line
  inflow[igrp] <- p$xa_fresh
  inflow[isub] <- p$sub_fresh[isub]
  inflow[p$sols] <- p$conc_fresh[p$sols]
  inflow[c('slurry_mass', 'slurry_load')] <- 1
  inflow['COD_load'] <- sum(inflow[p$grps], inflow[isub] * p$stoich['CH3COOH', isub], inflow['CH3COOH'])
  inflow <- inflow * p$slurry_prod_rate

  # VFA consumption rates (g/d) and growth
  # Rate of substrate utilization
  # For speed in indexing, use integer indices not names
  for (i in igrp) {
    # Extract group's stoichiometry
    st <- p$mstoich[i, ]
    # Utilization rate
    ri <- which(names(st[st == -1]) == names(y))
    ut <- qhat[i] * y[i] * prod(y[ri]/y['slurry_mass']) / (p$ksv[i] + prod(y[ri]/y['slurry_mass'])) * y['slurry_mass']
    # Growth minus death
    rxn[i] <- p$yield[i] * ut - p$dd_rate * y[i]
    # NTS: Death is lost--need to have a substrate pool for it
    # Substrate consumption and product formation 
    rxn[names(st)] <- st * ut + rxn[names(st)]
  }

  # Hydrolysis of particulate substrates and fermentation
  hydrol[isub] <- - alpha[isub] * y[isub]
  # Production of arbitrary products based on specified fermentation stoichiometry (can omit components)
  hydrol[rownames(p$stoich)] <- - p$stoich %*% hydrol[colnames(p$stoich)]
  
   # Add vectors to get derivatives
  # All elements in g/d as COD except 
  #   * slurry_mass (kg/d as fresh slurry mass)
  #   * CH4 (g/d as CH4 or C?)
  #   * solutes other than VFA (...)
  ders <- inflow + rxn + hydrol

  return(list(ders, c(CH4_emis_rate = rxn[['CH4']], temp_C = p$temp_C, pH = p$pH)))

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
