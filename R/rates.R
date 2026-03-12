# rates() and related functions

rates <- function(t, y, parms) {
    
  # Short name for parms to make indexing in code below simpler 
  p <- parms

  # Update any inhibition effects
  ## NTS: yc <- y[intersect(names(y), gsub('_conc', '', rownames(ic0)))] / y['slurry_mass']
  qhat <- update_inhib(p, y)

  # Initialize vectors with derivative components, all with same order of y elements
  inflow <- metab <- hyferm <- 0 * y

  # Inflow from slurry addition
  inflow[names(p$conc_fresh)] <- p$conc_fresh * p$slurry_prod_rate
  inflow[c('slurry_mass', 'slurry_load')] <- 1 * p$slurry_prod_rate
  inflow['COD_load'] <- sum(inflow[names(p$conc_fresh)])

  # Substrate matrix for utilization rate
  sm <- matrix(y[rownames(p$mstoich)], 
               nrow = nrow(p$mstoich), 
               ncol = ncol(p$mstoich), 
               byrow = FALSE)

  # Monod term
  monod <- sm / (sm + y['slurry_mass'] * p$ksmat)
  # Force to 1 for non-substrates
  monod[p$mstoich >= 0] <- 1
  # Product across multiple substrates
  monodprod <- apply(monod, 2, prod)

  # Utilization rate 
  rut <- qhat * monodprod * y[p$grps] / y['slurry_mass']

  # And consumption, growth, production all in one matrix operation
  metab[rownames(p$mstoich)] <- p$mstoich %*% rut * y['slurry_mass']

  # Convert CH4 from g COD / d to g C / d
  metab['CH4'] <- metab['CH4'] / p$COD_conv['CH4']

  # Hydrolysis and fermentation of particulate substrates
  # Consumption and production all in one line, including arbitrary products based on specified fermentation stoichiometry
  hyferm[rownames(p$fstoich)] <- p$fstoich %*% (p$alpha[p$subs] * y[p$subs])

  # Add vectors to get derivatives
  # All elements in g/d as COD except 
  #   * slurry_mass (kg/d as fresh slurry mass)
  #   * CH4 (g/d as CH4-C)
  #   * user-defined solutes other than VFA (as C, N, or S as described in documentation)
  ders <- inflow + metab + hyferm

  return(list(ders, c(CH4_emis_rate = metab[['CH4']], temp_C = p$temp_C, pH = p$pH)))

}

# Copy time-variable parameters into their normal par position for an interval
update_var_pars <- function(series, pars, y, i) {

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

