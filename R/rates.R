# rates() and related functions

rates <- function(t, y, parms) {
    
  # Short name for parms to make indexing in code below simpler 
  p <- parms

  # Update any inhibition effects
  ## NTS: yc <- y[intersect(names(y), gsub('_conc', '', rownames(ic0)))] / y['slurry_mass']
  qhat <- update_inhib(p, y)

  # Initialize vectors with derivative components, all with same order of y elements
  inflow <- metab <- death <- hyferm <- resp <- 0 * y

  # Inflow from slurry addition
  inflow[names(p$conc_fresh)] <- p$conc_fresh * p$slurry_prod_rate
  inflow[c('slurry_mass', 'slurry_load')] <- 1 * p$slurry_prod_rate
  inflow['COD_load'] <- sum(inflow[p$supercomps] * p$COD_conv[p$supercomps])

  # Substrate matrix for utilization rate
  sm <- matrix(y[rownames(p$mstoich)], 
               nrow = nrow(p$mstoich), 
               ncol = ncol(p$mstoich), 
               byrow = FALSE)

  # Names for debugging
  #dimnames(sm) <- dimnames(p$mstoich)

  # Monod term
  monod <- sm / (sm + y['slurry_mass'] * p$ksmat_t)
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

  # Biomass death
  death[rownames(p$dstoich)] <- p$dstoich %*% (p$d_rate * y[p$grps])

  # Hydrolysis and fermentation of particulate substrates
  # Consumption and production all in one line, including arbitrary products based on specified fermentation stoichiometry
  hyferm[rownames(p$fstoich)] <- p$fstoich %*% (p$h_rate[p$subs] * y[p$subs])

  # Surface respiration: aerobic oxidation limited by O2 surface flux
  # Monod term on VFA prevents negative values at low concentrations (ks = 0.05 g COD/m3)
  # O2 reduces (net) COD loading
  if (p$has_resp) {
    resp_rate <- p$O2_flux * p$area * y['VFA'] / (y['VFA'] + 0.05 * y['slurry_mass'])
    resp[names(p$rstoich)] <- resp_rate * p$rstoich
    resp['COD_load'] <- - resp_rate 
  }

  # Add vectors to get derivatives
  # All elements in g/d as COD except
  #   * slurry_mass (kg/d as fresh slurry mass)
  #   * CH4 (g/d as CH4-C)
  #   * user-defined solutes other than VFA (as C, N, or S as described in documentation)

  ders <- inflow + metab + death + hyferm + resp

  ## Debugging code below
  #if (y['VFA'] < 0) {
  #    print(xx <- c(t = t, y = y, der = ders))
  #    xxx <<- rbind(xxx, xx)
  #    if (t > 6) browser()
  #}

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
    if (length(newval) > 1 & !is.null(names(newval))) {
      pars[[pn]][names(newval)] <- newval # For things like grp_pars, where only some may be supplied
    } else {
      pars[[pn]] <- newval
    }
  }


  # Temperature and temperature-dependent pars could change
  # Convert temperature to K in case any C values were given in var
  pars <- tempsC2K(pars, cutoff = 200)

  # Calculate temperature-dependent par values
  pars <- calc_temp_pars(pars, y)

  return(pars)

}

