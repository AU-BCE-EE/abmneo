# Stoichometry functions

# Check COD balance
check_COD <- function(dat, pars, rtol = 0.001) {

  first <- unlist(dat[1, ])
  last <- unlist(dat[nrow(dat), ])

  supercomps <- pars$supercomps
  gases <- pars$gases
  COD_conv <- pars$COD_conv

  # Present at start
  CODstart <- sum(first[supercomps] * COD_conv[supercomps]) 
  # Loading
  CODin <- last[['COD_load']] 
  # Removed in effluent
  CODeff <- sum(last[paste0(supercomps, '_eff')] * COD_conv[supercomps])
  # Emitted
  CODemis <- last[[gases]] * COD_conv[[gases]]
  # Remaining
  CODrem <- sum(last[supercomps] * COD_conv[supercomps]) 

  bal <- CODstart + CODin - CODeff - CODemis - CODrem
  rbal <- bal / (CODstart + CODin)

  if (abs(rbal) > rtol) {
    warning('COD balance is off by ', signif(100 * rbal, 2), '%')
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

  # For convenience
  msc <- mstoich

  # R&M's fe
  femat <- matrix(1 - yield, nrow = nrow(msc), ncol = length(yield), byrow = TRUE)

  # Subtract microbial yields from product formation
  # Note that this is not simply subtraction because some products may not be in COD units!
  msc[msc > 0] <- (msc * femat)[msc > 0]

  # Subtract microbial yields from secondary reactant (e.g., electron acceptor)
  mscnovfa <- msc[rownames(msc) != 'VFA', ]
  fematnovfa <- femat[rownames(msc) != 'VFA', ]
  mscnovfa[mscnovfa < 0] <- (mscnovfa * fematnovfa)[mscnovfa < 0]

  msc[rownames(msc) != 'VFA', ] <- mscnovfa

  # Yield matrix
  ym <- diag(yield, ncol(msc), ncol(msc))
  rownames(ym) <- names(yield)
  combstoich <- rbind(msc, ym)

  return(combstoich)

}

fix_ksmat <- function(ksmat, mstoich) {

  ks <- 0 * mstoich + 1e10
  ks[rownames(ksmat), colnames(ksmat)] <- ksmat

  return(ks)

}

make_dstoich <- function(pars) {

  sm <- diag(-1, length(pars$grps), length(pars$grps))
  dstoich <- rbind(sm, rep(1, length(pars$grps)))
  rownames(dstoich) <- c(pars$grps, pars$xd)
  colnames(dstoich) <- pars$grps

  return(dstoich)

}
