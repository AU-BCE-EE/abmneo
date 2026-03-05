# Support functions (for dev)

# Figure out stoichiometry matrix of fermentation substrates from elemental formula
get_fstoich <- function(subs, convert = TRUE) {

  # Get molar stoichiometric coefficients
  # Vectorize, return matrix without substrate (for 1 mole substrate)
  res <- lapply(as.list(subs), predFerm)
  # Align names and sort before combining in matrix
  nn <- unique(unlist(lapply(res, names)))
  for (i in 1:length(res)) {
    res[[i]][nn[!nn %in% names(res[[i]])]] <- 0
    res[[i]] <- res[[i]][nn]
  }
  st <- matrix(unlist(res), ncol = length(subs), byrow = FALSE)
  rownames(st) <- names(res[[1]])
  colnames(st) <- subs

  # Drop 0
  st <- st[rowSums(st) != 0, , drop = FALSE]
  
  # Drop water (ignored, treated as conservative in system)
  st <- st[rownames(st) != 'H2O', , drop = FALSE] 

  # Adjust coefficients to COD mass, N mass, C mass, S mass, or total mass
  if (convert) {
    for (i in 1:nrow(st)) {
      ff <- rownames(st)[i]
      st[i, ] <- st[i, ] * get_mass_conv(ff)
    }
    
    for (i in 1:ncol(st)) {
      ff <- colnames(st)[i]
      st[, i] <- st[, i] * 1 /  get_mass_conv(ff)
    }
  }

  return(st)

}


