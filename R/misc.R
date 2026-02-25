
# Replace missing values in multiple data frame columns with interpolated values
interpm <- function(dat, x, ys, by = NA, ...) {

  if (is.na(by)) {
    for (i in ys) {
      rout <- which(is.na(dat[[i]])) 
      if (length(rout) > 0) {
        dat[[i]][rout] <- approx(dat[[x]][-rout], dat[[i]][-rout], xout = dat[[x]][rout], ...)$y 
      } 
    }
  } else {
    for (i in ys) {
      for (j in unique(dat[[by]])) {
        gr <- dat[[by]] == j
        rout <- is.na(dat[[i]])
        if (length(rout) > 0) {
          dat[[i]][gr & rout] <- approx(dat[[x]][gr & !rout], dat[[i]][gr & !rout], xout = dat[[x]][gr & rout], ...)$y 
        } 
      }
    }

  }

  return(dat)

}

logistic <- function(x) exp(x)/(1 + exp(x))

logit <- function(p) log(p/(1 - p))

# Fill values down, replacing missing values
# Suggestion from ChatGPT with a tweak (x[cumsum... -> x[idx][cumsum...)
fill_down <- function(x) {
  if (!anyNA(x)) return(x)
  
  idx <- !is.na(x)
  if (!any(idx)) return(x)  # all NA column
  
  x[!idx] <- x[idx][cumsum(idx)][!idx]

  return(x)
}

# When a column is a list of named vectors, e.g., for qhat_opt or even sub_fresh (could be length-one vectors)
# Another ChatGPT suggestion with tweaks
fill_down_list <- function(x) {
  last <- NULL

  for (i in seq_along(x)) {
    if (any(is.na(x[[i]]))) {
      x[[i]] <- last
    } else {
      last <- x[[i]]
    }
  }

  return(x)
}

fill_down_df <- function(df) {
  for (nm in names(df)) {
    if (is.atomic(df[[nm]])) {
      df[[nm]] <- fill_down(df[[nm]])
    } else {
      df[[nm]] <- fill_down_list(df[[nm]])
    }
  }

  return(df)

}

#fill_down_df <- function(df) {
#
#browser()
#i = 10
#x1 = df[, 1]
#class(x1)
#x = df[, i]
#x
#length(x)
#  x[[idx]][[cumsum(idx)]][[!idx]]
#
#class(x[1])
#y = x[idx][cumsum(idx)][!idx]
#class(y[1])
#y[1]
#class(x)
#  for (i in seq_len(ncol(df))) {
#    df[, i] <- fill_down(df[, i])
#  }
#
#  return(df)
#}
