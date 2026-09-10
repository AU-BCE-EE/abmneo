# Info from ABM pig pars 3.0

In

I see

```
arrh_pars_pig3.0 <- list(lnA =  29.96298,
                         E_CH4 = 81000, 
                         A = c(xa_dead= 3.61383 * 10^12, starch = 5 * 3.61383 * 10^12, Cfat = 0.1 * 3.61383 * 10^12, CPs = 3.61383 * 10^12, CPf = 3.61383 * 10^12, RFd = 3.61383 * 10^12, VSd = 3.61383 * 10^12, urea = 4.38*10^15), 
                         E = c(xa_dead= 81557, starch = 81557, Cfat = 81557, CPs = 81557, CPf = 81557, RFd = 81557, VSd = 81557, urea = 81559),  
                         R = 8.314,  
                         VS_CH4 = 6.67,
                         scale_alpha_opt = list(VSd = 0.46056, notVSd = 0.556886, CP = 2.065),
                         kl = c(NH3 = 336, NH3_floor = 0.36*337, H2S = 0.02),
                         scale_EF_NH3 = 0.0000048)
```

Applying these . . .

```
arrhenius <- function(temp_K, A, E, arr_max_temp_K = 313, R = 8.314){
  
  if (temp_K >= arr_max_temp_K) temp_K <- arr_max_temp_K

  y <- A * exp(-E/(R * temp_K))
  
  return(y)
}

0.46056 * arrhenius(273.15 + 20, A = 3.61383e12, E = 81557)
0.46056 * arrhenius(273.15 + 30, A = 3.61383e12, E = 81557)

```

```
> 0.46056 * arrhenius(273.15 + 20, A = 3.61383e12, E = 81557)
[1] 0.004881761
> 0.46056 * arrhenius(273.15 + 30, A = 3.61383e12, E = 81557)
[1] 0.01472195
```

I don't believe those.
Too low!
Lower than digestate pars.

my q10 func

```
q10 <- function(temp_K, yref, qpar, temp_ref = 273.15 + 20, max_temp = 313) {

  temp_K[temp_K >= max_temp] <- max_temp

  y <- yref * qpar^((temp_K - temp_ref) / 10)

  return(y)

}
```

Rearranging

```
y / yref = qpar^((temp_K - temp_ref) / 10)
log10(y / yref) = ((temp_K - temp_ref) / 10) * log10(qpar)

qpar = 10^( log10(y / yref) / ((temp_K - temp_ref) / 10) )
10^( log10(0.0147 / 0.00488) / ((303.15 - 293.15) / 10) )
```

Gives 3.0

