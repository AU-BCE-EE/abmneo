# Notes on par est

Fitting to Kasper AD tank data.
From help file in this package:

  A single set "dig_01" is presently available, and was fit to methane emission measurements from a set of three digestate storage tanks in Denmark.
  Measurements are described in Kasper et al. (2026) and original measurements were shared by first author P. Kasper.
  Note that the four objects are meant to be used together; the `resid_enrich` value in \code{storage_ex_dig_01} was fit to measurements along with the hydrolysis constant \code{h_rate_ref}, while other parameters were loosely based on ABM parameter sets or "expert judgement".
  
  Kasper, P.L., Lycus, P.M. & Holm, M., 2026. Field-scale assessment of ammonia and methane emissions from tent-covered slurry storages using high-resolution tracer decay measurements. Biosystems Engineering, 270, p. 104551. \doi{10.3168/jds.2012-5340}.

Using ABM as inspiration or quantitative values.

From `make-Rda/make_pars.R` at 3a4d47d1b69272a158bd52e5939477e245d0fb14

```
arrh_pars_digestate3.0 <- list(lnA = 30.33312,
                               E_CH4 = 81000, 
                               A = c(xa_dead= 22497134.8, starch = 5.86*10^18, Cfat = 0, CPs = 1.763612 * 10^12, CPf = 4.10669 * 10^11, RFd = 2.567944 * 10^8, VSd = 22497134.8, urea = 4.38*10^15), 
                               E = c(xa_dead= 55609, starch = 109400, Cfat = 0, CPs = 87515.5450, CPf = 69361.3, RFd = 57867, VSd = 55609, urea = 81559),  # CP is for GRA treatment in MILK project
                               R = 8.314,  
                               VS_CH4 = 6.67,
                               scale_alpha_opt = list(VSd = 2.2829985, notVSd = 0.23923, CP = 1.0255),
                               kl = c(NH3 = 1000, NH3_floor = 0.36*1000, H2S = 0.02),
                               scale_EF_NH3 = 0.0000048)
```


So, what is reference rate constants?
ABM 3a4d47d1b69272a158bd52e5939477e245d0fb14 has:

```
arrh_pars_digestate3.0 <- list(lnA = 30.33312,
                               E_CH4 = 81000, 
                               A = c(xa_dead= 22497134.8, starch = 5.86*10^18, Cfat = 0, CPs = 1.763612 * 10^12, CPf = 4.10669 * 10^11, RFd = 2.567944 * 10^8, VSd = 22497134.8, urea = 4.38*10^15), 
                               E = c(xa_dead= 55609, starch = 109400, Cfat = 0, CPs = 87515.5450, CPf = 69361.3, RFd = 57867, VSd = 55609, urea = 81559),  # CP is for GRA treatment in MILK project
                               R = 8.314,  
                               VS_CH4 = 6.67,
                               scale_alpha_opt = list(VSd = 2.2829985, notVSd = 0.23923, CP = 1.0255),
                               kl = c(NH3 = 1000, NH3_floor = 0.36*1000, H2S = 0.02),
                               scale_EF_NH3 = 0.0000048)

```

And don't forget

```
arrh_pars_digestate3.0 <- list(lnA = 30.33312,
                               E_CH4 = 81000, 
                               A = c(xa_dead= 22497134.8, starch = 5.86*10^18, Cfat = 0, CPs = 1.763612 * 10^12, CPf = 4.10669 * 10^11, RFd = 2.567944 * 10^8, VSd = 22497134.8, urea = 4.38*10^15), 
                               E = c(xa_dead= 55609, starch = 109400, Cfat = 0, CPs = 87515.5450, CPf = 69361.3, RFd = 57867, VSd = 55609, urea = 81559),  # CP is for GRA treatment in MILK project
                               R = 8.314,  
                               VS_CH4 = 6.67,
                               scale_alpha_opt = list(VSd = 2.2829985,

```

Applying these . . .

```
arrhenius <- function(temp_K, A, E, arr_max_temp_K = 313, R = 8.314){
  
  if (temp_K >= arr_max_temp_K) temp_K <- arr_max_temp_K

  y <- A * exp(-E/(R * temp_K))
  
  return(y)
}

2.283 * arrhenius(273.15 + 20, A = 22497134.8, E = 55609)
2.283 * arrhenius(273.15 + 30, A = 22497134.8, E = 55609)

```

That gives 0.00633 1/d for 20 C.

And 30 C gives 0.01344 1/d.

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
10^( log10(0.0134 / 0.00633) / ((303.15 - 293.15) / 10) )
```

Gives 2.122

Check

```
q10(273.15 + 20, 0.00633, 2.122)
q10(273.15 + 30, 0.00633, 2.122)
```

```
> q10(273.15 + 20, 0.00633, 2.122)
[1] 0.00633
> q10(273.15 + 30, 0.00633, 2.122)
[1] 0.01343226
```

Looks good.


