# Codebase Size Comparison: abmneo vs. ABM

### abmneo (R only)

| File | Lines | Functions |
|------|------:|----------:|
| `core.R` | 191 | 2 |
| `misc.R` | 96 | 6 |
| `pars.R` | 400 | 9 |
| `post.R` | 168 | 6 |
| `prep.R` | 28 | 1 |
| `rates.R` | 104 | 2 |
| `startup.R` | 45 | 1 |
| `stoich.R` | 95 | 5 |
| `timing.R` | 240 | 4 |
| **Total** | **1,367** | **36** |

### ABM (R + C++)

Auto-generated Rcpp files (`RcppExports.R`, `RcppExports.cpp`) excluded.

| File | Lines | Functions |
|------|------:|----------:|
| `abm.R` | 572 | 1 |
| `abm_regular.R` | 148 | 1 |
| `abm_variable.R` | 262 | 1 |
| `Arrh_func.R` | 5 | 1 |
| `checkGrpNames.R` | 13 | 1 |
| `coverfun.R` | 17 | 1 |
| `CTM.R` | 60 | 1 |
| `doy.R` | 24 | 1 |
| `emptyStore.R` | 50 | 1 |
| `et.R` | 47 | 3 |
| `graze_fun.R` | 15 | 1 |
| `H2SO4_titrat.R` | 31 | 1 |
| `hard_pars.R` | 99 | 1 |
| `logistic.R` | 3 | 2 |
| `makeConcFunc.R` | 26 | 1 |
| `makeTimeFunc.R` | 20 | 1 |
| `makeXaFreshFunc.R` | 28 | 1 |
| `pars_indices.R` | 19 | 1 |
| `predFerm.R` | 199 | 4 |
| `readFormula.R` | 81 | 1 |
| `stoich.R` | 86 | 1 |
| `stoich_species.R` | 114 | 1 |
| `tempsC2K.R` | 10 | 1 |
| **R subtotal** | **1,929** | **28** |
| `Arrh_func_cpp.cpp` | 16 | 1 |
| `call_int_cpp.cpp` | 12 | 1 |
| `CTM_cpp.cpp` | 42 | 1 |
| `extract_xa_cpp.cpp` | 10 | 1 |
| `rates_cpp.cpp` | 563 | 1 |
| **C++ subtotal** | **643** | **5** |
| **Total** | **2,572** | **33** |

### Summary

| | abmneo | ABM |
|---|--:|--:|
| R source files | 9 | 23 |
| C++ source files | 0 | 5 |
| **Total source files** | **9** | **28** |
| Lines of R | 1,367 | 1,929 |
| Lines of C++ | — | 643 |
| **Total lines** | **1,367** | **2,572** |
| **Functions** | **36** | **33** |
