# Substrate Complexity Comparison: ABM vs. abmneo

## Overview

Both models track substrates (organic matter pools) as state variables and use stoichiometry
to convert hydrolysis/fermentation into products consumed by methanogens and sulfate reducers.
The paths from user input to `rates()` are very different in complexity.

---

## ABM: substrate flow

### Step 1 — User provides named concentrations

The user fills in a fixed set of named slots in `conc_fresh`. The names are hardcoded into the
model; using any other name silently does nothing:

```
conc_fresh$starch    # rapidly fermentable carbohydrate
conc_fresh$RFd       # slowly degradable fiber (raw fiber degradable)
conc_fresh$CPs       # soluble crude protein
conc_fresh$CPf       # fiber-bound crude protein
conc_fresh$Cfat      # crude fat
conc_fresh$VSd       # generic volatile solids (alternative to the above five)
```

There is a hard constraint: `VSd` and the other five substrates are mutually exclusive
(`abm.R:292`). This is enforced with a `stop()`.

### Step 2 — Two parallel stoichiometry paths (one unused)

ABM contains a full Rittmann-McCarty stoichiometry pipeline:

```
readFormula()          # readFormula.R  — regex parser for chemical formulas
  └─ customOrgStoich() # predFerm.R     — derives half-reaction from element counts (R&M O-19)
       └─ predFerm()   # predFerm.R     — combines donor/acceptor/synthesis half-reactions
            └─ stoich_species()  # stoich_species.R — normalizes to gCOD, builds full matrix
```

`stoich_species()` produces a stoichiometry matrix (state variables × compounds) from the
chemical formulas `'C6H10O5'`, `'C51H98O6'`, `'C4H6.1O1.2N'`, etc.

**However, `stoich_species()` is never called in the main simulation.** Instead, `hard_pars()`
hardcodes the result as named numeric vectors with literal coefficients (`hard_pars.R:38-63`):

```r
pars$carb <- c(C6H10O5 = -1, C51H98O2 = 0, C4H6.1O1.2N = 0,
               NH3 = 0, H2O = -2.992714, C5H7O2N = 0,
               C2H4O2 = 2.003743, H2 = 3.985329, CO2 = 1.993114)
pars$pro  <- c(...)
pars$lip  <- c(...)
pars$OM   <- c(...)   # VSd path
```

There are also separate stoichiometry vectors for methanogenesis and sulfate reduction
(`pars$ace`, `pars$hyd`, `pars$ace_sr`, `pars$hyd_sr`) and for respiration
(`pars$carb_resp`, `pars$pro_resp`, `pars$lip_resp`). Each is a named numeric vector.
Commented-out alternatives ("with cell synthesis") remain in the file alongside the
active values, adding to the visual complexity.

### Step 3 — Substrates grouped by biochemical class in `rates_cpp.cpp`

Inside the C++ rates function, substrates are collapsed into three scalar mole-flows using
hardcoded molar-mass conversion factors and hardcoded position indices
(`rates_cpp.cpp:120–172`, `361–378`):

```cpp
double alpha_RFd   = alpha[5];   // position hard-coded
double alpha_starch = alpha[1];  // position hard-coded
// ...

double mol_carb = (alpha_RFd * RFd + alpha_starch * starch) * 0.005208333;
double mol_pro  = (alpha_CPf * CPf + alpha_CPs  * CPs ) * 0.00748503;
double mol_lip  =  alpha_Cfat * Cfat             * 0.0004194631;

for (int i = 0; i < carb.size(); ++i)
  ferm[i] = mol_carb * carb[i] + mol_pro * pro[i] + mol_lip * lip[i];
```

The molar mass factors (`0.005208333` = mol/gCOD for starch, etc.) are inline magic numbers
with no names. Adding a new substrate class would require modifying the C++ source.

### Step 4 — Per-class product extraction and secondary stoichiometry

From `ferm`, acetate (`C2H4O2`) and H₂ are extracted by position (`ferm[6]`, `ferm[7]`),
then scaled against separate stoichiometry vectors for methanogenesis and sulfate reduction.
This yields four more mole-flow scalars (`ace_ferm`, `hyd_ferm`, `ace_sr_ferm`, `hyd_sr_ferm`).

### Step 5 — State variable derivatives by hardcoded position

Each derivative is assigned to a hardcoded integer offset (`n_mic + k`):

```cpp
derivatives[n_mic+2]  = slurry_prod_rate * conc_fresh_RFd - alpha_RFd * RFd ...
derivatives[n_mic+6]  = slurry_prod_rate * conc_fresh_starch - alpha_starch * starch ...
derivatives[n_mic+10] = alpha_xa_dead * xa_dead + VFA_H2_ferm - sum_rut + ...
```

Adding, removing, or renaming a substrate requires updating the C++, rebuilding the package,
and adjusting any downstream index offsets.

### Summary of ABM substrate path

```
User: conc_fresh$starch, conc_fresh$RFd, ...  (6 fixed names)
         |
         |   [formula-parsing pipeline exists but is unused in simulation]
         |
hard_pars(): pars$carb, pars$pro, pars$lip  (hardcoded literal coefficients)
         |
rates_cpp.cpp:
  alpha[] indexed by position
  mol_carb / mol_pro / mol_lip  (hardcoded molar mass factors)
  ferm[] = weighted sum over classes
  derivatives[n_mic+k] = ...  (hardcoded integer offsets)
```

---

## abmneo: substrate flow

### Step 1 — User provides names and stoichiometry

The user defines substrates freely. There are no fixed names. Stoichiometry is supplied
directly as a matrix in gCOD/gCOD units — no chemical formulas or molar masses needed:

```r
sub_pars <- list(
  subs     = c('VSd'),
  mstoich  = matrix(c(-1, 1), nrow = 2,           # rows = reactants/products, cols = microbial groups
                    dimnames = list(c('VFA', 'VSd'), c('m0'))),
  sub_fresh = c(VSd = 30),
  sub_init  = c(VSd = 30),
  hydrol_opt = c(VSd = 0.15),
  T_opt_hyd  = c(VSd = 35), ...
)
```

`fstoich` (hydrolysis product stoichiometry) is optional; the default is that all substrates
hydrolyse to VFA.

### Step 2 — `pack_pars()` assembles the stoichiometry matrices once

`pack_pars()` (`pars.R`) does four short matrix operations at setup time:

```r
pars$fstoich <- comb_fstoich(pars$fstoich, pars$subs)  # add substrate self-consumption rows
pars$mstoich <- comb_mstoich(pars$mstoich, pars$yield) # add yield rows, scale fe into stoich
pars$dstoich <- make_dstoich(pars)                     # death: biomass → xd substrate
pars$ksmat   <- fix_ksmat(pars$ksmat, pars$mstoich)    # align ks dimensions with mstoich
```

Each of these is 5–15 lines of pure matrix arithmetic (`stoich.R`). The result is three
matrices — `mstoich`, `fstoich`, `dstoich` — stored in `pars` and reused on every ODE step.

### Step 3 — `calc_temp_pars()` computes temperature-adjusted rates

```r
pars$alpha[pars$subs] <- arrhenius(pars$temp_K, pars$arrA, pars$arrE, pars$arr_max_temp_K)
```

One vectorized call over all substrates. No class grouping, no magic number conversion factors.

### Step 4 — `rates()` uses three matrix multiplications

All metabolism, death, and hydrolysis/fermentation are each one line (`rates.R:40–50`):

```r
metab[rownames(p$mstoich)]  <- p$mstoich  %*% rut               * y['slurry_mass']
death[rownames(p$dstoich)]  <- p$dstoich  %*% (p$d_rate * y[p$grps])
hyferm[rownames(p$fstoich)] <- p$fstoich  %*% (p$alpha[p$subs] * y[p$subs])
```

All indexing is by name. Adding a new substrate means adding a row/column to the user-supplied
`mstoich`; no model code changes.

### Summary of abmneo substrate path

```
User: subs = c('VSd'), mstoich = matrix(...), sub_fresh = c(VSd = 30)
         |
pack_pars():
  comb_fstoich / comb_mstoich / make_dstoich  (matrix assembly, ~30 lines total)
         |
calc_temp_pars():
  pars$alpha[subs] <- arrhenius(...)  (one vectorized call)
         |
rates():
  metab  <- mstoich  %*% rut
  death  <- dstoich  %*% d_rate * xa
  hyferm <- fstoich  %*% alpha * y[subs]
```

---

## Side-by-side comparison

| Aspect | ABM | abmneo |
|---|---|---|
| Substrate names | 6 fixed (starch, RFd, CPs, CPf, Cfat, VSd) | User-defined, arbitrary |
| Stoichiometry source | Hardcoded literal vectors in `hard_pars.R` | User-supplied matrix in gCOD/gCOD |
| Formula-parsing pipeline | Present (`readFormula` → `predFerm` → `stoich_species`, ~370 lines) | Not needed |
| Formula-parsing used? | No — hardcoded values bypass it | — |
| Biochemical class grouping | carb / protein / lipid (hardcoded molar mass factors) | None — all substrates treated identically |
| Substrate indexing in rates | Integer position offsets (`alpha[5]`, `derivatives[n_mic+2]`) | Named (`pars$alpha['VSd']`, `y['VSd']`) |
| Language of rates function | C++ (563 lines) | R (77 lines) |
| Lines to add a new substrate | Must edit C++, rebuild, adjust offsets | Add row/column to `mstoich`; no code change |
| Stoichiometry matrices | Four separate named vectors per class | Three generic matrices (`mstoich`, `fstoich`, `dstoich`) |
| Key rates lines | ~50 lines of scalar arithmetic | 3 matrix multiplications |
