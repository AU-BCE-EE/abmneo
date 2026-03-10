# Assessment of ABM Issues

*Compiled by Claude (Anthropic)*

This document evaluates the list of issues with the ABM package identified in `abmneo_demo.Rmd`, with specific examples from the ABM source code.

## The list

1. complexity (some deliberate, some not)
2. hard-coded behavior and special cases
3. number of dependencies
4. lack of verification
5. size of codebase
6. use of Rcpp
7. model speed

## Assessment

### 1. Complexity — Fair

`rates_cpp.cpp` alone is 563 lines doing many things in sequence: extracting 27+ state variables, computing 5 inhibition types, two inhibition modes (overrule vs. multiplicative), emission rates, VFA utilization for two separate microbial metabolisms, fermentation stoichiometry, and derivatives for all 27 state variables. It is genuinely hard to hold in your head.

### 2. Hard-coded behavior and special cases — Fair, and arguably understated

The clearest examples:

- `hard_pars.R` lines 7–8: microbial groups are classified by *name prefix* — `grep('^[mp]', ...)` identifies methanogens, `grep('^sr', ...)` identifies sulfate reducers. These indices then control completely different code paths in `rates_cpp.cpp` (lines 260–296). If you named a group `"meth2"` thinking it was a second methanogen, it would work; if you named it `"acetoclastic"`, it would silently be treated as neither.
- The state vector in `rates_cpp.cpp` has a fixed layout with hard-coded names for all substrates (RFd, iNDF, ash, VSd, starch, CPs, CPf, Cfat, VFA, urea, TAN, sulfate, sulfide, ...). You cannot add or remove a substrate without editing the C++ source.
- `VSd` is special-cased for fermentation stoichiometry and Arrhenius scaling (checked explicitly in `abm.R` around line 290).
- pH threshold of 6.8 for H2S inhibition IC50 switching (`rates_cpp.cpp` line 220).
- HAC threshold of 0.05 g/kg switching inhibition model (line 215).
- Respiration requires depth ≥ 1 m, hard-coded (line 256).

### 3. Dependencies — Mostly fair, with one omission

DESCRIPTION lists deSolve, Rcpp, and dplyr — but `stoich_species.R` also imports `biogas`, which is not mentioned. Minor, but worth knowing.

### 4. Lack of verification — Fair

There is a `dev-tests/` directory with some scripts, but no formal testing framework (no `testthat`). More tellingly, `abm.R` has a large block of commented-out COD/VS/C flow calculations (lines ~311–416) with a note literally reading `"NTS: These need to be deleted"` — suggesting mass balance tracking was started and abandoned rather than completed.

### 5. Size of codebase — Fair

~1,950 R lines + ~738 C++ lines = ~2,688 total across 30 files.

### 6. Rcpp — Fair

Five C++ files totaling ~738 lines, with `rates_cpp.cpp` being 563 of those. The compilation requirement, debugging difficulty, and install complexity are real costs.

### 7. Speed — Fair as stated

The vignette is appropriately modest ("hope to not make it worse"). There is a `dev-examples/speed_test_dat/` directory, suggesting speed was actively tracked.

## One issue not on the list

`stoich.R` line 24 has a live `browser()` call — a debugging breakpoint that would pause any simulation hitting that code path and drop the user into an interactive debugger. This is a concrete fragility issue that arguably falls under point 2 but is specific enough to call out on its own.

## Overall

The list is fair and well-grounded. The main understatement is point 2 — the hard-coding runs deeper than just methanogens vs. sulfate reducers; the fixed state vector layout means the substrate list itself is hard-coded into the C++.
