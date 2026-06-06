# abmneo summary for PIGMET

By Sasha D. Hafner

A summary to go with a video for PIGMET meeting June 2026

# Overview

The ABM package grew out of a project with the two of us and Sven Sommer and Søren O. Petersen that must have originally started more than 10 years ago.
Since the first release that I probably mostly wrote, both Frederik and I have tweaked it but also added a lot of functionality. 
Complexity and size generally increased.
The new abmneo package is an alternative implementation of the ABM model that is substantially simpler.
It came from a ABM branch `simpler1` started in summer 2025.
I moved this work to a new abmneo package in early 2026 (first commit 18 Feb 2026), and further simplified.
At some point in 2026 I started using AI tools for the project, particularly Claude Code, but mainly for improving implementation of my ideas, checking for errors, or sorting out code for some sophisticated approaches. 
This summary explains some of the differences between abmneo and ABM.

# Motivation and progress

The abmneo package was developed to address some issues with ABM:

1. size of codebase, 
2. complexity (some deliberate, some not),
3. hard-coded behavior, special cases, and flexibility
4. number of software dependencies,
5. lack of verification,
6. use of Rcpp, and
7. model speed.

Only some of these are covered in detail in this summary.

Here is what I was able to do:

1. size of the codebase is half as large and the `rates()` function in particular is shorter and clearer,
2. model and codebase are both much simpler,
3. eliminated most hard-coding; even microbial groups can have any metabolic pathways and any names,
4. deSolve is the only dependency,
5. COD balance is built in
6. Rcpp is no longer used, but
7. the `abmneo::abmneo()` function is about as fast as `ABM::abm()` (faster in some cases)

Some of these changes are clear improvements (thinking about how microbial groups are handled now).
Others mean fewer features and processes.
A few are highlighted below.

# Codebase size

ABM still does more than abmneo, but core functionality is similar.
abmneo does it with much quite a bit less code.

|                         | abmneo    | ABM |
|---                      |--:        |--:  |
| R source files          | 9         | 23  |
| C++ source files        | 0         | 5   |
| **Total source files**  | **9**     | **28** |
| Lines of R              | 1,367     | 1,929 |
| Lines of C++            | —         | 643 |
| **Total lines**         | **1,367** | **2,572** |
| **Functions**           | **36**    | **33** |

Total lines of code (excluding auto-generated R code by Rcpp) is about 1400 lines versus about 2600 lines.
The number of functions is higher in abmneo, but the length of function definitions is much less.
Two drastic examples are 

* 104 lines for abmneo's `rates.R` vs. 563 for `rates_cpp.cpp`, and
* 191 lines for the core abmneo functions in `core.R` vs. 572 lines for `abm.R`.

The abmneo `rates()` function uses vectorized and matrix operations to reduce code complexity.
Shorter function definitions make understanding, maintaining, improving, extending, and debugging all easier, and reduce the risk of errors.

# Complexity, hard-coding, and flexibility

I tried to keep complexity down during abmneo development, because I struggle to understand some of the ABM code. 
I also tried to move decisions from the package code to input parameters.
In abmneo the reduction in complexity provided an increase in flexibility as well.
One extreme example is on substrate definitions and conversion in the model.
I worked hard on this!

In ABM, there are hard-coded substrate groups that are taken through multiple steps and somewhat convoluted steps for conversion. 
With C++ code, position is used for indexing.  
Adding or changing substrates means writing C++ code.
In abmneo, there are no hard-coded substrate groups.
Substrate are completely flexible and specified in parameters provided by the user.
Model code is based on a small number of matrix operations, with indexing automatically done by substrate name.
And stoichiometry is as simple as a 1:1 relationship betwee particulate substrates and VFAs, although there is an option to specify a stoichiometry matrix as an input parameter when other products are important, like NH3 or H2S.
Default parameter sets provide a way to replicate ABM behavior in abmneo (without effort from the user).

Microbial groups are another example.
In ABM, we had hard-coded group types: methanogens and sulfate reducers.
We used the names to figure out which was which, and have different code for the two.
In C++ this means using positional indexing.
In abmneo, groups types are defined in input parameters, and could include any stoichiometry, with flexible Monod kinetics as well.
The code relies on efficient matrix operations like substrate processing.
This doesn't make things challenging for users because normal microbial groups can be defined in default input parameter sets.

# A few other improvements
On verification, `abmneo()` automatically includes a check of the COD balance at the end of each simulation.
In theory it would catch problems with code and even some input parameters.

Rcpp is a package that support C++ code, which is compiled to run much faster than R code.
I was initially excited about it, but it makes the code much more complex and difficult to debug.
And adds dependencies for users.
It is not used in abmneo.

Other dependencies were eliminated as well by using base R functionality.
I am always happy to avoid tidyverse packages, because I have an unfair bias against them, sure, but also because they are a mess of dependencies and more prone to significant changes than most R packages.

# Missing functionality
abmneo does less than ABM.
One major design difference is in how time-variable inputs are used.
In ABM, several inputs variables, like temperature and pH, can vary smoothly over time.
Internally, this is done by creating interpolation functions that are then called from within the `rates()` function (or `rates_cpp()`)).
The approach is available to only some variables, and relies on hard-coding.
In general, this smooth approach slows down model runs, because more functions need to be called every time `rates()` is called.
Rcpp helps with this issue a lot!

In abmneo, time-variable inputs change in steps.
The user decides on the temporal resolution.
Output is less pleasing to look at, but methane emission predictions are nearly the same. 
And it is possible to change more input parameters--almost any input can vary over time.

Frederik had to pressure me a bit to include inhibition and respiration, mainly because they add complexity.
But now they are in and programmed efficiently I think (especially inhibition) without special cases.
I am sure there are other ABM processes not included; Frederik can list them.

# Comparison
Take a look at the two `comp` vignettes.
