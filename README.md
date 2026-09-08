# abmneo
R model for anaerobic microbial degradation (biodegradation) of organic matter with multiple microbial groups

# Description
The ABM model predicts conversion of animal manure or other high-moisture organic wastes to methane and carbon dioxide under anaerobic conditions.
The name comes from **a**naerobic **b**iodegradation **m**odel. 
With multiple microbial groups and group-specific parameters describing kinetics and yield, the model can predict realistic short- and long-term responses to temperature change.
ABM was originally implemented in the ABM R package, still maintained here: <https://github.com/AU-BCE-EE/ABM>.
The abmneo package is a new, simplified but more flexible version, with less hard-coding and fewer special cases in the code.

# Maintainers
Sasha D. Hafner (<https://au.dk/sasha.hafner@bce.au.dk>)

# Installation and getting started
This R package can be installed using `remotes::install_github("AU-BCE-EE/abmneo", build_vignettes = TRUE)`.
Load the package with `library(abmneo)`, and open up the vignette with `vignette("abmneo_start")`.

# Problems?
Report problems or ask questions using the Issues page.

# Research papers
See these papers for more information on the model, as well as examples of applications of the ABM package.

Dalby, F.R., Hafner, S.D., Petersen, S.O., Vanderzaag, A., Habtewold, H., Dunfield, K., Chantigny, M.H., Sommer, S.G., 2021. A mechanistic model of methane emission from animal slurry with a focus on microbial groups. PLOS ONE.  <https://doi.org/10.1371/journal.pone.0252881> ([download here](https://drive.google.com/file/d/17ECrisyZXCL2x7hV647ql2hG3hMMGmsD/view?usp=sharing))

Dalby, F.R., Hansen, M.J., Guldberg, L.B., Hafner, S.D., Feilberg, A., 2023. Simple Management Changes Drastically Reduce Pig House Methane Emission in Combined Experimental and Modeling Study. Environ. Sci. Technol. <https://doi.org/10.1021/acs.est.2c08891> ([download here](https://drive.google.com/file/d/1FsFzUrbM2O3GlTm3TYHjrlTfr_mFfLRC/view?usp=share_link))

Dalby, F.R., Ambrose, H.W.,  Poulsen, J.S., Nielsen, J.L., Adamsen, A.P.S., 2023. Pig slurry organic matter transformation and methanogenesis at ambient storage temperatures. JEQ. <https://doi.org/10.1002/jeq2.20512> ([download here](https://drive.google.com/file/d/1636ODKPMHjlmNP50sARJViEWv0JBPfRP/view?usp=sharing))

Dalby, F.R., Hafner, S.D., Ambrose, H.W., Adamsen, A.P.S., 2024. Pig manure degradation and carbon emission: Measuring and modeling combined aerobic-anaerobic transformations. JEQ. <https://doi.org/10.1002/jeq2.20603> ([download here](https://drive.google.com/file/d/1m56KUeJ5Rj9pN6ZeSjMQdNzgFENmZJ_l/view?usp=sharing))
