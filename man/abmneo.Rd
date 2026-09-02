\name{abmneo}
\alias{abmneo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  An ABM model
}
\description{
  An implementation of the ABM model for simulation of anaerobic biodegradation of organic substrates.
}
\usage{
abmneo(
  storage,
  days = 365,
  delta_t = 1,
  times = NULL,
  inf_pars = NULL,
  grp_pars = NULL,
  sub_pars = NULL,
  chem_pars = list(
    COD_conv = c(CH4 = 5.32, SO4 = -1.99, CO2 = 0, H2S = 0.498), 
    gases = c('CH4', 'CO2', 'H2S')
  ),
  ctrl_pars = list(
    rem_align_method = 'mid',
    var_fill_method = 'interp',
    par_key = '\\.',
    h_rate_max_temp = 313
  ),
  var_pars = list(var = NULL),
  add_pars = NULL,
  pars = NULL,
  startup = 0,
  starting = NULL,
  warn = TRUE,
  quiet = FALSE
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{storage}{
    A storage object.
  }
  \item{days}{
    Total number of days in simulation.
  }
  \item{delta_t}{
    Time step in output (d).
  }
  \item{times}{
    Alternative to \code{delta_t}, use to return specific times only.
    Only available with variable slurry mass.
  }
  \item{inf_pars}{
    Influent parameters.
  }
  \item{grp_pars}{
    Microbial group parameters.
  }
  \item{sub_pars}{
    Substrate parameters.
  }
  \item{chem_pars}{
    Chemical parameters.
    See default values.
  }
  \item{ctrl_pars}{
    Control parameters.
    See default values.
  }
  \item{var_pars}{
    Optional.
    Time-varying inputs.
  }
  \item{add_pars}{
    Use to add or change parameters.
  }
  \item{pars}{
    List of all parameters combined together.
    Expected to only be used internally with recursive calls.
  }
  \item{startupt}{
    Number of startup runs.
    Integer.
  }
  \item{starting}{
    Starting conditions.
  }
  \item{warn}{
    Set to \code{FALSE} to suppress warnings.
  }
  \item{quiet}{
    Set to \code{TRUE} to suppress some messages.
  }
}

\details{
  The function simulates anaerobic biodegradation.
}

\value{
  A data frame with output and some input variables.
  Each row is separated by the value set with \code{delta_t}.
  Typically the primary variable of interest is methane emission, which is returned as a total (g) and rate (g/d), overall or normalized to COD or VS loading.
  These columns all start with \code{CH4}.
  Microbial biomass (g) is given in columns with names that match those used for the names of the groups.
  For more details, see the vignette.
}

\references{
Dalby, F.R., Hafner, S.D., Petersen, S.O., Vanderzaag, A., Habtewold, J., Dunfield, K., Chantigny, M.H., Sommer, S.G., 2021. A mechanistic model of methane emission from animal slurry with a focus on microbial groups. PLOS ONE 16, e0252881. \doi{10.1371/journal.pone.0252881}

Dalby, F.R., Hansen, M.J., Guldberg, L.B., Hafner, S.D., Feilberg, A., 2023. Simple Management Changes Drastically Reduce Pig House Methane Emission in Combined Experimental and Modeling Study. Environ. Sci. Technol. \doi{10.1021/acs.est.2c08891}

Dalby, F.R., Ambrose, H.W.,  Poulsen, J.S., Nielsen, J.L., Adamsen, A.P.S., 2023. Pig slurry organic matter transformation and methanogenesis at ambient storage temperatures. JEQ. \doi{10.1002/jeq2.20512}

Dalby, F.R., Hafner, S.D., Ambrose, H.W., Adamsen, A.P.S., 2024. Pig manure degradation and carbon emission: Measuring and modeling combined aerobic–anaerobic transformations. JEQ. \doi{10.1002/jeq2.20603}

Sommer, S.G, Hafner S.D., Laubach, J., van der Weerden, T.J., Pacholski, A., 2022. Model for calculating ammonia emission from stored animal liquid manure. 
Biosyst Eng \doi{10.1016/j.biosystemseng.2022.08.007}
} 

\author{
  Sasha D. Hafner and Frederik R. Dalby
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{

out_dig <- abmneo(
  storage = storage_ex_dig_01,
  365,
  inf_pars = inf_pars_dig_01,
  grp_pars = grp_pars_dig_01,
  sub_pars = sub_pars_dig_01
)

plot(CH4 ~ time, data = out_dig, type = 'l')

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
