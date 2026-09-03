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
    Influent and initial slurry parameters, as a named list with elements:
    \describe{
      \item{\code{VFA_fresh}, \code{VFA_init}}{
        Numeric scalars giving influent and initial concentration of the intermediate VFA pool, in the same concentration units as \code{xa_fresh}/\code{xa_init} in \code{grp_pars}.
      }
      \item{\code{pH}}{
        Numeric scalar giving slurry pH. Can also be supplied as a time-varying forcing variable via \code{var_pars}, in which case the \code{var_pars} value overrides this default at each time step.
      }
      \item{\code{dens}}{
        Numeric scalar giving slurry density (kg/m\eqn{^3}{^3}), used with \code{storage} depth and area to determine slurry mass.
      }
    }
    See \code{\link{inf_pars_dig_01}} for a complete example.
  }
  \item{grp_pars}{
    Parameters for microbial groups, as a named list with elements:
    \describe{
      \item{\code{grps}}{
        Character vector of microbial group names, e.g. \code{c('m1', 'm2', 'm3')} for methanogen groups 1, 2, and 3.
        Defines the number, names, and order of groups used by all other elements below.
      }
      \item{\code{yield}, \code{xa_fresh}, \code{xa_init}, \code{d_max}}{
        Named numeric vectors, one value per group in \code{grps}.
        \code{yield} is microbial biomass yield in g biomass COD / g consumed substrate COD.
        \code{xa_init} and \code{xa_fresh} are initial and influent microbial biomass concentration in g biomass COD / kg.
        And \code{d_max} is the maximum death/decay rate in 1/d.
        The name \code{"default"} can be used instead of group names to apply a single value to some or all groups, e.g. \code{c(default = 0.05)} or \code{c(m1 = 0.04, default = 0.05)} .
      }
      \item{\code{qhat_opt}, \code{T_opt}, \code{T_min}, \code{T_max}}{
        Named numeric vectors giving maximum substrate utilization rate (\code{qhat_opt}) and Cardinal Temperature Model parameters (maximum specific substrate uptake rate, and optimum/minimum/maximum temperature, respectively) for each group.
        As with the elements above, the name \code{"default"} can be used instead of group names to apply a single value to some or all groups.
      }
      \item{\code{ksmat}}{
        Numeric matrix of half-saturation constants, with one row per limiting reactant and one column per group.
        Row and column names must match reactant and group (\code{grps}) names. For example:
        \preformatted{
ksmat <- matrix(
  c(1.1, 1.1, 1.1, 1.1, 1.1),
  nrow = 1,
  byrow = TRUE,
  dimnames = list(
    c('VFA'),
    c('m1', 'm2', 'm3', 'm4', 'm5')
  )
)}
      }
      \item{\code{mstoich}}{
        Numeric matrix of microbial reaction stoichiometry, with one row per reactant/product and one column per group.
        Negative values indicate consumption, positive values production.
        Row and column names must match reactant/product and group (\code{grps}) names.
        Unlike the elements above, \code{mstoich} is \emph{not} expanded automatically to match \code{grps}, so it must already include a column for every group. For example:
        \preformatted{
mstoich <- matrix(
  c(
    -1, -1, -1, -1, -1,
     1,  1,  1,  1,  1
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    c('VFA', 'CH4'),
    c('m1', 'm2', 'm3', 'm4', 'm5')
  )
)}
      }
    }
    See \code{\link{grp_pars_dig_01}} for a complete example.
  }
  \item{sub_pars}{
    Substrate parameters, as a named list with elements:
    \describe{
      \item{\code{subs}}{
        Character vector of particulate substrate names, e.g. \code{c('PS')} for a single group particulate substrate.
        Defines the number, names, and order of substrates used by the elements below.
        \code{xd} (see below) is added to \code{subs} automatically if not already included.
      }
      \item{\code{xd}}{
        Character scalar giving the name of the substrate pool that receives dead/decayed microbial biomass, e.g. \code{'xd'}.
        Should match an element of \code{subs}; otherwise it is added automatically.
      }
      \item{\code{sub_enrich}}{
        Character vector, a subset of \code{subs}, naming substrates that are enriched (rather than diluted in proportion to residual slurry mass) during a removal event. See \code{resid_enrich} in the \code{storage} object.
      }
      \item{\code{h_rate_ref}, \code{h_rate_q10}}{
        Named numeric vectors, one value per substrate in \code{subs} (including \code{xd}).
        \code{h_rate_ref} is the hydrolysis rate constant at a reference temperature of 20 C (1/d), and \code{h_rate_q10} is the associated Q10 temperature coefficient.
        As with \code{grp_pars} vector elements, the name \code{"default"} can be used instead of substrate names to apply a single value to some or all substrates.
      }
      \item{\code{sub_fresh}, \code{sub_init}}{
        Named numeric vectors of influent and initial substrate concentration, one value per substrate in \code{subs} (including \code{xd}).
      }
    }
    See \code{\link{sub_pars_dig_01}} for a complete example.
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
