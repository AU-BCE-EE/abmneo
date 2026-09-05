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
    A named list that describes the storage structure.
    Units are kg for mass, m for length, and d for time.
    The first element \code{type} (character, \code{"regular"} or \code{"series"}) selects between two kinds of storage, which take different additional elements.
    Elements common to both types:
    \describe{
      \item{\code{storage_depth}}{
        Maximum (design) slurry depth in the storage.
        Combined with \code{area} and slurry density (\code{dens} in \code{inf_pars}) to get maximum slurry mass.
      }
      \item{\code{area}}{
        Storage area. Really simply horizontal area used to relate slurry mass and depth.
      }
      \item{\code{temp_C}}{
        Slurry temperature.
      }
      \item{\code{resid_enrich}}{
        Controls preferential retention of particulate matter (relative to the elements in \code{sub_enrich} in \code{sub_pars}) in the slurry remaining after a removal event, on a logit scale.
        A value of \code{0} means no enrichment, i.e., all components are removed in proportion to the fraction of slurry mass removed.
      }
    }
    For \code{type = "regular"}, storage is emptied at fixed intervals, and additional elements are required:
    \describe{
      \item{\code{slurry_mass}}{
        Initial slurry mass in storage.
      }
      \item{\code{slurry_prod_rate}}{
        Constant rate of slurry mass addition between emptying events.
      }
      \item{\code{empty_int}}{
        Interval between complete emptying events.
        Use \code{0} or \code{NA} for a storage that is never emptied.
      }
      \item{\code{resid_depth}}{
        Depth of slurry remaining immediately after an emptying event.
        Combined with \code{area} and \code{dens} to get residual slurry mass.
      }
    }
    For example:
    \preformatted{
stor_reg_ref <- list(
  type = 'regular',
  slurry_prod_rate = 5,
  storage_depth = 2,
  area = 0.65,
  temp_C = 20,
  resid_enrich = 0.6,
  slurry_mass = 6.5,
  resid_depth = 0.01,
  empty_int = 28
)}
    For \code{type = "series"}, slurry mass over time is instead given directly, and the additional required element is:
    \describe{
      \item{\code{dat}}{
        A data frame with columns \code{time} and \code{slurry_mass}, giving a (typically measured) slurry mass trajectory.
        A decrease in \code{slurry_mass} between rows is interpreted as a removal event; an increase reflects slurry addition.
        See \code{rem_align_method} in \code{ctrl_pars} for how the exact timing of a removal within an interval is determined.
      }
    }
    For example:
    \preformatted{
mass_series <- data.frame(
  time =        c( 0, 121, 122, 152,    153, 182, 214, 244, 365),
  slurry_mass = c(67,  80,   2,  13,   0.01,  11,  22,  32,  67)
)

stor_ser_ref <- list(
  type = 'series',
  dat = mass_series,
  storage_depth = 2,
  area = 0.7,
  temp_C = 18,
  resid_enrich = 0.6
)}
  }
  \item{days}{
    Total number of days in simulation.
  }
  \item{delta_t}{
    Time step in output (d).
  }
  \item{times}{
    Alternative to \code{delta_t}; use to return specific times only.
  }
  \item{inf_pars}{
    Influent and initial slurry parameters, as a named list with elements:
    \describe{
      \item{\code{VFA_fresh}, \code{VFA_init}}{
        Numeric values with influent and initial concentration of the intermediate VFA pool, in the same concentration units as \code{xa_fresh}/\code{xa_init} in \code{grp_pars}.
      }
      \item{\code{pH}}{
        Surry pH. Can also be supplied as a time-varying forcing variable via \code{var_pars}, in which case the \code{var_pars} value overrides this default at each time step.
      }
      \item{\code{dens}}{
        Slurry density (kg/m3), used with \code{storage} area and slurry mass to determine slurry depth.
      }
    }
    See \code{\link{inf_pars_dig_01}} for an example.
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
        The keyword name \code{default} can be used instead of group names to apply a single value to some or all groups, e.g. \code{c(default = 0.05)} or \code{c(m1 = 0.04, default = 0.05)} .
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
      \item{\code{ic0}, \code{ic100}}{
        Optional.
        A pair of numeric matrices, with matching dimensions, that add inhibition of microbial activity by one or more inhibitors (inhibiting components).
        Each matrix has one row per inhibitor and one column per microbial group (matching \code{grps}). \code{ic0} gives the inhibitor level at which inhibition begins (0\% effect), and \code{ic100} the level at which it is complete (100\% effect).
        Inhibition is linear between these values.
        A row name ending in \code{"_conc"} (e.g. \code{"VFA_conc"}) refers to that component's concentration in the simulated slurry and any other row name (e.g. \code{"pH"}) is matched directly to a scalar parameter of the same name (from \code{inf_pars}, \code{storage}, etc.), which can be constant or time-varying via \code{var_pars}.
        A row name can be repeated to give the same inhibitor two independent effects, e.g. inhibiting at both low and high values, by reversing the order of the \code{ic0}/\code{ic100} pair (a decreasing rather than increasing effect with the inhibitor level) for the second occurrence (the only obvious application is for pH).
        Effects of multiple inhibitors on a group are combined multiplicatively.
        If omitted (the default), no inhibition is applied. 
        For example, inhibition by VFA concentration, and by pH at both low and high values:
        \preformatted{
ic0 <- matrix(
  c(8, 8, 8,
    8, 8, 8,
    6, 6, 6),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(
    c('VFA_conc', 'pH', 'pH'),
    c('m1', 'm2', 'm3')
  )
)

ic100 <- matrix(
  c(15, 15, 15,
     9,  9,  9,
     5,  5,  5),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(
    c('VFA_conc', 'pH', 'pH'),
    c('m1', 'm2', 'm3')
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
    See \code{\link{sub_pars_dig_01}} for an example.
  }
  \item{chem_pars}{
    Chemical parameters, as a named list with elements:
    \describe{
      \item{\code{COD_conv}}{
        Named numeric vector giving the COD-equivalent weight (g COD / g compound) of each component.
        Used to convert between the COD-based mass balance used internally (via \code{mstoich} etc.) and the reported mass of that substance.
        A value of \code{0} is appropriate for a fully oxidized product (e.g. \code{CO2}), and a negative value for an oxidized reactant (e.g. \code{SO4}).
      }
      \item{\code{gases}}{
        Character vector naming which components are tracked as gases.
        These get their own cumulative output columns and are excluded from removal during a storage emptying event (see \code{storage}), unlike solutes and particulate components.
      }
    }
    See default values above.
    Only gases actually produced (i.e., present as a row in \code{mstoich}) are kept.
  }
  \item{ctrl_pars}{
    Control parameters, affecting model mechanics rather than the represented system, as a named list with elements:
    \describe{
      \item{\code{rem_align_method}}{
        Character, one of \code{"mid"} (default), \code{"late"}, or \code{"early"}.
        Relevant only for \code{type = "series"} \code{storage} input, where the exact time of a removal event within the interval spanned by two mass observations is not directly given.
        \code{"early"} places the removal at the earlier of the two times; \code{"late"} at the later time; \code{"mid"} instead inserts an additional row at the midpoint time, so slurry accumulates up to the midpoint, is emptied there, and accumulates again up to the next given time.
      }
      \item{\code{var_fill_method}}{
        Character, \code{"interp"} (default) or \code{"forward"}.
        Determines how values in \code{var_pars} (and in \code{storage} for \code{type = "series"}) are filled in for simulation times that fall between the times explicitly provided.
        \code{"interp"} linearly interpolates; \code{"forward"} carries the last observed value forward.
        \code{"interp"} is not possible for a time-varying input that requires more than one value per time step (e.g. one value per microbial group), in which case \code{"forward"} is used instead, with a warning.
      }
      \item{\code{par_key}}{
        Character, treated as a regular expression (default matches a literal period).
        Used to split compound names supplied through \code{add_pars}, e.g. \code{qhat_opt.m1}, into a parameter name and an element (group or substrate) name.
      }
      \item{\code{h_rate_max_temp}}{
        Numeric.
        Maximum temperature (K) used when calculating the Q10 temperature dependence of hydrolysis rate (\code{h_rate_ref}/\code{h_rate_q10} in \code{sub_pars}); temperatures above this are capped, to avoid unrealistically high rates at high temperature.
      }
    }
  }
  \item{var_pars}{
    Optional time-varying inputs, as a named list of one or more data frames.
    The main (and typically only) element, \code{var}, is a data frame with a \code{time} column plus one column per input to vary over time; column names should match another parameter name (see \code{inf_pars}, \code{storage}, etc.) to override that value at the given times, with values in between determined by \code{var_fill_method} in \code{ctrl_pars}.
    Most scalar inputs can, in principle, be varied this way, although not all combinations have been tested.
    Additional named elements besides \code{var} can be used for inputs that need more than one column per time (e.g., a value for each microbial group); these are combined into \code{var} automatically.
    These data frames are included with a list element name that matches the parameter name, e.g., `list(var = ..., qhat_opt = data.frame(time = 1:2, m1 = c(0.1, 0.2), m2 = c(0, 0.1))` 
    See the vignette for an example with temperature.
  }
  \item{add_pars}{
    Use to add or change parameters.
    A named list that uses parameter names.
    For parameters that are normally nested, e.q., microbial or substrate parameters, use a dot `.` to separate the main and element names, e.g. `list(qhat_opt.m1 = 0.01)` for the microbial group `m1`.
  }
  \item{pars}{
    List of all parameters combined together.
    Expected to only be used internally with recursive calls.
  }
  \item{startup}{
    Number of startup runs.
    Integer.
  }
  \item{starting}{
    Starting conditions supplied as output directly from a previous \code{abmneo} run.
    Only concentrations are extracted and used.
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
