\name{grp_pars_dig_01}
\docType{data}
\alias{inf_pars_dig_01}
\alias{sub_pars_dig_01}
\alias{storage_ex_dig_01}
\title{ Default parameters and example storage input for \code{\link{abmneo}}. }
\description{
  Default microbial group, influent, and substrate parameters, plus an example \code{storage} input, for the \code{\link{abmneo}} function.
  A single set "dig_01" is presently available, and was fit to methane emission measurements from a set of three digestate storage tanks in Denmark.
  Measurements are described in Kasper et al. (2026) and original measurements were shared by first author P. Kasper.
  Note that the four objects are meant to be used together; the `resid_enrich` value in \code{storage_ex_dig_01} was fit to measurements along with the hydrolysis constant \code{h_rate_ref}, while other parameters were loosely based on ABM parameter sets or "expert judgement".
}

\usage{
grp_pars_dig_01
inf_pars_dig_01
sub_pars_dig_01
storage_ex_dig_01
}
\format{ All parameter objects, and the storage object, are lists. See \code{\link{abmneo}} for required elements. }
\references{
Kasper, P.L., Lycus, P.M. & Holm, M., 2026. Field-scale assessment of ammonia and methane emissions from tent-covered slurry storages using high-resolution tracer decay measurements. Biosystems Engineering, 270, p. 104551. \doi{10.3168/jds.2012-5340}.
}

\examples{
out <- abmneo(
  storage = storage_ex_dig_01,
  365,
  inf_pars = inf_pars_dig_01,
  grp_pars = grp_pars_dig_01,
  sub_pars = sub_pars_dig_01,
  startup = 1
)
plot(CH4 ~ time, data = out, type = 'l')
}

\keyword{datasets}
