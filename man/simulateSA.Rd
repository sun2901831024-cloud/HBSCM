\name{simulateSA}
\alias{simulateSA}
\title{Simulate a panel with staggered treatment adoption}
\description{Generate the panel design used in the H-BSCM simulations.}
\usage{
simulateSA(N, TT, r, tr.threshold = 0.5, tr.start = 10, p = 10,
  beta = c(4, 4, 2, 2, rep(0, 6)), alpha = NULL, xi = NULL,
  mu = 0, force = 3, factor = NULL, Rtype = "n", Ftype = "drift",
  Fsize = 1, eff.size = 1, eff.noise = 1, tr.noise = 0.3,
  time.invariant = TRUE, ar.coef = 0.7, error.type = "n", seed = NULL)
}
\arguments{
\item{N, TT}{Numbers of units and periods.}
\item{r}{Number of factors; the treatment rule uses the first two loadings.}
\item{tr.threshold, tr.start}{Adoption quantiles and matching start periods. For
multiple cohorts, later entries overwrite treatment onset for nested upper-tail units
as in the replication generator; supply increasing thresholds with earlier starts.}
\item{p, beta}{Number of covariates and corresponding common coefficients.}
\item{alpha, xi}{Optional coefficient matrices of size N by p and TT by p.}
\item{mu}{Common intercept.}
\item{force}{Add no varying intercepts (0), unit (1), time (2) or both (3).}
\item{factor}{Optional TT by r factor matrix.}
\item{Rtype}{Use uniform factor loadings with "u" or "uniform"; default normal.}
\item{Ftype, Fsize}{Factor processes ("ar1", "drift" or white noise) and multipliers.}
\item{eff.size, eff.noise, tr.noise}{Effect increment, effect noise SD and adoption noise SD.}
\item{time.invariant}{Hold covariates constant within each unit.}
\item{ar.coef}{AR coefficient used by the relevant simulated time series.}
\item{error.type}{Outcome noise: normal ("n"), standardized gamma ("g") or Student t with three degrees of freedom ("t").}
\item{seed}{Optional seed passed to \code{set.seed}.}
}
\value{A long data frame including unit \code{id}, period \code{time}, outcome
\code{Y}, treatment \code{D}, true effect \code{eff}, error, adoption information
and generated covariates, effects and factors.}
