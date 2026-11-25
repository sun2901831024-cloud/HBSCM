\name{HBSCM}
\alias{HBSCM}
\alias{bpCausalgroup}
\alias{bpCausaladp}
\title{Fit a Bayesian synthetic control model}
\description{Estimate posterior predictive untreated outcomes in a panel with
unit-specific effects, time-varying effects and latent factors.}
\usage{
HBSCM(data, index, Yname, Dname, Xname = NULL, Zname = NULL,
  Aname = NULL, re = "unit", ar1 = TRUE, r = 2L,
  niter = 10000L, burn = 5000L, xhorseshoe = TRUE,
  zhorseshoe = TRUE, ahorseshoe = TRUE, fhorseshoe = TRUE)
bpCausalgroup(data, index, Yname, Dname, Xname = NULL, Zname = NULL,
  Aname = NULL, re = "unit", ar1 = TRUE, r = 2L,
  niter = 10000L, burn = 5000L, xlasso = TRUE, zlasso = TRUE,
  alasso = TRUE, flasso = TRUE, groupindex = NULL, group = "group")
bpCausaladp(data, index, Yname, Dname, Xname = NULL, Zname = NULL,
  Aname = NULL, re = "unit", ar1 = TRUE, r = 2L,
  niter = 10000L, burn = 5000L, xlasso = TRUE, zlasso = TRUE,
  alasso = TRUE, flasso = TRUE, a1 = 0.001, a2 = 0.001,
  b1 = 0.001, b2 = 0.001, c1 = 0.001, c2 = 0.001,
  p1 = 0.001, p2 = 0.001)
}
\arguments{
\item{data}{Long-format panel with one row per unit and period.}
\item{index}{Names of the unit and time columns, in that order.}
\item{Yname, Dname}{Names of the numeric outcome and binary treatment columns.}
\item{Xname, Zname, Aname}{Names of numeric covariates with common, unit-varying
and time-varying coefficients, respectively. NULL omits a covariate block.}
\item{re}{One of \code{"none"}, \code{"unit"}, \code{"time"}, \code{"both"}.
Adds the corresponding varying intercepts. A common intercept is always included.}
\item{ar1}{Estimate AR(1) coefficients for time-varying terms and factors.
If FALSE, the transition coefficients are fixed at zero.}
\item{r}{Nonnegative integer number of latent factors.}
\item{niter, burn}{Total iterations and initial iterations to discard.}
\item{xhorseshoe, zhorseshoe, ahorseshoe, fhorseshoe}{Enable horseshoe
shrinkage for common coefficients, unit scales, time scales and factor scales, respectively.}
\item{xlasso, zlasso, alasso, flasso}{Corresponding switches for the extension
scripts. The grouped model uses horseshoe updates; the adaptive model uses lasso updates.}
\item{group}{Name of a group column, constant within unit.}
\item{groupindex}{Optional legacy vector of integer group codes, one per retained
unit in sorted unit order. By default codes are derived from the group column.}
\item{a1, a2, b1, b2, c1, c2, p1, p2}{Positive shape and rate hyperparameters
for adaptive lasso mixing parameters in the four coefficient blocks.}
}
\details{Rows with missing model variables and units without untreated observations
are removed. At least one treated unit with an untreated baseline and at least
one untreated observation per retained time period are needed. Covariates are
divided by their sample standard deviation without centering. Fitting uses
untreated observations. Predictions include observation noise and span the
complete observed histories of the retained treated units.

For fidelity to the replication scripts, the horseshoe factor-global update
uses the number of time-varying columns in the shape and refreshes the
time-block auxiliary instead of the factor auxiliary. The grouped version
uses the supplied C++ conditional-mean update for group coefficients and loadings.}
\value{A list. Columns of \code{yct} are retained predictive draws and rows
follow treated-unit histories in unit and time order. \code{yo_t} gives the
corresponding observations. \code{beta} and \code{sigma2} store common
coefficient and residual variance draws. If present, \code{alpha}, \code{xi},
\code{gamma} and \code{f} are unit (or group) coefficients, time coefficients,
loadings and factors with draws on the last axis. \code{wa}, \code{wxi} and
\code{wg} contain their scales. Active shrinkage blocks return \code{sb2},
\code{swa2}, \code{swxi2} and \code{swg2}. AR draws are in \code{p},
\code{pxi} and \code{pf}. \code{unit.labels}, \code{time.labels} and
\code{covariate.sd} describe the encoding and standardization.}
\seealso{\code{\link{heffSummary}}, \code{\link{simulateSA}}}
\examples{
d <- simulateSA(N = 14, TT = 10, r = 2, tr.start = 7, p = 2,
                beta = c(1, 0), seed = 106)
set.seed(92)
fit <- HBSCM(d, c("id", "time"), "Y", "D", Xname = c("X1", "X2"),
             re = "both", r = 2, niter = 40, burn = 20)
heffSummary(fit)$est.avg
}
