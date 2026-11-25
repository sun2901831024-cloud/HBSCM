\name{heffSummary}
\alias{heffSummary}
\alias{coefSummary}
\title{Summarize treatment effects and coefficients}
\description{Summarize draws by posterior mean and equal-tailed 95 percent interval.}
\usage{
heffSummary(x, usr.id = NULL, burn = 0, cumu = FALSE, rela.period = TRUE)
coefSummary(x, burn = 0)
}
\arguments{
\item{x}{A fit from \code{\link{HBSCM}} or either extension.}
\item{usr.id}{Optional vector of original treated unit identifiers.}
\item{burn}{Additional draws to discard from the already trimmed fit.}
\item{cumu}{Include cumulative effects indexed by relative period.}
\item{rela.period}{Use treatment-relative periods. The last untreated period is
zero and the first treated period is one. FALSE uses common panel time indices and
adds original labels in \code{period_label}.}
}
\value{\code{heffSummary} returns \code{est.eff} with observed outcomes,
counterfactuals, effects, interval endpoints and counts by period; \code{eff_i}
contains the effect draws by period. \code{est.avg} summarizes the mean over
all treated observations, giving each unit-period equal weight. Optional
\code{est.cumu} gives the number of included post-treatment periods times
the observation-count-weighted mean effect through that period.

\code{coefSummary} returns \code{est.beta} and, when present, \code{est.alpha},
\code{est.xi}, \code{est.phi.xi}, and \code{est.phi.f}. Unit and time
coefficients incorporate their corresponding scale draws. Coefficients use the
standardized covariate scale.}
