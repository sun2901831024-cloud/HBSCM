\name{effPlot}
\alias{effPlot}
\title{Plot posterior predictive outcomes or treatment effects}
\description{Build a plot from a treatment-effect summary.}
\usage{
effPlot(x, type, plotlength = 5, main = NULL, legend = TRUE, xlim = NULL,
  ylim = NULL, xlab = NULL, ylab = NULL, cex.xlab = 13, cex.ylab = 13,
  cex.xaxis = 10, cex.yaxis = 10, cex.main = 15, cex.legend = 10,
  x1.pos = 0, x2.pos = NULL, y.pos = 0, legend.pos = NULL,
  legend.labs = NULL, CI = TRUE)
}
\arguments{
\item{x}{Output of \code{\link{heffSummary}}.}
\item{type}{\code{"outcome"}, \code{"eff"}, or \code{"cumu"}. The last
requires a summary constructed with \code{cumu = TRUE}.}
\item{plotlength}{Reserved for compatibility with the original scripts.}
\item{main, xlab, ylab}{Plot title and axis labels.}
\item{legend, legend.pos, legend.labs}{Show the outcome legend, set its
position or supply labels for observed, predicted and optionally the interval.}
\item{xlim, ylim}{Optional display ranges.}
\item{cex.xlab, cex.ylab, cex.xaxis, cex.yaxis, cex.main, cex.legend}{Text sizes in points.}
\item{x1.pos, x2.pos, y.pos}{Vertical and horizontal reference positions. NULL omits a line.}
\item{CI}{Show the equal-tailed 95 percent interval.}
}
\value{A ggplot object.}
