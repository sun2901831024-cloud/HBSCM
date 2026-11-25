effPlot <- function(x, type, plotlength = 5, main = NULL, legend = TRUE,
                    xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
                    cex.xlab = 13, cex.ylab = 13, cex.xaxis = 10, cex.yaxis = 10,
                    cex.main = 15, cex.legend = 10, x1.pos = 0, x2.pos = NULL,
                    y.pos = 0, legend.pos = NULL, legend.labs = NULL, CI = TRUE) {
  type <- match.arg(type, c("outcome", "eff", "cumu"))
  tab <- x$est.eff
  if (type == "cumu") {
    if (is.null(x$est.cumu)) stop("Use heffSummary(..., cumu = TRUE) first.")
    tab <- rbind(data.frame(mean = 0, ci_l = 0, ci_u = 0, count = 0, time = 0), x$est.cumu)
  }
  tab$value <- switch(type, outcome = tab$estimated_counterfactual, eff = tab$estimated_ATT, cumu = tab$mean)
  tab$lower <- switch(type, outcome = tab$counterfactual_ci_l, eff = tab$estimated_ATT_ci_l, cumu = tab$ci_l)
  tab$upper <- switch(type, outcome = tab$counterfactual_ci_u, eff = tab$estimated_ATT_ci_u, cumu = tab$ci_u)
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$time))
  if (CI) p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                                          fill = "grey60", alpha = 0.35)
  if (type == "outcome") {
    if (is.null(legend.labs)) legend.labs <- c("Observed", "Estimated Y(0)", "95% interval")
    p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$observed, linetype = "observed")) +
      ggplot2::geom_line(ggplot2::aes(y = .data$value, linetype = "counterfactual")) +
      ggplot2::scale_linetype_manual(name = NULL, breaks = c("observed", "counterfactual"),
                  values = c(observed = "solid", counterfactual = "dashed"), labels = legend.labs[1:2])
    if (CI && length(legend.labs) >= 3L) p <- p + ggplot2::labs(caption = legend.labs[3L])
  } else {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$value))
  }
  for (pos in c(x1.pos, x2.pos)) p <- p + ggplot2::geom_vline(xintercept = pos, linetype = "dashed", colour = "grey55")
  if (!is.null(y.pos)) p <- p + ggplot2::geom_hline(yintercept = y.pos, linetype = "dashed", colour = "grey55")
  if (is.null(legend.pos)) legend.pos <- if (legend) "bottom" else "none"
  if (is.null(main)) main <- switch(type, outcome = "Observed and counterfactual outcomes",
                                    eff = "Treatment effects", cumu = "Cumulative treatment effects")
  if (is.null(ylab)) ylab <- switch(type, outcome = "Outcome", eff = "ATT", cumu = "Cumulative ATT")
  if (is.null(xlab)) xlab <- "Time"
  p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
    ggplot2::labs(title = main, x = xlab, y = ylab) + ggplot2::theme_bw() +
    ggplot2::theme(legend.position = legend.pos, legend.text = ggplot2::element_text(size = cex.legend),
      axis.title.x = ggplot2::element_text(size = cex.xlab), axis.title.y = ggplot2::element_text(size = cex.ylab),
      axis.text.x = ggplot2::element_text(size = cex.xaxis), axis.text.y = ggplot2::element_text(size = cex.yaxis),
      plot.title = ggplot2::element_text(size = cex.main))
}
