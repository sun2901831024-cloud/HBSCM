retained <- function(x, burn) {
  n <- ncol(x$sigma2)
  if (is.null(n) || !is_count(burn) || burn >= n)
    stop("burn must leave at least one retained draw.")
  seq.int(burn + 1L, n)
}

interval_rows <- function(draws) {
  limits <- t(apply(draws, 1L, quantile, probs = c(0.025, 0.975), names = FALSE))
  data.frame(mean = rowMeans(draws), ci_l = limits[, 1L], ci_u = limits[, 2L])
}

heffSummary <- function(x, usr.id = NULL, burn = 0, cumu = FALSE, rela.period = TRUE) {
  draws <- retained(x, burn)
  if (!is_switch(cumu) || !is_switch(rela.period)) stop("Invalid summary flag.")
  if (cumu) rela.period <- TRUE
  ids <- as.vector(x$raw.id.tr)
  if (!is.null(usr.id) && (!length(usr.id) || !all(usr.id %in% ids)))
    stop("usr.id must identify treated units in this fit.")
  rows <- if (is.null(usr.id)) seq_along(ids) else which(ids %in% usr.id)
  relative <- x$rela.time.tr[rows]
  period <- if (rela.period) relative else as.vector(x$time.tr)[rows]
  if (all(is.na(period))) stop("No eligible periods for this summary.")
  observed <- as.vector(x$yo_t)[rows]
  yct <- x$yct[rows, draws, drop = FALSE]
  periods <- sort(unique(period[!is.na(period)]))
  groups <- split(seq_along(period), factor(period, levels = periods))
  count <- lengths(groups)
  obs <- vapply(groups, function(g) mean(observed[g]), numeric(1))
  pred <- do.call(rbind, lapply(groups, function(g) colMeans(yct[g, , drop = FALSE])))
  cf <- interval_rows(pred)
  eff <- as.numeric(obs) - pred
  effect <- interval_rows(eff)
  by_period <- data.frame(observed = as.numeric(obs), estimated_counterfactual = cf$mean,
                          counterfactual_ci_l = cf$ci_l, counterfactual_ci_u = cf$ci_u,
                          estimated_ATT = effect$mean, estimated_ATT_ci_l = effect$ci_l,
                          estimated_ATT_ci_u = effect$ci_u, time = periods, count = as.integer(count))
  if (!rela.period && !is.null(x$time.labels))
    by_period$period_label <- x$time.labels[periods]
  post <- which(relative > 0)
  if (!length(post)) stop("No post-treatment observations in this selection.")
  avg <- colMeans(observed[post] - yct[post, , drop = FALSE])
  result <- list(est.eff = by_period, est.avg = as.matrix(interval_rows(matrix(avg, nrow = 1L))),
                 eff_i = eff)
  if (cumu) {
    keep <- which(periods > 0)
    post_eff <- eff[keep, , drop = FALSE]
    weights <- count[keep]
    cumulative <- post_eff
    for (j in seq_along(keep)) {
      through <- seq_len(j)
      cumulative[j, ] <- j * colSums(weights[through] * post_eff[through, , drop = FALSE]) /
        sum(weights[through])
    }
    result$est.cumu <- cbind(interval_rows(cumulative), count = as.integer(weights), time = periods[keep])
  }
  result
}

coefSummary <- function(x, burn = 0) {
  draws <- retained(x, burn)
  result <- list(est.beta = interval_rows(x$beta[, draws, drop = FALSE]))
  scaled_summary <- function(values, weights) {
    n <- dim(values)[1L]
    lapply(seq_len(dim(values)[2L]), function(j) {
      mat <- matrix(values[, j, draws, drop = FALSE], nrow = n)
      interval_rows(sweep(mat, 2L, weights[j, draws], "*"))
    })
  }
  if (!is.null(x$alpha)) {
    grouped <- !is.null(x$group.index)
    ids <- if (grouped) x$group.labels else x$unit.labels
    if (is.null(ids)) ids <- seq_len(dim(x$alpha)[1L])
    treated <- if (grouped) unique(x$group.index[x$tr.unit.pos]) else x$tr.unit.pos
    result$est.alpha <- lapply(scaled_summary(x$alpha, x$wa), function(tab) {
      tab$id <- ids
      tab$tr <- factor(as.integer(seq_along(ids) %in% treated), levels = 0:1)
      tab
    })
  }
  if (!is.null(x$xi)) {
    result$est.xi <- lapply(scaled_summary(x$xi, x$wxi), function(tab) {
      tab$time <- seq_len(nrow(tab))
      tab
    })
  }
  if (!is.null(x$pxi)) result$est.phi.xi <- interval_rows(x$pxi[, draws, drop = FALSE])
  if (!is.null(x$pf)) result$est.phi.f <- interval_rows(x$pf[, draws, drop = FALSE])
  result
}
