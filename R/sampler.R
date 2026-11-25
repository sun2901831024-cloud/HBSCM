# Initial scales and auxiliary variables in the replication implementation.
initial_scale <- function(k) {
  list(variance = rep(100, k), local_aux = rep(1, k),
       global = 1, global_aux = 1, lambda = rep(10, k))
}

update_scale <- function(value, state, hyper, shape_dimension, factor_block = FALSE) {
  k <- length(value)
  if (!is.null(hyper)) {
    for (j in seq_len(k)) {
      mu <- sqrt(state$lambda[j] / value[j]^2)
      state$variance[j] <- 1 / rrinvgauss(mu, state$lambda[j])
      state$lambda[j] <- sampleG(hyper[1] + 1, 1 / (hyper[2] + state$variance[j] / 2))
    }
    return(state)
  }
  local <- 1 / rgamma(k, 1, 0.5 * as.numeric(value)^2 / state$global + 1 / state$local_aux)
  state$local_aux <- 1 / rgamma(k, 1, 1 / local + 1)
  state$global <- 1 / rgamma(1, (1 + shape_dimension) / 2,
                             0.5 * sum(as.numeric(value)^2 / local) + 1 / state$global_aux)
  state$variance <- local * state$global
  if (!factor_block) state$global_aux <- 1 / rgamma(1, 1, 1 + 1 / state$global)
  state
}

sample_panel <- function(design, pred, y, id, time, id_tr, time_tr,
                         members, grouped, r, ar1, active, adaptive, niter, burn) {
  N <- length(members)
  TT <- max(time)
  G <- max(members)
  X <- design$X
  Z <- design$Z
  A <- design$A
  k1 <- ncol(X)
  k2 <- ncol(Z)
  k3 <- ncol(A)
  sizes <- c(k1, k2, k3, r)
  state <- lapply(sizes, initial_scale)
  flags <- c("sb2", "swa2", "swxi2", "swg2")
  keep <- niter - burn
  beta <- matrix(rnorm(k1), k1, 1L)
  alpha <- matrix(rnorm(G * k2), G, k2)
  xi <- matrix(rnorm(TT * k3), TT, k3)
  gamma <- matrix(rnorm(G * r), G, r)
  factors <- matrix(rnorm(TT * r), TT, r)
  wa <- rep(1, k2)
  wxi <- rep(1, k3)
  wg <- rep(1, r)
  phi <- matrix(0, k3 + r, 1L)
  tau_old <- cbind(xi, factors)
  by_time <- order(time, id)
  id_break <- breakID(as.matrix(id))
  time_break <- breakID(as.matrix(time[by_time]))
  expand_unit <- function(x) x[members, , drop = FALSE]
  scaled <- function(x, w) x * matrix(rep(w, each = nrow(x)), nrow(x), length(w))
  re_fit <- function(x, coef, w, index) {
    if (!ncol(x)) return(matrix(0, length(index), 1L))
    getREfit(x, scaled(coef, w), as.matrix(index))
  }
  factor_fit <- function(index, when) {
    if (r == 0L) return(matrix(0, length(index), 1L))
    getFactorFit(scaled(expand_unit(gamma), wg), factors, as.matrix(index), as.matrix(when))
  }
  xfit <- X %*% beta
  zfit <- re_fit(Z, expand_unit(alpha), wa, id)
  afit <- re_fit(A, xi, wxi, time)
  ffit <- factor_fit(id, time)
  sigma2 <- sampleSigmaE2(y - xfit - zfit - afit - ffit, 0.001, 0.001)
  out <- list(yct = matrix(0, length(id_tr), keep), sigma2 = matrix(0, 1L, keep),
              beta = matrix(0, k1, keep))
  if (k2) {
    out$alpha <- array(0, c(G, k2, keep))
    out$wa <- matrix(0, k2, keep)
  }
  if (k3) {
    out$xi <- array(0, c(TT, k3, keep))
    out$wxi <- matrix(0, k3, keep)
  }
  if (r) {
    out$gamma <- array(0, c(G, r, keep))
    out$f <- array(0, c(TT, r, keep))
    out$wg <- matrix(0, r, keep)
  }
  if (ar1 && k3 + r > 0) out$p <- matrix(0, k3 + r, keep)
  for (block in which(active)) out[[flags[block]]] <- matrix(0, sizes[block], keep)

  for (i in seq_len(niter)) {
    beta <- sampleN(X, y - zfit - afit - ffit, diag(state[[1]]$variance, k1), sigma2)
    xfit <- X %*% beta
    if (active[1]) state[[1]] <- update_scale(beta, state[[1]], adaptive[, 1], k1)
    if (k2 + r > 0) {
      tz <- genTildeZ(Z, factors, as.matrix(c(wa, wg)), as.matrix(time))
      if (grouped) {
        block <- sampleAlphag(tz, y - xfit - afit, diag(k2 + r), id_break,
                              as.matrix(members), N, r, G, sigma2)
      } else {
        block <- sampleAlpha(tz, y - xfit - afit, diag(k2 + r), id_break, N, r, sigma2)
      }
      alpha <- block[, seq_len(k2), drop = FALSE]
      gamma <- block[, k2 + seq_len(r), drop = FALSE]
      zfit <- re_fit(Z, expand_unit(alpha), wa, id)
    }
    if (k3 + r > 0) {
      ta <- genTildeZ(A, expand_unit(gamma), as.matrix(c(wxi, wg)), as.matrix(id))
      tau <- sampleXi(ta[by_time, , drop = FALSE], (y - xfit - zfit)[by_time, , drop = FALSE],
                      tau_old, phi, time_break, TT, sigma2)
      tau_old <- tau
      xi <- tau[, seq_len(k3), drop = FALSE]
      factors <- tau[, k3 + seq_len(r), drop = FALSE]
    }
    if (k2 + k3 + r > 0) {
      design_omega <- genTildeTau(Z, A, expand_unit(alpha), xi, expand_unit(gamma),
                                  factors, as.matrix(id), as.matrix(time))
      variance <- unlist(lapply(state[-1], "[[", "variance"), use.names = FALSE)
      omega <- sampleN(design_omega, y - xfit, diag(variance, length(variance)), sigma2)
      if (k2) {
        wa <- omega[seq_len(k2), , drop = FALSE]
        if (active[2]) state[[2]] <- update_scale(wa, state[[2]], adaptive[, 2], k2)
        signed <- permute(wa, alpha)
        wa <- signed$omega
        alpha <- signed$Xi
      }
      if (k3) {
        wxi <- omega[k2 + seq_len(k3), , drop = FALSE]
        if (active[3]) state[[3]] <- update_scale(wxi, state[[3]], adaptive[, 3], k3)
        signed <- permute(wxi, xi)
        wxi <- signed$omega
        xi <- signed$Xi
      }
      if (r) {
        wg <- omega[k2 + k3 + seq_len(r), , drop = FALSE]
        if (active[4]) {
          # The supplied and publicly released code use k3 in this factor
          # update and refresh the time-block auxiliary. Preserve that chain.
          state[[4]] <- update_scale(wg, state[[4]], adaptive[, 4], k3, factor_block = TRUE)
          if (is.null(adaptive))
            state[[3]]$global_aux <- 1 / rgamma(1, 1, 1 + 1 / state[[3]]$global)
        }
        for (j in seq_len(r)) {
          sign_draw <- runif(1)
          if (sign_draw > 0.25 && sign_draw <= 0.5) {
            gamma[, j] <- -gamma[, j]
            factors[, j] <- -factors[, j]
          } else if (sign_draw > 0.5 && sign_draw <= 0.75) {
            gamma[, j] <- -gamma[, j]
            wg[j] <- -wg[j]
          } else if (sign_draw > 0.75) {
            factors[, j] <- -factors[, j]
            wg[j] <- -wg[j]
          }
        }
      }
    }
    if (ar1 && k3 + r > 0) phi <- samplePhi(tau, diag(100, k3 + r))
    zfit <- re_fit(Z, expand_unit(alpha), wa, id)
    afit <- re_fit(A, xi, wxi, time)
    ffit <- factor_fit(id, time)
    sigma2 <- sampleSigmaE2(y - xfit - zfit - afit - ffit, 0.001, 0.001)
    yct <- pred$X %*% beta + re_fit(pred$Z, expand_unit(alpha), wa, id_tr) +
      re_fit(pred$A, xi, wxi, time_tr) + factor_fit(id_tr, time_tr) +
      matrix(rnorm(length(id_tr), sd = sqrt(sigma2)), ncol = 1L)
    if (i <= burn) next
    draw <- i - burn
    out$yct[, draw] <- yct
    out$sigma2[1, draw] <- sigma2
    out$beta[, draw] <- beta
    if (k2) {
      out$alpha[, , draw] <- alpha
      out$wa[, draw] <- wa
    }
    if (k3) {
      out$xi[, , draw] <- xi
      out$wxi[, draw] <- wxi
    }
    if (r) {
      out$gamma[, , draw] <- gamma
      out$f[, , draw] <- factors
      out$wg[, draw] <- wg
    }
    if (ar1 && k3 + r > 0) out$p[, draw] <- phi
    for (block in which(active)) out[[flags[block]]][, draw] <- state[[block]]$variance
  }
  if (ar1 && k3) out$pxi <- out$p[seq_len(k3), , drop = FALSE]
  if (ar1 && r) out$pf <- out$p[k3 + seq_len(r), , drop = FALSE]
  out
}
