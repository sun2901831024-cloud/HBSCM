library(HBSCM)

data <- simulateSA(14, 10, 2, tr.start = 7, p = 2, beta = c(1, 0),
                   time.invariant = FALSE, seed = 106)
data$group <- data$id %% 3L + 1L
args <- list(data = data, index = c("id", "time"), Yname = "Y", Dname = "D",
             Xname = c("X1","X2"), Zname = c("X1","X2"), Aname = c("X1","X2"),
             r = 2, niter = 16, burn = 5, re = "both")
for (model in list(HBSCM, bpCausalgroup, bpCausaladp)) {
  set.seed(92)
  fit <- do.call(model, args)
  stopifnot(ncol(fit$yct) == 11, all(is.finite(fit$yct)), all(fit$sigma2 > 0))
  change <- heffSummary(fit, cumu = TRUE)
  post <- fit$rela.time.tr > 0
  gap <- as.vector(fit$yo_t)[post] - fit$yct[post, , drop = FALSE]
  stopifnot(isTRUE(all.equal(unname(change$est.avg[1, 1]), mean(gap))),
            all(change$est.eff$counterfactual_ci_l <= change$est.eff$counterfactual_ci_u),
            nrow(coefSummary(fit)$est.alpha[[1L]]) == dim(fit$alpha)[1L])
  for (kind in c("outcome", "eff", "cumu"))
    stopifnot(inherits(effPlot(change, kind), "ggplot"))
}

# Two post periods formerly left the second cumulative estimate missing.
mini <- list(sigma2 = matrix(1, 1, 2), yo_t = matrix(c(10, 12, 16)),
             yct = matrix(c(8, 9, 10, 9, 11, 12), 3), raw.id.tr = rep("a", 3),
             rela.time.tr = 0:2, time.tr = 1:3)
h <- heffSummary(mini, cumu = TRUE)
stopifnot(identical(unname(h$est.cumu$mean), c(2, 7)),
          identical(unname(h$est.eff$counterfactual_ci_l),
                    apply(mini$yct, 1, quantile, 0.025, names = FALSE)))
# Single retained draw, null covariates, no latent factors and unused ID levels.
data$id <- factor(data$id, levels = c(0, sort(unique(data$id))))
set.seed(9)
minimal <- HBSCM(data, c("id","time"), "Y","D", re = "none", r = 0, niter = 3, burn = 2)
stopifnot(ncol(minimal$beta) == 1, nrow(coefSummary(minimal)$est.beta) == 1,
          all(is.finite(heffSummary(minimal)$est.avg)))
stopifnot(inherits(try(HBSCM(rbind(data, data[1,]), c("id","time"), "Y", "D"),
                       silent = TRUE), "try-error"))
