# Public signatures retain the argument names used by the replication scripts.
HBSCM <- function(data, index, Yname, Dname, Xname = NULL, Zname = NULL,
                  Aname = NULL, re = "unit", ar1 = TRUE, r = 2L,
                  niter = 10000L, burn = 5000L, xhorseshoe = TRUE,
                  zhorseshoe = TRUE, ahorseshoe = TRUE, fhorseshoe = TRUE) {
  panel_fit(data, index, Yname, Dname, Xname, Zname, Aname, re, ar1, r,
            niter, burn, list(xhorseshoe, zhorseshoe, ahorseshoe, fhorseshoe))
}

bpCausalgroup <- function(data, index, Yname, Dname, Xname = NULL, Zname = NULL,
                          Aname = NULL, re = "unit", ar1 = TRUE, r = 2L,
                          niter = 10000L, burn = 5000L, xlasso = TRUE,
                          zlasso = TRUE, alasso = TRUE, flasso = TRUE,
                          groupindex = NULL, group = "group") {
  panel_fit(data, index, Yname, Dname, Xname, Zname, Aname, re, ar1, r,
            niter, burn, list(xlasso, zlasso, alasso, flasso),
            group = group, groupindex = groupindex)
}

bpCausaladp <- function(data, index, Yname, Dname, Xname = NULL, Zname = NULL,
                        Aname = NULL, re = "unit", ar1 = TRUE, r = 2L,
                        niter = 10000L, burn = 5000L, xlasso = TRUE,
                        zlasso = TRUE, alasso = TRUE, flasso = TRUE,
                        a1 = 0.001, a2 = 0.001, b1 = 0.001, b2 = 0.001,
                        c1 = 0.001, c2 = 0.001, p1 = 0.001, p2 = 0.001) {
  hyper <- list(a1, a2, b1, b2, c1, c2, p1, p2)
  if (!all(vapply(hyper, function(x) is.numeric(x) && length(x) == 1L &&
                  is.finite(x) && x > 0, logical(1)))) {
    stop("Adaptive lasso hyperparameters must be positive finite numbers.")
  }
  panel_fit(data, index, Yname, Dname, Xname, Zname, Aname, re, ar1, r,
            niter, burn, list(xlasso, zlasso, alasso, flasso),
            adaptive = matrix(unlist(hyper), 2L, 4L))
}

is_count <- function(x, minimum = 0L) {
  is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= minimum && x <= .Machine$integer.max && x == floor(x)
}

is_switch <- function(x) {
  (is.logical(x) || is.numeric(x)) && length(x) == 1L && x %in% c(0, 1)
}

panel_fit <- function(data, index, Yname, Dname, Xname, Zname, Aname, re,
                      ar1, r, niter, burn, switches, group = NULL,
                      groupindex = NULL, adaptive = NULL) {
  if (!is_count(niter, 1L) || !is_count(burn) || burn >= niter)
    stop("Require integer iteration counts with 0 <= burn < niter.")
  if (!is_count(r)) stop("r must be a nonnegative integer.")
  if (!all(vapply(c(switches, list(ar1)), is_switch, logical(1))))
    stop("Shrinkage and ar1 flags must be TRUE, FALSE, 1 or 0.")
  re <- match.arg(re, c("none", "unit", "time", "both"))
  if (!is.data.frame(data)) stop("data must be a data frame.")
  roles <- c(index, Yname, Dname, group)
  if (length(index) != 2L || length(Yname) != 1L || length(Dname) != 1L ||
      length(group) > 1L || anyDuplicated(roles) ||
      !is.character(roles)) stop("Specify distinct unit, time, outcome, treatment and group columns.")
  covars <- unique(c(Xname, Zname, Aname))
  if (any(covars %in% roles)) stop("Covariates must be distinct from panel index, outcome, treatment and group.")
  columns <- unique(c(roles, covars))
  absent <- setdiff(columns, names(data))
  if (length(absent)) stop("Unknown columns: ", paste(absent, collapse = ", "))
  data <- as.data.frame(data)[, columns, drop = FALSE]
  data <- data[stats::complete.cases(data), , drop = FALSE]
  if (!nrow(data)) stop("No complete observations.")
  if (anyDuplicated(data[index])) stop("Duplicate unit-time observations.")
  data <- data[order(data[[index[1]]], data[[index[2]]]), , drop = FALSE]
  d <- data[[Dname]]
  if (!(is.numeric(d) || is.logical(d)) || !all(d %in% 0:1))
    stop("Treatment must be coded 0 or 1.")
  data <- data[data[[index[1]]] %in% data[[index[1]]][d == 0], , drop = FALSE]
  units <- unique(data[[index[1]]])
  periods <- sort(unique(data[[index[2]]]))
  id <- match(data[[index[1]]], units)
  time <- match(data[[index[2]]], periods)
  d <- data[[Dname]]
  untreated <- which(d == 0)
  treated_units <- sort(unique(id[d == 1]))
  if (!length(treated_units)) stop("At least one treated unit needs both untreated and treated observations.")
  if (length(unique(time[untreated])) != length(periods))
    stop("Each retained period needs an untreated observation.")
  for (nm in c(Yname, covars)) {
    if (!is.numeric(data[[nm]]) || any(!is.finite(data[[nm]])))
      stop("Expected finite numeric values in ", nm, ".")
  }
  scales <- vapply(data[covars], sd, numeric(1))
  if (any(!is.finite(scales) | scales <= 0)) stop("Covariates must have positive sample standard deviations.")
  for (nm in covars) data[[nm]] <- data[[nm]] / scales[[nm]]

  make_design <- function(names, intercept) {
    columns <- as.matrix(data[, names, drop = FALSE])
    if (intercept) columns <- cbind("(Intercept)" = rep(1, nrow(data)), columns)
    storage.mode(columns) <- "double"
    columns
  }
  design <- list(X = make_design(Xname, TRUE),
                 Z = make_design(Zname, re %in% c("both", "unit")),
                 A = make_design(Aname, re %in% c("both", "time")))
  active <- as.logical(unlist(switches)) & c(length(Xname) > 0L,
               length(Zname) > 0L, length(Aname) > 0L, r > 0L)
  members <- seq_along(units)
  group_labels <- NULL
  if (!is.null(group)) {
    labels <- data[[group]][match(members, id)]
    if (any(labels[id] != data[[group]])) stop("Group membership must be constant within each unit.")
    group_labels <- sort(unique(labels))
    members <- match(labels, group_labels)
    if (!is.null(groupindex) && !identical(as.numeric(groupindex), as.numeric(members)))
      stop("groupindex must match group codes in sorted retained-unit order.")
  }
  selected <- which(id %in% treated_units)
  reltime <- rep(NA_real_, nrow(data))
  for (unit in treated_units) {
    rows <- which(id == unit)
    observed <- integer(length(periods))
    status <- integer(length(periods))
    observed[time[rows]] <- 1L
    status[time[rows]] <- d[rows]
    reltime[rows] <- get_term(status, observed)[time[rows]]
  }
  select_design <- function(rows) lapply(design, function(x) x[rows, , drop = FALSE])
  Xtr <- select_design(selected)
  fit <- sample_panel(select_design(untreated), Xtr, as.matrix(data[[Yname]][untreated]),
                      id[untreated], time[untreated], id[selected], time[selected],
                      members, !is.null(group), r, ar1, active, adaptive, niter, burn)
  c(fit, list(raw.id.tr = as.matrix(data[[index[1]]][selected]), id.tr = as.matrix(id[selected]),
              time.tr = as.matrix(time[selected]), tr.unit.pos = treated_units,
              rela.time.tr = reltime[selected], yo_t = as.matrix(data[[Yname]][selected]),
              X.tr = Xtr$X, A.tr = Xtr$A, Z.tr = Xtr$Z,
              unit.labels = units, time.labels = periods, group.labels = group_labels,
              group.index = if (!is.null(group)) members, covariate.sd = scales))
}
