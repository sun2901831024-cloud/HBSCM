library(HBSCM)

smoking <- read.csv(system.file("extdata", "smoking.csv", package = "HBSCM"))
smoking$D <- as.integer(smoking$state == "California" & smoking$year >= 1989)
stopifnot(any(smoking$D == 1))

# Short example of the no-covariate specification in a4.R.
set.seed(1234567)
fit <- HBSCM(smoking, c("state", "year"), "cigsale", "D", re = "both",
             r = 10, niter = 120, burn = 60)
effects <- heffSummary(fit, cumu = TRUE)
print(effects$est.avg)
print(head(heffSummary(fit, rela.period = FALSE)$est.eff))
print(effPlot(effects, "eff", xlab = "Periods relative to adoption",
              ylab = "Cigarette sales per capita"))
