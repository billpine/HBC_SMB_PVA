# 04_demographic_decomposition.R
# CONTROLLED COMPARISON OF DEMOGRAPHIC ASSUMPTIONS
# Run AFTER 01-03, in the same R session.
# Common reference: LCR demographics. Change one factor at a time.

if (!exists("pops") || !exists("controlled_survival_lh"))
  stop("Run 01_fixed_predation.R first, in the same R session.")

reference <- pops$LCR
DECOMP_RATES <- seq(0, 1, by = 0.1)
decomp_rows <- list()

# Collect year-20/year-100 QE risk for one controlled scenario.
add_case <- function(group, label, pop, lh) {
  pop <- calibrate_population(pop, lh)
  rows <- lapply(DECOMP_RATES, function(p) {
    r <- simulate_fixed(pop, p, n_sim = N_SIM,
                        n_years = N_YEARS, seed = SEED, lh = lh)
    data.frame(
      experiment = group, scenario = label,
      N0 = pop$N0, CR = pop$CR, R0 = pop$R0,
      fixed_mortality = p,
      p_qe_20 = r$risk$p_qe[match(20L, r$years)],
      p_qe_100 = tail(r$risk$p_qe, 1L),
      first_year_50pct_risk = first_year_at_risk(
        r$risk$p_qe, r$years, target = 0.5)
    )
  })
  do.call(rbind, rows)
}

# A. Starting abundance: life history and CR are held at LCR values.
for (N0 in c(500L, 5000L, 35000L)) {
  pop <- reference
  pop$N0 <- N0
  decomp_rows[[length(decomp_rows) + 1L]] <- add_case(
    "initial_abundance", paste0("N0=", N0), pop,
    build_life_history(reference))
}

# B. Compensation: initial abundance and life history held at LCR values.
for (cr in c(2, 3, 4)) {
  pop <- reference
  pop$CR <- cr
  decomp_rows[[length(decomp_rows) + 1L]] <- add_case(
    "compensation", paste0("CR=", cr), pop,
    build_life_history(reference))
}

# C. Survival schedule: fecundity stays at LCR values, survival changes.
for (i in seq_len(2L)) {
  K <- c(0.36, 0.34)[i]
  M <- c(0.598, 0.371)[i]
  lh <- controlled_survival_lh(reference, K_survival = K, M_survival = M)
  decomp_rows[[length(decomp_rows) + 1L]] <- add_case(
    "survival_schedule", sprintf("K=%.2f, M=%.3f", K, M),
    reference, lh)
}

decomposition <- do.call(rbind, decomp_rows)
write.csv(decomposition, "results/04_demographic_decomposition.csv",
          row.names = FALSE)
cat("04 complete: results/04_demographic_decomposition.csv\n")
