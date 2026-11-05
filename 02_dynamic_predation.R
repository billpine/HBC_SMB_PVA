# 02_dynamic_predation.R
# DYNAMIC SMALLMOUTH BASS PREDATION
# Run AFTER 01_fixed_predation.R, in the same R session.
# Uses the model functions and recruitment values defined in 01.

if (!exists("pops") || !exists("simulate_dynamic"))
  stop("Run 01_fixed_predation.R first, in the same R session.")

F_LEVELS <- c(0, 0.5, 1.0, 1.5, 2.0)
PACKET_PRIMARY <- c(DNM = 1L, LCR = 50L, WGC = 50L)

dynamic_risk_rows <- list()
dynamic_settled_rows <- list()

for (pk in names(pops)) {
  pop <- pops[[pk]]
  for (Fmax in F_LEVELS) {
    cat(sprintf("02 dynamic: %s, Fmax=%.1f, packet=%d\n",
                pk, Fmax, PACKET_PRIMARY[[pk]]))
    r <- simulate_dynamic(pop, F_max = Fmax,
                          observed_recruits = SMB_RECRUITS,
                          packet_size = PACKET_PRIMARY[[pk]],
                          n_sim = N_SIM, n_years = N_YEARS, seed = SEED)
    for (yr in c(20L, 100L)) {
      k <- match(yr, r$years)
      dynamic_risk_rows[[length(dynamic_risk_rows) + 1L]] <- data.frame(
        population = pk, Fmax = Fmax, packet_size = r$packet_size,
        peak_mortality = r$peak_mortality,
        mean_realized_mortality = r$weighted_realized_mortality,
        year = yr, p_qe = r$risk$p_qe[k],
        p_zero_adults = r$risk$p_zero_adults[k],
        p_modeled_extinction = r$risk$p_modeled_extinction[k]
      )
    }
    avg <- settled_adults(r)
    dynamic_settled_rows[[length(dynamic_settled_rows) + 1L]] <- data.frame(
      population = pk, Fmax = Fmax,
      packet_size = r$packet_size,
      peak_mortality = r$peak_mortality,
      mean_realized_mortality = r$weighted_realized_mortality,
      settled_mean = mean(avg),
      settled_lo = unname(quantile(avg, 0.025)),
      settled_hi = unname(quantile(avg, 0.975)),
      p_qe_100 = tail(r$risk$p_qe, 1L)
    )
  }
}

dynamic_risk <- do.call(rbind, dynamic_risk_rows)
dynamic_settled <- do.call(rbind, dynamic_settled_rows)
write.csv(dynamic_risk, "results/02_dynamic_risk.csv", row.names = FALSE)
cat("02 complete: results/02_dynamic_risk.csv\n")
