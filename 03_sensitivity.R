# 03_sensitivity.R
# SENSITIVITY TO COMPENSATION, NATURAL MORTALITY, AND PACKET SIZE
# Run AFTER 01 and 02, in the same R session.
# Same parameter ranges as the manuscript, but fewer predation levels for
# a manageable public example. All populations are evaluated.

if (!exists("pops") || !exists("SMB_RECRUITS"))
  stop("Run scripts 01 and 02 first, in the same R session.")

SENS_FIXED_RATES <- c(0.4, 0.5, 0.6, 0.7, 0.8)
SENS_FMAX <- c(1.5, 2.0)
sens_rows <- list()

for (pk in names(pops)) {
  original <- pops[[pk]]

  # A. Vary compensation (CR), keeping calibrated R0 unchanged.
  for (cr in c(2, 3, 4)) {
    pop <- original
    pop$CR <- cr
    for (p in SENS_FIXED_RATES) {
      r <- simulate_fixed(pop, p, n_sim = N_SIM,
                          n_years = N_YEARS, seed = SEED)
      sens_rows[[length(sens_rows) + 1L]] <- data.frame(
        population = pk, test = "compensation_CR", value = cr,
        predation = p, p_qe_100 = tail(r$risk$p_qe, 1L))
    }
  }

  # B. Vary mortality scaling and recalibrate R0 to preserve adult N0.
  for (M in c(0.16, 0.25, 0.371, 0.598)) {
    pop <- original
    pop$Madult <- M
    lh <- build_life_history(pop)
    pop <- calibrate_population(pop, lh)
    for (p in SENS_FIXED_RATES) {
      r <- simulate_fixed(pop, p, n_sim = N_SIM,
                          n_years = N_YEARS, seed = SEED, lh = lh)
      sens_rows[[length(sens_rows) + 1L]] <- data.frame(
        population = pk, test = "natural_mortality_M", value = M,
        predation = p, p_qe_100 = tail(r$risk$p_qe, 1L))
    }
  }

  # C. Vary predation packet size only; natural survival stays individual.
  for (m in c(1L, 10L, 50L)) {
    for (Fmax in SENS_FMAX) {
      r <- simulate_dynamic(original, F_max = Fmax,
                            observed_recruits = SMB_RECRUITS,
                            packet_size = m, n_sim = N_SIM,
                            n_years = N_YEARS, seed = SEED)
      sens_rows[[length(sens_rows) + 1L]] <- data.frame(
        population = pk, test = "packet_size", value = m,
        predation = Fmax, p_qe_100 = tail(r$risk$p_qe, 1L))
    }
  }
  cat("03 sensitivity complete for ", pk, "\n", sep = "")
}

sensitivity <- do.call(rbind, sens_rows)
write.csv(sensitivity, "results/03_sensitivity.csv", row.names = FALSE)
cat("03 complete: results/03_sensitivity.csv\n")
