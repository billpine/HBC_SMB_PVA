# 05_settled_abundance.R
# LONG-TERM ADULT ABUNDANCE RELATIVE TO NO-PREDATION CASES
# Run LAST, after scripts 01-04 in the same R session.
# Uses results already computed by 01 and 02; no new simulations.

if (!exists("fixed_settled") || !exists("dynamic_settled"))
  stop("Run scripts 01 and 02 first, in the same R session.")

# Each settled_mean is the across-replicate mean of adult-female abundance
# averaged over years 51-100. Dividing by the corresponding no-predation
# scenario describes declines even if the QE threshold was never crossed.
# 'realized_mortality' is weighted by the number of juveniles exposed.

fixed_out <- fixed_settled
fixed_out$type <- "fixed"
fixed_out$scenario <- fixed_out$fixed_mortality
fixed_out$peak_mortality <- NA_real_
fixed_out$realized_mortality <- fixed_out$weighted_realized_mortality
fixed_out$packet_size <- NA_integer_

dynamic_out <- dynamic_settled
dynamic_out$type <- "dynamic"
dynamic_out$scenario <- dynamic_out$Fmax
dynamic_out$realized_mortality <- dynamic_out$mean_realized_mortality

all_settled <- rbind(
  fixed_out[, c("type", "population", "scenario", "packet_size",
                "peak_mortality", "realized_mortality", "settled_mean",
                "settled_lo", "settled_hi", "p_qe_100")],
  dynamic_out[, c("type", "population", "scenario", "packet_size",
                  "peak_mortality", "realized_mortality", "settled_mean",
                  "settled_lo", "settled_hi", "p_qe_100")]
)

all_settled$relative_to_no_predation <- NA_real_
all_settled$deterministic_expected <- NA_real_
for (mode in c("fixed", "dynamic")) {
  for (pk in names(pops)) {
    subset <- all_settled$type == mode & all_settled$population == pk
    baseline <- all_settled$settled_mean[subset & all_settled$scenario == 0]
    if (length(baseline) != 1L || baseline <= 0)
      stop("Missing no-predation reference for ", pk, "/", mode)
    all_settled$relative_to_no_predation[subset] <-
      all_settled$settled_mean[subset] / baseline
  }
}

# Beverton-Holt deterministic expectation for fixed predation only.
is_fixed <- all_settled$type == "fixed"
CR <- vapply(pops, function(p) p$CR, numeric(1))
all_settled$deterministic_expected[is_fixed] <- pmax(
  0, (CR[all_settled$population[is_fixed]] *
        (1 - all_settled$scenario[is_fixed]) - 1) /
       (CR[all_settled$population[is_fixed]] - 1)
)

write.csv(all_settled, "results/05_settled_abundance.csv",
          row.names = FALSE)
cat("05 complete: results/05_settled_abundance.csv\n")
