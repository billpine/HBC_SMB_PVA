# 01_fixed_predation.R
# HUMPBACK CHUB PVA — MODEL DEFINITION AND FIXED PREDATION
# Run this script FIRST. It defines the model functions used by 02-05.
# All inputs are included below. Requires base R only.

# CHANGE THIS ONE NUMBER FOR MORE MONTE CARLO REPLICATES:
N_SIM   <- 100L    # quick example; use 1000L for manuscript-scale runs
N_YEARS <- 100L
SEED    <- 42L

# These are the fixed juvenile mortality rates examined by this example.
# Includes the paper's 10% grid plus 25% and 75% trajectory scenarios.
FIXED_RATES <- sort(unique(c(seq(0, 1, by = 0.10), 0.25, 0.75)))

# Middle Green River SMB recruitment values digitized from Bruckerhoff et al.
# (Figure 8c). Included here to eliminate the need for an input data file.
SMB_RECRUITS <- c(
  2751, 5981, 1555, 1077, 478, 957, 239,
  5981, 1316, 2153, 1555, 2392, 2990, 7416,
  239, 9689, 13756, 3828
)

dir.create("results", showWarnings = FALSE)

# =====================================================================
# MODEL FUNCTIONS (shared with scripts 02-05)
# =====================================================================

# Humpback Chub population viability analysis: model functions
# ----------------------------------------------------------------
# Model equations adapted from the authors' manuscript code. Base R only.
# Female-only age-structured model, with age-1 recruitment, ages 1-30,
# age-4+ adults, Beverton-Holt density dependence, and stochastic survival.

hbc_populations <- function() {
  populations <- list(
    DNM = list(name = "Dinosaur National Monument", short = "DNM",
               N0 = 500L, CR = 2, K = 0.36, Linf_mm = 318,
               Madult = 0.598),
    LCR = list(name = "Little Colorado River", short = "LCR",
               N0 = 5000L, CR = 4, K = 0.36, Linf_mm = 318,
               Madult = 0.598),
    WGC = list(name = "Western Grand Canyon", short = "WGC",
               N0 = 35000L, CR = 4, K = 0.34, Linf_mm = 347,
               Madult = 0.371)
  )
  lapply(populations, function(pop) {
    pop$Wmat <- (200 / pop$Linf_mm)^3
    calibrate_population(pop, build_life_history(pop))
  })
}

build_life_history <- function(pop, max_age = 30L,
                               maturity_age = 4L, fecundity_scale = 600) {
  stopifnot(max_age >= maturity_age, maturity_age >= 2L,
            pop$K > 0, pop$Linf_mm > 0, pop$Madult >= 0)
  ages <- seq_len(max_age)
  rel_length <- 1 - exp(-pop$K * ages)
  rel_weight <- rel_length^3
  Wmat <- (200 / pop$Linf_mm)^3

  # The 200-mm value is the fecundity reference size; maturity is fixed at age 4.
  fecundity <- pmax(0, rel_weight - Wmat) * fecundity_scale
  fecundity[ages < maturity_age] <- 0

  # The original model parameterizes Lorenzen-type instantaneous mortality
  # with a 0.66 multiplier and relative length as the size proxy.
  # "Madult" is retained as the source code's parameter name; it is NOT
  # the actual annual instantaneous mortality of every adult age class.
  annual_survival <- pmin(exp(-pop$Madult / rel_length * 0.66), 1)

  survivorship <- numeric(max_age)
  survivorship[1L] <- 1
  survivorship[2L:max_age] <- cumprod(annual_survival[1L:(max_age - 1L)])

  list(ages = ages, rel_length = rel_length, rel_weight = rel_weight,
       fecundity = fecundity, annual_survival = annual_survival,
       survivorship = survivorship, maturity_age = maturity_age)
}

calibrate_population <- function(pop, lh) {
  # N0 denotes adult females age 4+, not all age-1+ females.
  adult_per_recruit <- sum(lh$survivorship[lh$ages >= lh$maturity_age])
  if (!is.finite(adult_per_recruit) || adult_per_recruit <= 0) {
    stop("Cannot calibrate equilibrium age-1 recruitment")
  }
  pop$R0 <- max(1L, round(pop$N0 / adult_per_recruit))
  pop$ext_thresh <- max(1L, round(0.05 * pop$N0))
  pop
}

beverton_holt_parameters <- function(pop, lh) {
  # Eggs per age-1 recruit from survivorship and fecundity at age.
  epr <- max(sum(lh$survivorship * lh$fecundity), 1e-12)
  list(alpha = pop$CR / epr,
       beta = (pop$CR - 1) / (pop$R0 * epr),
       eggs_per_recruit = epr)
}

# Return the fraction of simulations that have crossed each endpoint at
# least once at/before each year (cumulative event probability).
cumulative_risk <- function(N_adults, N_age1plus, threshold) {
  if (!identical(dim(N_adults), dim(N_age1plus))) {
    stop("Adult and age-1+ matrices must have matching dimensions")
  }
  running_minimum <- function(x) {
    # vapply preserves a years x replicates matrix with 1 replicate too.
    vapply(seq_len(ncol(x)), function(s) cummin(x[, s]),
           numeric(nrow(x)))
  }
  adults_min <- running_minimum(N_adults)
  all_ages_min <- running_minimum(N_age1plus)
  list(p_qe = rowMeans(adults_min <= threshold),
       p_zero_adults = rowMeans(adults_min <= 0),
       p_modeled_extinction = rowMeans(all_ages_min <= 0))
}

# Predation on juveniles in correlated groups of m. Each *complete* packet
# survives or dies as a group; the remaining (<m) fish survive independently.
# This reproduces the original dynamic-model algorithm, including its two
# binomial calls (even when the remainder is zero).
apply_packet_predation <- function(J_pre, pred_prob, m = 50L) {
  if (length(J_pre) != length(pred_prob)) {
    stop("Juvenile counts and predation probabilities have different lengths")
  }
  if (length(m) != 1L || !is.finite(m) || m < 1L ||
      m != as.integer(m)) stop("Packet size must be a positive integer")
  if (any(!is.finite(J_pre)) || any(J_pre < 0) ||
      any(!is.finite(pred_prob)) || any(pred_prob < 0 | pred_prob > 1)) {
    stop("Invalid juvenile counts or predation probabilities")
  }
  m <- as.integer(m)
  J_pre <- as.integer(round(J_pre))
  n_packets <- J_pre %/% m
  remainder <- J_pre %% m
  packet_survivors <- rbinom(length(J_pre), n_packets, 1 - pred_prob) * m
  remainder_survivors <- rbinom(length(J_pre), remainder, 1 - pred_prob)
  as.integer(packet_survivors + remainder_survivors)
}

# The SMB index in the submitted model is an illustrative, replicate-specific
# bootstrap of the Middle Green recruitment series defined below.
# Abundance recursion: A[t+1] = 0.5*A[t] + R[t], initialized at mean(R).
# Each replicate's 100-year abundance is normalized by ITS OWN maximum.
build_smb_index <- function(n_years, observed_recruits, seed = 42L) {
  if (n_years < 2L || length(observed_recruits) < 2L ||
      any(!is.finite(observed_recruits)) || any(observed_recruits < 0)) {
    stop("Invalid SMB recruitment series or projection length")
  }
  set.seed(seed)
  recruits <- sample(observed_recruits, size = n_years, replace = TRUE)
  smb_abundance <- numeric(n_years)
  smb_abundance[1L] <- mean(observed_recruits)
  for (t in seq_len(n_years - 1L)) {
    smb_abundance[t + 1L] <- 0.5 * smb_abundance[t] + recruits[t]
  }
  list(smb_abundance = smb_abundance,
       bootstrapped_recruits = recruits,
       smb_index = smb_abundance / max(smb_abundance))
}

build_smb_predation <- function(n_years, n_sim, observed_recruits,
                                F_max, seed = 42L) {
  stopifnot(n_years >= 2L, n_sim >= 1L, F_max >= 0)
  pred_mat <- matrix(NA_real_, nrow = n_years, ncol = n_sim)
  for (s in seq_len(n_sim)) {
    smb <- build_smb_index(n_years, observed_recruits, seed = seed + s)
    pred_mat[, s] <- 1 - exp(-F_max * smb$smb_index)
  }
  pred_mat
}

# Core projection engine. Fixed scenarios use independent binomial mortality;
# dynamic scenarios use packet predation and year x replicate probability.
# Like the source code, set.seed(seed) resets the HBC random stream for
# each scenario. Initializations/recruitment anomalies match across scenarios;
# later survival draws can diverge when fish counts differ.
simulate_hbc <- function(pop, lh = build_life_history(pop),
                         n_years = 100L, n_sim = 100L, seed = 42L,
                         fixed_predation = NULL, dynamic_predation = NULL,
                         packet_size = 50L, recruitment_sd = 0.50) {
  fixed <- !is.null(fixed_predation)
  dynamic <- !is.null(dynamic_predation)
  if (fixed == dynamic) stop("Supply exactly one predation mode")
  if (length(n_years) != 1L || n_years < 1L ||
      length(n_sim) != 1L || n_sim < 1L ||
      recruitment_sd < 0) stop("Invalid simulation settings")
  n_years <- as.integer(n_years)
  n_sim <- as.integer(n_sim)
  if (fixed) {
    if (length(fixed_predation) != 1L || !is.finite(fixed_predation) ||
        fixed_predation < 0 || fixed_predation > 1) {
      stop("Fixed predation must be a probability in [0, 1]")
    }
  } else {
    if (!is.matrix(dynamic_predation) ||
        !identical(dim(dynamic_predation), c(n_years, n_sim)) ||
        any(!is.finite(dynamic_predation)) ||
        any(dynamic_predation < 0 | dynamic_predation > 1)) {
      stop("Dynamic predation must be a valid n_years x n_sim matrix")
    }
  }

  max_age <- length(lh$ages)
  maturity_age <- lh$maturity_age
  set.seed(seed)
  bh <- beverton_holt_parameters(pop, lh)

  # Year-0 age structure: rmultinom on total age-1+ females, using stable
  # survivorship proportions. Calibrated expected adults ~= specified N0.
  init_props <- lh$survivorship / sum(lh$survivorship)
  adult_prop <- sum(init_props[lh$ages >= maturity_age])
  N0_age1plus <- as.integer(round(pop$N0 / adult_prop))
  N <- rmultinom(n_sim, size = N0_age1plus, prob = init_props)

  rec_anom <- matrix(
    exp(rnorm(n_years * n_sim, 0, recruitment_sd) -
          0.5 * recruitment_sd^2),
    nrow = n_years, ncol = n_sim
  )

  N_adults <- matrix(0, nrow = n_years + 1L, ncol = n_sim)
  N_age1plus <- matrix(0, nrow = n_years + 1L, ncol = n_sim)
  N_age_mean <- matrix(0, nrow = n_years + 1L, ncol = max_age)
  J_pre <- matrix(0, nrow = n_years, ncol = n_sim)
  J_post <- matrix(0, nrow = n_years, ncol = n_sim)

  N_adults[1L, ] <- colSums(N[maturity_age:max_age, , drop = FALSE])
  N_age1plus[1L, ] <- colSums(N)
  N_age_mean[1L, ] <- rowMeans(N)

  for (step in seq_len(n_years)) {
    eggs <- colSums(N * lh$fecundity)
    expected_juveniles <- (bh$alpha * eggs * rec_anom[step, ]) /
      (1 + bh$beta * eggs)
    J_pre[step, ] <- pmax(0L, round(expected_juveniles))

    if (fixed) {
      J_post[step, ] <- rbinom(n_sim, J_pre[step, ], 1 - fixed_predation)
    } else {
      J_post[step, ] <- apply_packet_predation(
        J_pre[step, ], dynamic_predation[step, ], m = packet_size)
    }

    N_next <- matrix(0L, nrow = max_age, ncol = n_sim)
    N_next[1L, ] <- J_post[step, ]
    for (a in seq_len(max_age - 1L)) {
      N_next[a + 1L, ] <- rbinom(n_sim, N[a, ], lh$annual_survival[a])
    }
    # Fish of the maximum modeled age do not enter an age-31+ group.
    N <- N_next
    state_row <- step + 1L
    N_adults[state_row, ] <- colSums(N[maturity_age:max_age, , drop = FALSE])
    N_age1plus[state_row, ] <- colSums(N)
    N_age_mean[state_row, ] <- rowMeans(N)
  }

  risk <- cumulative_risk(N_adults, N_age1plus, pop$ext_thresh)
  valid <- J_pre > 0
  mean_realized <- if (any(valid)) {
    mean((J_pre[valid] - J_post[valid]) / J_pre[valid])
  } else NA_real_
  weighted_realized <- if (sum(J_pre) > 0) {
    sum(J_pre - J_post) / sum(J_pre)
  } else NA_real_

  list(population = pop$short, pop = pop, years = 0:n_years,
       transition_years = seq_len(n_years), N_adults = N_adults,
       N_age1plus = N_age1plus, N_age_mean = N_age_mean,
       J_pre = J_pre, J_post = J_post, risk = risk,
       predation_type = if (fixed) "fixed" else "dynamic",
       fixed_predation = if (fixed) fixed_predation else NA_real_,
       packet_size = if (dynamic) as.integer(packet_size) else NA_integer_,
       mean_imposed_predation = if (fixed) fixed_predation else mean(dynamic_predation),
       mean_realized_mortality = mean_realized,
       weighted_realized_mortality = weighted_realized)
}

simulate_fixed <- function(pop, pred_rate, n_sim = 100L,
                           n_years = 100L, seed = 42L,
                           lh = build_life_history(pop)) {
  simulate_hbc(pop, lh = lh, n_sim = n_sim, n_years = n_years,
               seed = seed, fixed_predation = pred_rate)
}

simulate_dynamic <- function(pop, F_max, observed_recruits, packet_size,
                             n_sim = 100L, n_years = 100L, seed = 42L,
                             lh = build_life_history(pop)) {
  pm <- build_smb_predation(n_years, n_sim, observed_recruits, F_max, seed)
  result <- simulate_hbc(pop, lh = lh, n_sim = n_sim,
                         n_years = n_years, seed = seed,
                         dynamic_predation = pm, packet_size = packet_size)
  result$F_max <- F_max
  result$peak_mortality <- 1 - exp(-F_max)
  result
}

summarize_trajectory <- function(result) {
  mat <- result$N_adults
  data.frame(
    year = result$years,
    mean_adults = rowMeans(mat),
    median_adults = apply(mat, 1L, median),
    lower_025 = apply(mat, 1L, quantile, probs = 0.025, names = FALSE),
    upper_975 = apply(mat, 1L, quantile, probs = 0.975, names = FALSE),
    p_qe = result$risk$p_qe,
    p_zero_adults = result$risk$p_zero_adults,
    p_modeled_extinction = result$risk$p_modeled_extinction
  )
}

# Per-replicate mean abundance during years 51-100 in a 100-year run.
# Relativize population means against the *same type* no-predation baseline.
settled_adults <- function(result, years = 51:100) {
  if (any(!years %in% result$years)) stop("Settled years outside simulation")
  colMeans(result$N_adults[match(years, result$years), , drop = FALSE])
}

first_year_at_risk <- function(risk_by_year, model_years, target = 0.5) {
  hit <- which(risk_by_year >= target)
  if (length(hit) == 0L) NA_integer_ else as.integer(model_years[hit[1L]])
}

# Controlled decomposition: preserve reference fecundity while replacing
# ONLY the survival schedule using a candidate K/M combination.
controlled_survival_lh <- function(reference_pop, K_survival, M_survival) {
  lh <- build_life_history(reference_pop)
  rel_length <- 1 - exp(-K_survival * lh$ages)
  annual_survival <- pmin(exp(-M_survival / rel_length * 0.66), 1)
  survivorship <- numeric(length(lh$ages))
  survivorship[1L] <- 1
  survivorship[2L:length(lh$ages)] <-
    cumprod(annual_survival[1L:(length(lh$ages) - 1L)])
  lh$annual_survival <- annual_survival
  lh$survivorship <- survivorship
  lh
}


# =====================================================================
# 01. FIXED JUVENILE PREDATION (after Beverton-Holt recruitment)
# =====================================================================

pops <- hbc_populations()
# Simple audit against the manuscript's equilibrium recruitment values.
expected_R0 <- c(DNM = 2713, LCR = 27132, WGC = 49014)
stopifnot(all(vapply(pops, function(x) x$R0, numeric(1)) == expected_R0))
cat("\nCalibrated equilibrium age-1 female recruits (R0):\n")
for (pk in names(pops)) {
  cat(sprintf("  %s: N0=%d adults, R0=%d recruits, QE threshold=%d\n",
              pk, pops[[pk]]$N0, pops[[pk]]$R0, pops[[pk]]$ext_thresh))
}

fixed_risk_rows <- list()
fixed_settled_rows <- list()
fixed_trajectory_rows <- list()

for (pk in names(pops)) {
  pop <- pops[[pk]]
  for (p in FIXED_RATES) {
    cat(sprintf("01 fixed: %s, juvenile mortality = %.0f%%\n", pk, 100 * p))
    r <- simulate_fixed(pop, pred_rate = p, n_sim = N_SIM,
                        n_years = N_YEARS, seed = SEED)
    for (yr in c(20L, 100L)) {
      k <- match(yr, r$years)
      fixed_risk_rows[[length(fixed_risk_rows) + 1L]] <- data.frame(
        population = pk, fixed_mortality = p, year = yr,
        p_qe = r$risk$p_qe[k],
        p_zero_adults = r$risk$p_zero_adults[k],
        p_modeled_extinction = r$risk$p_modeled_extinction[k]
      )
    }

    # Per-simulation mean adult abundance in years 51-100 (used by 05).
    avg <- settled_adults(r)
    fixed_settled_rows[[length(fixed_settled_rows) + 1L]] <- data.frame(
      population = pk, fixed_mortality = p,
      weighted_realized_mortality = r$weighted_realized_mortality,
      settled_mean = mean(avg),
      settled_lo = unname(quantile(avg, 0.025)),
      settled_hi = unname(quantile(avg, 0.975)),
      p_qe_100 = tail(r$risk$p_qe, 1L)
    )

    # Compact trajectories for the rates shown in the manuscript figures.
    if (p %in% c(0, 0.25, 0.50, 0.75)) {
      t <- summarize_trajectory(r)
      t$population <- pk
      t$fixed_mortality <- p
      fixed_trajectory_rows[[length(fixed_trajectory_rows) + 1L]] <- t
    }
  }
}

fixed_risk <- do.call(rbind, fixed_risk_rows)
fixed_settled <- do.call(rbind, fixed_settled_rows)
fixed_trajectories <- do.call(rbind, fixed_trajectory_rows)
write.csv(fixed_risk, "results/01_fixed_risk.csv", row.names = FALSE)
write.csv(fixed_trajectories, "results/01_fixed_trajectories.csv", row.names = FALSE)
cat("01 complete: results/01_fixed_risk.csv and 01_fixed_trajectories.csv\n")
