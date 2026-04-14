library(tidyverse)
library(here)
library(janitor)


# Read in the data
spawning_data <- read_csv(here("input", "raw_data.csv")) |>
  clean_names() |>
  rename(
    spawn_date_julian = spawn_date
  ) |>
  mutate(
    date = as.Date(date, format = "%d/%m/%Y"),
    # Convery julian date to standard date format
    spawn_date = date - spawn_date_julian
  ) |>
  select(
    lab_id,
    country,
    region,
    length_mm,
    wet_wt_g,
    survey_date = date,
    parent_site,
    reach,
    site_code = site_code_x,
    survey_lat,
    survey_long,
    temp_f,
    fish_age,
    pelagic_duration = pld,
    spawn_date,
    spawn_date_julian,
  ) |>
  mutate(nursery_duration = fish_age - pelagic_duration)

# Select welsh data ------------------------------------------------------

welsh_data <- spawning_data |>
  filter(country == "Wales")

# ── Exploratory visualisation ───────────────────────────────────────────
# Show the bias before any correction: in later surveys the raw mean spawn
# date drifts later because early-spawned fish have been dying longer.

welsh_data |>
  group_by(survey_date) |>
  summarise(
    mean_spawn = mean(spawn_date_julian, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) |>
  ggplot(aes(x = lubridate::yday(survey_date), y = mean_spawn)) +
  geom_point(aes(size = n)) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    title = "Bias: raw mean spawn date drifts later with survey date (Wales)",
    x = "Survey day of year",
    y = "Raw mean spawn date (Julian day)",
    size = "n fish"
  ) +
  theme_classic(base_size = 14)


# ── Wales reference distribution ─────────────────────────────────────────
# All Wales fish pooled = the "true" spawning distribution we correct toward.
# Simple Normal fit (mean, sd) is sufficient for n = 138.

mu_wales <- mean(welsh_data$spawn_date_julian, na.rm = TRUE)
sd_wales <- sd(welsh_data$spawn_date_julian, na.rm = TRUE)

# Global median PLD used to define the per-survey truncation point L.
# L = survey_day - median_pld is the latest possible spawn date for a fish
# that JUST settled (nursery_duration ≈ 0) on the survey day.  Fish spawned
# before L had longer to die and are progressively absent from the sample.
median_pld <- median(spawning_data$pelagic_duration, na.rm = TRUE)

message(sprintf(
  "Wales reference: mu = %.1f (≈ Julian day), sd = %.1f | median PLD = %.0f days",
  mu_wales,
  sd_wales,
  median_pld
))


# ── Method 1: Mortality-weighted mean ───────────────────────────────────
# Each fish is upweighted by 1/survival so that early-spawned (more-depleted)
# fish count more.  W = 1 / (Sp * Sn) where:
#   Sp = exp(-Mp * pelagic_duration)   pelagic-phase survival
#   Sn = exp(-Mn * nursery_duration)   nursery-phase survival
#
# Mp is fixed at the plaice prior (0.079).  Mn is *optimised per survey date*
# on the Wales data so that the weighted mean matches the Wales pooled mean,
# then applied to every survey date in the full dataset via a smooth model.

Mp_fixed <- 0.079

# Helper: find Mn that makes the mortality-weighted mean equal mu_true
fit_Mn <- function(df_survey, mu_true, Mp = Mp_fixed, lower = 0, upper = 1) {
  obj <- function(Mn) {
    Sn <- exp(-Mn * df_survey$nursery_duration)
    Sp <- exp(-Mp * df_survey$pelagic_duration)
    W <- 1 / (Sp * Sn)
    mu_hat <- weighted.mean(df_survey$spawn_date_julian, w = W, na.rm = TRUE)
    (mu_hat - mu_true)^2
  }
  optimize(obj, interval = c(lower, upper))$minimum
}

# Apply per Wales survey date (≥ 3 fish)
wales_mn_fits <- welsh_data |>
  group_by(survey_date) |>
  filter(n() >= 3) |>
  group_modify(
    ~ {
      mn_hat <- fit_Mn(.x, mu_true = mu_wales)
      tibble(
        survey_day = lubridate::yday(.y$survey_date),
        Mn_hat = mn_hat,
        n = nrow(.x)
      )
    }
  ) |>
  ungroup()

# Smooth Mn_hat over survey day so we can predict for any survey date
# (including non-Wales countries).  A simple linear model on log(Mn) works
# well for the range of data; adjust if the relationship looks non-linear.
mn_model <- lm(log(Mn_hat) ~ survey_day, data = wales_mn_fits)

# Convenience function: predicted Mn for a given survey day-of-year
predict_Mn <- function(survey_day) {
  exp(predict(mn_model, newdata = data.frame(survey_day = survey_day)))
}

# Compute corrected weighted mean per survey date for all countries
method1_results <- spawning_data |>
  mutate(survey_day = lubridate::yday(survey_date)) |>
  group_by(country, survey_date, survey_day) |>
  group_modify(
    ~ {
      Mn_use <- predict_Mn(.y$survey_day)
      Sp <- exp(-Mp_fixed * .x$pelagic_duration)
      Sn <- exp(-Mn_use * .x$nursery_duration)
      W <- 1 / (Sp * Sn)
      tibble(
        n = nrow(.x),
        Mn_used = Mn_use,
        mean_spawn_raw = mean(.x$spawn_date_julian, na.rm = TRUE),
        mean_spawn_corrected = weighted.mean(
          .x$spawn_date_julian,
          w = W,
          na.rm = TRUE
        )
      )
    }
  ) |>
  ungroup()

# Validation plot: Wales corrected means should be flat near mu_wales
method1_results |>
  filter(country == "Wales") |>
  tidyr::pivot_longer(
    c(mean_spawn_raw, mean_spawn_corrected),
    names_to = "type",
    values_to = "mean_spawn"
  ) |>
  ggplot(aes(x = survey_day, y = mean_spawn, colour = type)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = mu_wales, linetype = "dashed") +
  labs(
    title = "Method 1 validation (Wales): corrected means should flatten to dashed line",
    x = "Survey day of year",
    y = "Mean spawn date (Julian day)"
  ) +
  theme_classic(base_size = 14)


# ── Method 2: LOQ / Left-truncated Normal MLE ───────────────────────────
# The LOQ (Limit of Quantification) is the *minimum detectable spawn date*
# for a given survey.  For a survey on day D with median PLD days of larval
# life, a fish that just settled has spawn_date ≈ D - median_pld — so fish
# spawned earlier than L = D - median_pld are progressively absent.
#
# Treating absence as *left-truncation* (not censoring), we fit a Normal
# distribution via MLE to the observed spawn dates, adjusting for the fact
# that observations below L are absent.  The truncated-Normal likelihood is:
#
#   log L = Σ log φ(Sᵢ | μ, σ)  −  n · log(1 − Φ((L − μ)/σ))
#
# This IS "optimising the CDF": optim() maximises this log-likelihood, which
# depends on pnorm() (the CDF) to compute the truncation survival probability.

# LOQ_nd: minimum detectable nursery duration (days).  Set to 0 (conservative).
# Increase to e.g. 7 or 14 if very young fish can't be reliably aged.
LOQ_nd <- 0

# Truncation point as a function of survey day
L_for_survey <- function(survey_day, pld_med = median_pld, loq_nd = LOQ_nd) {
  survey_day - pld_med - loq_nd
}

# Negative log-likelihood for a left-truncated Normal
fit_truncnorm_mle <- function(S_obs, L, init_mu = NULL, init_sd = NULL) {
  S_obs <- S_obs[is.finite(S_obs)]

  # Safe initialisations: fall back to sd_wales if sd is NA/0
  if (is.null(init_mu)) {
    init_mu <- mean(S_obs)
  }
  if (is.null(init_sd) || !is.finite(init_sd) || init_sd <= 0) {
    init_sd <- sd_wales
  }

  # L must be strictly below the minimum observation, otherwise the
  # truncated-normal is undefined for this sample.  Shift it down if needed.
  L <- min(L, min(S_obs) - 0.5)

  nll <- function(par) {
    mu <- par[1]
    sigma <- exp(par[2]) # always > 0
    p_survive <- 1 - pnorm(L, mean = mu, sd = sigma) # P(S ≥ L)
    if (!is.finite(p_survive) || p_survive < 1e-10) {
      return(1e15)
    }
    ll <- sum(dnorm(S_obs, mean = mu, sd = sigma, log = TRUE)) -
      length(S_obs) * log(p_survive)
    -ll
  }

  opt <- tryCatch(
    optim(
      par = c(init_mu, log(init_sd)),
      fn = nll,
      method = "BFGS"
    ),
    error = function(e) list(par = c(init_mu, log(init_sd)), convergence = 1)
  )

  list(
    mu = opt$par[1],
    sd = exp(opt$par[2]),
    L = L,
    converged = opt$convergence == 0
  )
}

# Apply to each Wales survey date
fits_m2 <- welsh_data |>
  mutate(survey_day = lubridate::yday(survey_date)) |>
  group_by(survey_date) |>
  filter(n() >= 3) |>
  group_modify(
    ~ {
      S_obs <- .x$spawn_date_julian
      survey_day <- lubridate::yday(.y$survey_date)
      # Per-survey diagnostic: use observed minimum as L so it is always
      # below every observation in this sample.
      L <- min(S_obs, na.rm = TRUE) - 0.5
      # L CHANGER
      # If you wanjt to apply a global L, comment out line above and uncomment line below
      # Large L = more correction, smaller L = less correction
      # L <- 87
      fit <- fit_truncnorm_mle(
        S_obs,
        L = L,
        init_sd = sd(S_obs, na.rm = TRUE)
      )
      p_obs <- 1 - pnorm(L, mean = fit$mu, sd = fit$sd)
      tibble(
        survey_day = survey_day,
        n = length(S_obs),
        L = L,
        mu_hat = fit$mu,
        sd_hat = fit$sd,
        CF_M2 = 1 / p_obs,
        converged = fit$converged
      )
    }
  ) |>
  ungroup()

fits_m2

# Validation: mu_hat per survey date should cluster near mu_wales
fits_m2 |>
  ggplot(aes(x = survey_day, y = mu_hat)) +
  geom_point(aes(size = n)) +
  geom_hline(yintercept = mu_wales, linetype = "dashed", colour = "steelblue") +
  geom_hline(
    yintercept = mean(welsh_data$spawn_date_julian),
    linetype = "dotted"
  ) +
  labs(
    title = "Method 2 (truncated-Normal MLE): recovered mu per Wales survey date",
    subtitle = "Dashed = Wales pooled mean (target)",
    x = "Survey day of year",
    y = "Estimated true mean spawn date"
  ) +
  theme_classic(base_size = 14)

# Apply to all countries.  L for each survey date is the *observed minimum
# spawn date* in that sample minus a small offset (data-driven truncation
# point).  This is more reliable than a theoretical formula because:
#  - it never exceeds any observation (preventing degenerate likelihoods), and
#  - it naturally increases for later surveys as early fish are depleted by
#    mortality — exactly capturing the bias we want to correct.
#
# The mean correction uses the inverse Mills ratio (λ) to compute the
# expected mean of the *observed* (truncated) sample under the Wales
# reference distribution, then shifts the raw observed mean by the same
# amount:
#
#   mean_trunc  = μ_W + σ_W · φ(α)/Φ̄(α)   where α = (L − μ_W)/σ_W
#   mean_corrected = mean_raw + (μ_W − mean_trunc)
#
# Positive correction → raw mean too high (later survey, early fish depleted)
# Negative correction → raw mean too low  (unlikely in this dataset)

# Compute L per survey date and then the truncation correction in one step
method2_results <- spawning_data |>
  group_by(country, survey_date) |>
  mutate(
    survey_day = lubridate::yday(survey_date),
    # Data-driven L: observed minimum for this survey, shifted just below it
    L = min(spawn_date_julian, na.rm = TRUE) - 0.5,
    # L CHANGER
    # If you wanjt to apply a global L, comment out line above and uncomment line below
    # Large L = more correction, smaller L = less correction
    # L <- 87
    alpha = (L - mu_wales) / sd_wales,
    p_obs = 1 - pnorm(alpha), # P(S ≥ L | Wales distribution)
    CF_M2 = 1 / p_obs,
    mean_trunc = mu_wales + sd_wales * dnorm(alpha) / (1 - pnorm(alpha))
  ) |>
  summarise(
    survey_day = first(survey_day),
    n = n(),
    L = first(L),
    CF_M2 = first(CF_M2),
    mean_spawn_raw = mean(spawn_date_julian, na.rm = TRUE),
    mean_spawn_corrected = mean_spawn_raw + (mu_wales - first(mean_trunc)),
    .groups = "drop"
  )


# ── Method 2.1: Combined true spawning distribution from per-survey MLE ─────
# Method 2 only corrects the *mean* of each survey sample.  Method 2.1 goes
# further: it runs the truncated-Normal MLE for every survey date across ALL
# countries to estimate the full true (mu_hat, sd_hat) for that survey, then
# combines the per-survey distributions into a single country-level spawning
# distribution as a mixture of Normals weighted by sample size.
#
# Why this is useful:
#   - It recovers a full distribution, not just a corrected mean
#   - Different surveys cover different parts of the spawning season, so
#     combining them weights each part of the season by how many fish it
#     contributed
#   - The mixture can be asymmetric or multi-modal in ways a simple mean hides

# Step 1: fit truncated Normal per survey date for every country (≥ 3 fish)
fits_m2_all <- spawning_data |>
  group_by(country, survey_date) |>
  filter(n() >= 3) |>
  group_modify(
    ~ {
      S_obs <- .x$spawn_date_julian
      L <- min(S_obs, na.rm = TRUE) - 0.5
      fit <- fit_truncnorm_mle(S_obs, L = L, init_sd = sd(S_obs, na.rm = TRUE))
      tibble(
        survey_day = lubridate::yday(.y$survey_date),
        n = length(S_obs),
        L = L,
        mu_hat = fit$mu,
        sd_hat = fit$sd,
        converged = fit$converged
      )
    }
  ) |>
  ungroup()

# Step 2: evaluate the mixture density over a spawn-date grid for each country.
# Each survey contributes dnorm(x, mu_hat, sd_hat) weighted by n / n_country_total.
# Using outer() for a vectorised calculation — no loops.

spawn_grid <- seq(-50, 250, by = 1)

mixture_density_m2_1 <- fits_m2_all |>
  group_by(country) |>
  group_modify(
    ~ {
      w <- .x$n / sum(.x$n)
      # density_matrix rows = spawn_grid values, columns = survey components
      density_matrix <- outer(
        spawn_grid,
        seq_along(.x$mu_hat),
        function(d, j) dnorm(d, mean = .x$mu_hat[j], sd = .x$sd_hat[j])
      )
      tibble(
        spawn_date_julian = spawn_grid,
        density = as.vector(density_matrix %*% w)
      )
    }
  ) |>
  ungroup()

# Step 3: derive mean and sd of each country's mixture distribution.
# Mixture mean  = weighted mean of component means
# Mixture var   = weighted mean of (sigma_i^2 + mu_i^2) − (mixture mean)^2
#                 (law of total variance)

mixture_summary_m2_1 <- fits_m2_all |>
  group_by(country) |>
  summarise(
    n_surveys = n(),
    n_total = sum(n),
    mean_M2_1 = weighted.mean(mu_hat, w = n),
    var_M2_1 = weighted.mean(sd_hat^2 + mu_hat^2, w = n) -
      weighted.mean(mu_hat, w = n)^2,
    .groups = "drop"
  ) |>
  mutate(sd_M2_1 = sqrt(var_M2_1))

print(mixture_summary_m2_1)

# Plot: overlay corrected distributions for all countries plus the raw Wales
# reference, so differences in spawning timing become visually apparent
raw_wales_ref <- tibble(
  spawn_date_julian = spawn_grid,
  density = dnorm(spawn_grid, mean = mu_wales, sd = sd_wales),
  country = "Wales (raw reference)"
)

mixture_density_m2_1 |>
  bind_rows(raw_wales_ref) |>
  ggplot(aes(x = spawn_date_julian, y = density, colour = country)) +
  geom_line(linewidth = 0.9) +
  geom_vline(
    data = mixture_summary_m2_1,
    aes(xintercept = mean_M2_1, colour = country),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  labs(
    title = "Method 2.1: corrected spawning distribution per country",
    subtitle = "Mixture of per-survey truncated-Normal fits; dashed = mixture mean",
    x = "Spawn date (Julian day of year)",
    y = "Probability density",
    colour = "Country"
  ) +
  theme_classic(base_size = 14)


# ── Comparison and summary across countries ──────────────────────────────
# Combine corrected means from all three methods and compare.

summary_by_country <- spawning_data |>
  group_by(country) |>
  summarise(
    n = n(),
    mean_spawn_raw = mean(spawn_date_julian, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    method1_results |>
      group_by(country) |>
      summarise(
        mean_M1 = weighted.mean(mean_spawn_corrected, w = n),
        .groups = "drop"
      ),
    by = "country"
  ) |>
  left_join(
    method2_results |>
      group_by(country) |>
      summarise(
        mean_M2 = weighted.mean(mean_spawn_corrected, w = n),
        .groups = "drop"
      ),
    by = "country"
  ) |>
  left_join(
    mixture_summary_m2_1 |> select(country, mean_M2_1),
    by = "country"
  )

print(summary_by_country)

# Per-survey-date comparison for Wales (validation)
comparison_wales <- method2_results |>
  filter(country == "Wales") |>
  select(
    survey_date,
    survey_day,
    n,
    mean_spawn_raw,
    corrected_M2 = mean_spawn_corrected
  ) |>
  left_join(
    method1_results |>
      filter(country == "Wales") |>
      select(survey_date, corrected_M1 = mean_spawn_corrected),
    by = "survey_date"
  ) |>
  left_join(
    fits_m2_all |>
      filter(country == "Wales") |>
      select(survey_date, corrected_M2_1 = mu_hat),
    by = "survey_date"
  )

# Comparison plot for Wales (all methods)
comparison_wales |>
  tidyr::pivot_longer(
    c(mean_spawn_raw, corrected_M1, corrected_M2, corrected_M2_1),
    names_to = "method",
    values_to = "mean_spawn"
  ) |>
  ggplot(aes(x = survey_day, y = mean_spawn, colour = method)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = mu_wales, linetype = "dashed") +
  scale_colour_manual(
    values = c(
      mean_spawn_raw = "grey50",
      corrected_M1 = "#E69F00",
      corrected_M2 = "#56B4E9",
      corrected_M2_1 = "#CC79A7"
    )
  ) +
  labs(
    title = "Correction methods vs. raw mean spawn date (Wales)",
    subtitle = "Dashed line = Wales pooled mean (target ≈ Julian day 125)",
    x = "Survey day of year",
    y = "Mean spawn date (Julian day)",
    colour = "Method"
  ) +
  theme_classic(base_size = 14)


# ── Overall combined mean spawn date (all countries except Belgium) ───────
# Pool England, France, and Wales together and compute a single overall
# corrected mean under each method, weighted by sample size.
#
# Method 1  : weighted mean of per-survey corrected means
# Method 2  : weighted mean of per-survey corrected means
# Method 2.1: mixture mean of all per-survey truncated-Normal fits
#             (uses the same law-of-total-variance pooling as above, but
#              now across countries rather than within a single country)

overall_M1 <- method1_results |>
  filter(country != "Belgium") |>
  summarise(
    n_total = sum(n),
    overall_mean = weighted.mean(mean_spawn_corrected, w = n)
  )

overall_M2 <- method2_results |>
  filter(country != "Belgium") |>
  summarise(
    n_total = sum(n),
    overall_mean = weighted.mean(mean_spawn_corrected, w = n)
  )

overall_M2_1 <- fits_m2_all |>
  filter(country != "Belgium") |>
  summarise(
    n_total = sum(n),
    overall_mean = weighted.mean(mu_hat, w = n),
    overall_var = weighted.mean(sd_hat^2 + mu_hat^2, w = n) -
      weighted.mean(mu_hat, w = n)^2
  ) |>
  mutate(overall_sd = sqrt(overall_var))

overall_summary <- tibble(
  method = c("Raw (excl. Belgium)", "Method 1", "Method 2", "Method 2.1"),
  n = c(
    sum(spawning_data$country != "Belgium"),
    overall_M1$n_total,
    overall_M2$n_total,
    overall_M2_1$n_total
  ),
  overall_mean = c(
    mean(
      spawning_data$spawn_date_julian[spawning_data$country != "Belgium"],
      na.rm = TRUE
    ),
    overall_M1$overall_mean,
    overall_M2$overall_mean,
    overall_M2_1$overall_mean
  )
)

print(overall_summary)
