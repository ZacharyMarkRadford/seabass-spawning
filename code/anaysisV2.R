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
      # L <- 10
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
    method2_results |>
      group_by(country) |>
      summarise(
        mean_M2 = weighted.mean(mean_spawn_corrected, w = n),
        .groups = "drop"
      ),
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
  )

# Three-method comparison plot for Wales
comparison_wales |>
  tidyr::pivot_longer(
    c(mean_spawn_raw, corrected_M2),
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
      corrected_M3 = "#009E73"
    )
  ) +
  labs(
    title = "Three correction methods vs. raw mean spawn date (Wales)",
    subtitle = "Dashed line = Wales pooled mean (target ≈ Julian day 125)",
    x = "Survey day of year",
    y = "Mean spawn date (Julian day)",
    colour = "Method"
  ) +
  theme_classic(base_size = 14)
