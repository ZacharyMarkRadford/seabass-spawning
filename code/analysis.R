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


# Quick plot -------------------------------------------------------------

spawning_data |>
  mutate(month = month(spawn_date)) |>
  ggplot(aes(x = month, fill = country)) +
  geom_bar() +
  labs(
    title = "Age Distribution of Welsh Samples",
    x = "SpawnDate",
    y = "Count"
  ) +
  theme_minimal() +
  facet_wrap(~country) +
  xlim(1, 12)


spawning_data |>
  filter(country == "Wales") |>
  group_by(survey_date) |>
  summarise(n_fish = n_distinct(lab_id)) |>
  ggplot(aes(x = lubridate::yday(survey_date), y = n_fish)) +
  geom_line()

# Begin  analysis --------------------------------------------------------

# Merge all the welsh samples to gety an overall distribution for wales
# Call that the true distribution
# Get a good pelagic phase mortality from Jenny
# Try and fit the nursery mortality to match the true distribution for each of the sampling dates
# Optimise the mean

# Pull welsh data

welsh_data <- spawning_data |>
  filter(country == "Wales")

# Calculate % that survive the pelagic duration and the weight for each fish
welsh_data <- welsh_data |>
  mutate(
    Mp_fixed = 0.50, # Daily mortality rate for pelagic duration
    survival_pelagic = exp(-Mp_fixed * pelagic_duration) # % that survive the pelagic duration
  )


# Calculate % that survive the pelagic duration and the weight for each fish
welsh_data <- welsh_data |>
  mutate(
    Mp_fixed = 0.50, # Daily mortality rate for pelagic duration
    survival_pelagic = exp(-Mp_fixed * pelagic_duration) # % that survive the pelagic duration
  )

# Try using plaice Mn and Mp to see how much it shifts the mean spawn date

# Mp 0.079
# Mn 0.011

plaice_weights <- welsh_data |>
  mutate(
    Mp_fixed = 0.079, # Daily mortality rate for pelagic duration
    Mn_fixed = 0.011, # Daily mortality rate for nursery duration

    Sp = exp(-Mp_fixed * pelagic_duration), # Pelagic survival
    Sn = exp(-Mn_fixed * nursery_duration), # Nursery survival
    St = Sp * Sn, # Total survival
    W = 1 / St, # Weights are inverse of survival
  )

plaice_weights |>
  group_by(survey_date) |>
  summarise(
    mean_spawn_date_raw = mean(spawn_date_julian, na.rm = TRUE),
    mean_spawn_date_weighted = weighted.mean(
      spawn_date_julian,
      w = W,
      na.rm = TRUE
    ),
    n = n_distinct(lab_id),
    .groups = "drop"
  )

plaice_weights |>
  group_by(country) |>
  summarise(
    mean_spawn_date_raw = mean(spawn_date_julian, na.rm = TRUE),
    mean_spawn_date_weighted = weighted.mean(
      spawn_date_julian,
      w = W,
      na.rm = TRUE
    ),
    n = n_distinct(lab_id),
    .groups = "drop"
  )

# Optimise Mn ------------------------------------------------------------

# Optimise Mn
true_means <- welsh_data |>
  group_by(country) |>
  summarise(
    mu_true = mean(spawn_date_julian, na.rm = TRUE),
    sd_true = sd(spawn_date_julian, na.rm = TRUE),
    .groups = "drop"
  )

# Opimise for each spawning date in Welsh sample to fit mean to this 125 average spawn date
# Get correction factors time since first settlement

# Remove survey dates that have only 0 nursery duration fish as they will not provide a good estimate of the mean spawn date after weighting

Mn_optimisation_data <- welsh_data

# filter(nursery_duration > 0)

# Remove survey dates with less than 3 fish as they will not provide a good estimate of the mean spawn date
Mn_optimisation_data <- Mn_optimisation_data |>
  group_by(survey_date) |>
  filter(n_distinct(lab_id) >= 3) |>
  ungroup()

Mn_optimisation_data |>
  group_by(survey_date) |>
  summarise(n = n_distinct(lab_id))

# Create the optimisation function ---------------------------------------

fit_Mn <- function(df_year, mu_true, Mp_fixed = 0.079, lower = 0, upper = 1) {
  obj <- function(Mn) {
    # combined weights: inverse survival across pelagic + nursery
    Sn <- exp(-Mn * df_year$nursery_duration) # Nursery survival
    Sp <- exp(-Mp_fixed * df_year$pelagic_duration) # Pelagic survival
    St = Sn * Sp # Total survival
    W = 1 / St # Weights are inverse of survival
    mu_hat <- weighted.mean(df_year$spawn_date_julian, w = W, na.rm = TRUE) # Weighted mean spawn date
    (mu_hat - mu_true)^2 # Objective is to minimize squared difference from true mean
  }

  opt <- optimize(obj, interval = c(lower, upper))
  opt$minimum
}

mn_by_survey_date <- Mn_optimisation_data |>
  # Add in ther tue mean spawn date for wales to compare to
  left_join(true_means, by = "country") |>
  mutate(Mp_fixed = 0.079) |>
  # For each survey date, optimise the Mn to fit the mean spawn date to the true mean spawn date for wales
  group_by(survey_date) |>
  mutate(
    Mn_hat = fit_Mn(
      pick(everything()),
      mu_true = mu_true[1],
      Mp_fixed = Mp_fixed
    )
  ) |>
  ungroup() |>
  mutate(
    Sp = exp(-Mp_fixed * pelagic_duration), # Pelagic survival
    Sn = exp(-Mn_hat * nursery_duration), # Nursery survival
    St = Sp * Sn, # Total survival
    W = 1 / St # Weights are inverse of survival
  ) |>
  group_by(survey_date) |>
  summarise(
    mean_spawn_date_raw = mean(spawn_date_julian, na.rm = TRUE),
    mean_spawn_date_weighted = weighted.mean(
      spawn_date_julian,
      w = 1 / exp(-Mn_hat * nursery_duration),
      na.rm = TRUE
    ),
    mean_spawn_date_true = mu_true[1],
    n = n_distinct(lab_id),
    .groups = "drop"
  )

mn_by_survey_date

# Method 2: LOQ ----------------------------------------------------------

# OTher way - Fittong LOQ method
# Need to look at shape of distribution
# Fiut the distribution, then truncate
# Know true mean and variance
# Can look at impact of truncation on other samples

# Examine current distribution -------------------------------------------

welsh_data |>
  ggplot(aes(x = spawn_date_julian)) +
  geom_histogram(binwidth = 7) +
  theme_minimal()

# Define current distribution --------------------------------------------
mu_raw <- mean(welsh_data$spawn_date_julian)
sd_raw <- sd(welsh_data$spawn_date_julian)

theoretical_spawning_dist <-
  expand_grid(spawn_date_julian = 1:365, mean = mu_raw, sd = sd_raw) |>
  mutate(probability = dnorm(spawn_date_julian, mean, sd))

observed_spawning_dist <- welsh_data |>
  group_by(survey_date) |>
  summarise(
    n = n(),
    mean = mean(spawn_date_julian),
    sd = sd(spawn_date_julian),
    .groups = "drop"
  ) |>
  filter(!n < 6) |>
  expand_grid(spawn_date_julian = 1:365) |>
  mutate(probability = dnorm(spawn_date_julian, mean, sd))

spawn_date_distributions <- theoretical_spawning_dist |>
  mutate(survey_date = "Overall") |>
  # select(survey_date, spawn_date_julian, probability) |>
  bind_rows(
    observed_spawning_dist |>
      select(-n) |>
      mutate(survey_date = as.character(survey_date))
  )

ggplot() +
  geom_bar(
    data = spawn_date_distributions |> filter(survey_date == "Overall"),
    aes(
      x = spawn_date_julian,
      y = probability,
      fill = survey_date,
      colour = survey_date
    ),
    stat = "identity"
  ) +
  theme_minimal()

observed_spawning_dist |>
  select(survey_date, mean, sd) |>
  distinct() |>
  mutate(survey_date_julian = yday(survey_date)) |>
  ggplot(aes(x = survey_date_julian, y = mean)) +
  geom_point() +
  geom_line() +
  geom_smooth(se = FALSE)


# Method B: reconstruct the *true* spawning distribution from late samples
# by fitting a left-truncated Normal distribution (MLE).
#
# Interpretation:
#   - We assume true spawning dates S ~ Normal(mu, sd).
#   - For a late survey date, we only observe fish with S >= L (early spawners are missing).
#   - That is left-*truncation* (not censoring): missing values are absent, not recorded as "< L".
#   - The MLE adjusts for truncation via the factor P(S >= L) = 1 - Phi((L - mu)/sd).

fit_truncnorm_mle <- function(
  S_obs,
  L,
  init_mu = mean(S_obs),
  init_sd = sd(S_obs)
) {
  # Keep only finite observations
  S_obs <- S_obs[is.finite(S_obs)]

  # Negative log-likelihood (NLL) for a left-truncated Normal:
  #   log L = sum(log dnorm(S_i | mu, sd)) - n * log(1 - pnorm(L | mu, sd))
  # The second term corrects for the fact that values below L are *not present* in the sample.
  nll <- function(par) {
    mu <- par[1]
    sigma <- exp(par[2]) # optimise log(sigma) so sigma is always > 0

    # Truncation survival probability: P(S >= L) under Normal(mu, sigma)
    p_survive <- 1 - pnorm(L, mean = mu, sd = sigma)

    # If p_survive is ~0, likelihood is undefined (would imply essentially no probability of observing the sample)
    if (!is.finite(p_survive) || p_survive <= 0) {
      return(Inf)
    }

    # Log-likelihood = standard normal log-density contributions minus truncation penalty
    ll <- sum(dnorm(S_obs, mean = mu, sd = sigma, log = TRUE)) -
      length(S_obs) * log(p_survive)

    # optim() minimises, so return negative log-likelihood
    -ll
  }

  # Run optimisation to find mu and sigma that maximise the truncated-normal likelihood.
  # Parameterisation: par = c(mu, log(sigma)).
  opt <- optim(
    par = c(init_mu, log(init_sd)),
    fn = nll,
    method = "BFGS"
  )

  # Extract fitted parameters
  mu_hat <- opt$par[1]
  sd_hat <- exp(opt$par[2])

  # Return fitted mu, sd and the truncation point used
  list(mu = mu_hat, sd = sd_hat, L = L, opt = opt)
}

# Apply per survey_date (each survey date has its own truncated sample)
# L is the left-truncation point for that sample; simplest choice is min observed spawn DOY.
# (Biologically, L could instead be derived from "time since first settlement" if you have that rule.)
fits_B <- spawning_data |>
  filter(country == "Wales") |>
  mutate(S = spawn_date_julian, L = min(S, na.rm = TRUE)) |>
  group_by(survey_date) |>
  filter(n() >= 3) |>
  group_modify(
    ~ {
      S_obs <- .x$S
      L <- unique(.x$L) # operational truncation point for this survey-date sample

      fit <- fit_truncnorm_mle(S_obs, L = L)

      tibble(
        n = length(S_obs), # number of fish contributing to this survey-date fit
        L = L, # truncation point used
        mu_hat = fit$mu, # estimated mean of the *untruncated* (true) spawning distribution
        sd_hat = fit$sd # estimated SD of the *untruncated* (true) spawning distribution
      )
    }
  ) |>
  ungroup()

fits_B

# Plot for each survey date ----------------------------------------------

raw_params <- welsh_data |>
  group_by(survey_date) |>
  summarise(
    mu_hat = mean(spawn_date_julian, na.rm = TRUE),
    sd_hat = sd(spawn_date_julian, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(method = "raw")

hist_data <- fits_B |>
  select(survey_date, mu_hat, sd_hat) |>
  mutate(method = "truncnorm_mle") |>
  bind_rows(raw_params) |>
  group_by(survey_date) |>
  filter(n_distinct(method) == 2) |> # Only keep survey dates where we have both raw and fitted estimates
  group_by(method, survey_date) |>
  expand_grid(day_of_year = seq(1, 365, by = 1)) |>
  mutate(
    prob = dnorm(day_of_year, mean = mu_hat, sd = sd_hat),
  )

# Add in the true mean spawn date for wales to compare to
true_p <- true_means |>
  expand_grid(day_of_year = seq(1, 365, by = 1)) |>
  mutate(
    prob = dnorm(day_of_year, mean = mu_true, sd = sd_true),
    method = "true"
  ) |>
  expand_grid(survey_date = unique(hist_data$survey_date))


hist_data <- hist_data |>
  bind_rows(true_p)

ggplot(hist_data, aes(x = day_of_year, y = prob, colour = method)) +
  geom_line() +
  facet_wrap(~survey_date, scales = "free_y") +
  labs(
    title = "Fitted true spawning distribution per survey date",
    x = "Day of Year",
    y = "Probability Density"
  ) +
  theme_classic(base_size = 18) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# LOQ corrections --------------------------------------------------------

# Start with mean of means
# Look at tiem between fiorst settlement we haev recorded and sampling dat

fits_B |>
  ggplot(aes(x = ))


welsh_cf_overall <- fits_B |>
  mutate(
    p_obs = 1 - pnorm(L, mean = mu_hat, sd = sd_hat),
    CF = 1 / p_obs
  ) |>
  summarise(
    CF_overall_n_weighted = weighted.mean(CF, w = n),
    CF_overall_simple = mean(CF),
    .groups = "drop"
  )

welsh_cf_overall

LOQ <- fits_B |>
  summarise(
    mu_hat = weighted.mean(mu_hat, w = n),
    sd_hat = weighted.mean(sd_hat, w = n),
    L = weighted.mean(L, w = n)
  ) |>
  mutate(
    alpha = (L - mu_hat) / sd_hat,
    p_obs = 1 - pnorm(alpha),
    CF_total = 1 / p_obs
  )
LOQ
