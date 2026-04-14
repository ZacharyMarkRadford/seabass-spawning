# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the analysis

Open `code/analysis.R` in RStudio or run it from the R console:

```r
source(here::here("code", "analysis.R"))
```

The script uses `here::here()` for all paths, so it must be run from the repo root (or via RStudio with the project open).

**Formatting:** The project uses [Air](https://posit-dev.github.io/air/) for R code formatting. In VSCode, format-on-save is enabled via `Posit.air-vscode`. To format from the terminal: `air format code/analysis.R`.

## Architecture

This is a single-script R analysis project. The data flow is:

1. **Input**: `input/raw_data.csv` — individual fish records with otolith-derived ages (`fish_age`), pelagic larval duration (`PLD`/`pld`), body size, survey location/date, and raw spawn date (as julian day offset from survey date).

2. **Data prep** (`code/analysis.R`, lines 1–35): Reads and cleans the CSV, renames columns, converts the julian-day spawn date offset to a calendar date, and derives `nursery_duration = fish_age - pelagic_duration`.

3. **Two analytical methods** to correct for sampling bias — early-spawned fish settle earlier and are underrepresented in late surveys:

   - **Method 1 — Mn optimisation** (lines 131–211): For each survey date in the Welsh sample, numerically optimises the nursery mortality rate (`Mn`) so that the mortality-weighted mean spawn date matches the overall Welsh mean. Uses `optimize()` on a squared-error objective. Output: `mn_by_survey_date`.

   - **Method 2 — Left-truncated Normal MLE** (lines 213–442): Fits a left-truncated Normal distribution to each survey-date subsample via `optim()`/BFGS, treating fish with spawn dates below the minimum observed as absent (not censored). Computes a correction factor (CF) as the reciprocal of the probability of observing a fish given the truncation point. Output: `fits_B`, `LOQ`.

4. **Reference literature**: `misc/` contains three PDFs (Fox et al. 2007, Nash & Geffen 2012, Geffen et al. 2011) providing biological context on mortality rates and the LOQ method.

## Key domain concepts

- `spawn_date_julian`: day-of-year the fish was spawned (derived from survey date minus age)
- `pelagic_duration` (PLD): days spent as a larva before settlement
- `nursery_duration`: days spent in nursery habitat after settlement (`fish_age - pelagic_duration`)
- `Mp` / `Mn`: daily instantaneous mortality rates for pelagic and nursery phases respectively
- The Wales pooled sample is treated as the "true" spawning distribution; per-survey-date estimates are corrected to match it
- Plaice mortality priors used as a starting point: `Mp = 0.079`, `Mn = 0.011`
