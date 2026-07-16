# Non-consumptive effects in *Littorina littorea*

Data and R code for a study on how parasitic infection and predator cues
jointly shape antipredator behaviour in the periwinkle snail
*Littorina littorea*.

## Overview

We measured three behavioural responses to predator cues under infection,
across 60 snails × 9 trials (540 trials total):

- **Climbing height** — how high a snail climbs (continuous)
- **Time to emerge** — latency to emerge from the shell (continuous)
- **Emergence** — whether a snail emerged at all (binary)

## Analysis approach

Two behavioural measures involve a physical ceiling, handled with
right-censored models (`brms`):

- **Time to emerge** — trials with no emergence within 20 min are
  right-censored at 20 (lognormal model)
- **Climbing height** — trials climbing past the 30 cm cylinder are
  right-censored at 30 (Gaussian model); trials that did not climb
  (value 0) are kept as observed zeros
- **Emergence** — binomial GLMM, an anchor independent of the
  censoring assumptions

Earlier non-censored models (`lmer` / `rptR`) are retained in each script
as sensitivity analyses, labelled `SENSITIVITY (OLD)`.

## Repository structure

| File | Purpose |
|------|---------|
| `00_merge_clean.R` | Data cleaning |
| `01_setup.R` | Load data, build censoring variables, colours, helpers |
| `02_climbing_height.R` | Climbing height: censored (main) + lmer/rptR (sensitivity) |
| `03_time_to_emerge.R` | Time to emerge: full censored analysis |
| `04_emergence.R` | Emergence: binomial GLMM |
| `05_syndrome.R` | Behavioural syndrome: height x time-to-emerge, height x emergence |
| `06_figures.R` | Figures |

Analysis scripts source `01_setup.R` and can be run independently.

## Data

- `data/raw/` — raw trial data
- `data/processed/` — cleaned data produced by `00_merge_clean.R`

## Requirements

R (≥ 4.2) with `brms`, `lme4`, `rptR`, and `tidyverse`.

## Author

Haiyan He
