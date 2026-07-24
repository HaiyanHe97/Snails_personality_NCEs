# Non-consumptive effects in *Littorina littorea*

Data and R code for a study on how parasitic infection and predator cues
jointly shape antipredator behaviour in the periwinkle snail
*Littorina littorea*.

## Overview

We measured three behavioural responses to predator cues under infection,
across 60 snails × 9 trials (540 trials total):

- **Climbing height** — how high a snail climbs (continuous)
- **Time to emerge** — latency to leave the water (continuous)
- **Emergence** — whether a snail left the water at all (binary)

Each snail was tested three times under each of three cue treatments:
filtered seawater (CON), a cue from the invasive crab *Hemigrapsus
sanguineus* (INV), and a cue from the native crab *Carcinus maenas* (NAT).

## Analysis approach

Two of the three measures are limited by the apparatus and are handled
with right-censored models in `brms`:

- **Climbing height** — trials in which a snail climbed past the 30 cm
  cylinder are right-censored at 30 (Gaussian model); trials in which a
  snail did not climb (value 0) are kept as observed zeros
- **Time to emerge** — trials with no emergence within the 20 min window
  are right-censored at 20 (lognormal model)
- **Emergence** — binomial GLMM, an anchor independent of the censoring
  assumptions

Earlier non-censored models (`lmer` / `rptR`) are retained in each script
as sensitivity analyses, labelled `SENSITIVITY (OLD)`. They are not
reported in the manuscript.

Body condition is checked once, as a possible confound rather than a
result: weight at a given shell size is compared between infection groups.
This is the only part of the analysis fitted outside `brms` apart from the
emergence models and the covariate screening.

## Repository structure

| File | Purpose |
|------|---------|
| `00_merge_clean.R` | Data cleaning |
| `01_setup.R` | Load data, build censoring variables, colours, sampler settings, helpers |
| `02_climbing_height.R` | Climbing height: effects, repeatability, random slopes, distributional model |
| `03_time_to_emerge.R` | Time to emerge: full censored analysis |
| `04_emergence.R` | Emergence: binomial GLMM |
| `05_syndrome.R` | Behavioural syndrome: height × time to emerge, height × emergence |
| `06_figures.R` | Figures |
| `07_body_condition.R` | Shell size (PCA) and body condition check; Figures S1a and S1b |

Analysis scripts source `01_setup.R` and can be run independently.
Fitted models are written to `outputs/models/` so that the figure scripts
do not refit them.

## Data

- `data/raw/` — raw trial data
- `data/processed/master_snails.csv` — cleaned data produced by
  `00_merge_clean.R`; contains one row per trial, including the shell
  measurements and weight repeated across each snail's nine trials

## Requirements

R (≥ 4.2; analyses were run under 4.4.2) with `brms`, `lme4`, `lmerTest`,
`car`, `emmeans`, `rptR`, `loo`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`
and `psych`.

## Author

Haiyan He
