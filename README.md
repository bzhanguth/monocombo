# monocombo

Design, simulation, and analysis of seamless clinical trials that start with
monotherapy and roll in a combination arm mid-trial (OR-Platform, Platform,
2-stage MAMS, and 3-stage MAMS) with binary endpoints and strong control of
the family-wise error rate.

## Installation

```r
# install.packages("remotes")
remotes::install_github("bzhanguth/monocombo", build_vignettes = TRUE)
```

## Main functions

- `mc_Design()`: critical values and per-arm sample size
- `mc_OC()`: operating characteristics by simulation
- `mc_Analyze()`: stage-wise decision analysis on observed data

See `vignette("monocombo-examples", package = "monocombo")` for worked examples.
