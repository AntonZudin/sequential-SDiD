
# sequential-SDiD

Sequential Synthetic Difference-in-Differences (sequential SDiD)
estimator for event studies with staggered treatment adoption.

This package estimates average treatment effect, particularly when the
parallel trends assumption fails.

## Installation

The package can be installed from github using devtools.

``` r
install.packages("devtools")
devtools::install_github("AntonZudin/sequential-SDiD")
```

## Usage

``` r
library(seq.sdid)
data(CHC)
```

The package uses 2 functions to prepare the dataset .

### 1. `to_wide`

- Converts panel from long to wide format.
- Creates Y, W and X dataframes with adoption date column being the
  first one.

### 2. `prepare_wide`

- If the level is `cohort`, aggregates on the data on the cohort level.
- If the level is `unit`, removes the adoption date column.

``` r
panel <- to_wide(CHC,
  unit = 1, time = 3, outcome = 2,
  treatment = 4, contr_covs = c(5), treat_covs = c(),
  population = 6
)

panel$X$popwt <- panel$X$popwt / 10000

panel_avg <- prepare_wide(panel, level = "cohort")
```

``` r
set.seed(42)
s2 <- estimate_s2(panel$Y[, -1], panel$W[, -1], panel$X$popwt)
est <- sequential_estimator(
  panel_avg,
  s2 = s2,
  type = "sdid")
est
```

    ## Synthetic DiD:
    ## -4.74 -5.68 -6.56 -9.44 -13.22 -13.56 -16.10 -11.43 -13.07 -12.85 -11.87 -13.72 -12.59 -10.72 -12.35 -15.61 -14.19 -12.45 -9.56 -12.38 -0.60 -2.39 -75.64 -64.37 
    ## 
    ## s2: 2728.60

``` r
se <- vcov(est, panel)
se
```

    ##  [1]  3.442402  3.517237  4.417285  4.447733  5.227896  5.278684  5.962290
    ##  [8]  6.283429  6.569504  7.342846  6.974452  9.543372  9.972621 10.263685
    ## [15] 10.968570 11.736318 11.854961 12.786417 14.031863 15.222889 14.862728
    ## [22] 16.469635 48.961334 54.116311

``` r
plot(est, se)
```

<img src="man/figures/unnamed-chunk-4-1.png" width="672" />

#### References

Dmitry Arkhangelsky, Aleksei Samkov. <b>Sequential Synthetic Difference
in Differences</b>, 2025.
\[<a href="https://arxiv.org/abs/2404.00164">arxiv</a>\]
