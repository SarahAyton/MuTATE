# MuTATE

Collection of functions to recursively partition datasets on binary splits across multiple targets.

The goal of MuTATE is to create automated, explainable, and comprehensive models across multiple dependent variables of interest.

It provides a collection of functions to recursively partition data on binary splits across multiple targets having different dependent variable types. This overcomes single-target limitations of traditional decision trees, without losing model interpretability, and while handling continuous, categorical, count, and survival outcome variables. This suite of functions also includes a number of parameters for model customization, dependent variable weights, parameter tuning functions, and visualization tools.

## Installation

You can install the development version of MuTATE from [GitHub](https://github.com/SarahAyton/MuTATE) with:

``` r
# install.packages("devtools")
devtools::install_github("SarahAyton/MuTATE")
```

## Usage

```r
library(MuTATE)
data(mutate_example)

features     <- c("age", "sex", "biomarker")
outcomes     <- c("response", "tumor_size", "ae_count", "OS_definition_time_status")
outcome_defs <- c("Cat", "Cont", "Count", "Surv")

tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
               depth = 2, nodesize = 30)
PlotTree(tree)
```

See `vignette("mutate-getting-started", package = "MuTATE")` for a full walkthrough covering fitting, pruning, prediction on new data, and visualization.

## License

GPL (>= 3) © Sarah Ayton
