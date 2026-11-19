# m6APrediction

## Overview

`m6APrediction` is an R package for predicting potential N6-methyladenosine (m6A)
modification sites in RNA sequences.

It uses a pre-trained Random Forest model to estimate whether a candidate RNA site
is likely to be modified. The prediction is based on several features, including GC
content, RNA type, RNA region, exon length, distance to the nearest splice junction,
evolutionary conservation, and the surrounding 5-mer sequence.

The package provides two main prediction functions: one for predicting multiple sites
from a feature table and another for predicting an individual site from its features.

## Installation

You can install the package directly from GitHub using either `devtools` or `remotes`.

Using `devtools`:

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("468636/m6APrediction")
```

Alternatively, using `remotes`:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("468636/m6APrediction")
```

## Usage

First, load the package and its pre-trained model:

```r
library(m6APrediction)

rf_fit <- load_rf_model()
```

### 1. Predict multiple RNA sites

Use `prediction_multiple()` when you have a feature table containing several candidate
sites.

The package includes an example dataset that can be used to try the function:

```r
example_df <- read.csv(
  system.file(
    "extdata",
    "m6A_input_example.csv",
    package = "m6APrediction"
  ),
  stringsAsFactors = FALSE
)

results <- prediction_multiple(rf_fit, example_df)

head(results)
```

The returned table includes two additional columns:

- `predicted_m6A_prob`: the predicted probability that a site is modified.
- `predicted_m6A_status`: the predicted classification (`Positive` or `Negative`).

### 2. Predict a single RNA site

Use `prediction_single()` when you want to predict one candidate site by providing
its feature values individually.

For example:

```r
prediction_single(
  rf_fit,
  gc_content = 0.66,
  RNA_type = "lncRNA",
  RNA_region = "intron",
  exon_length = 0,
  distance_to_junction = 10.67,
  evolutionary_conservation = 0.016,
  DNA_5mer = "GGACC"
)
```

The function returns the prediction for the individual site, including its probability
and predicted status.

## Model Performance

The Random Forest model was evaluated using the m6A feature dataset. Its performance
is summarised using ROC and precision-recall curves.

![ROC curve](man/figures/roc.png)

![Precision-Recall curve](man/figures/prc.png)

The reported performance values are:

- ROC AUC: **0.996**
- PRC AUC: **0.988**

## Features Used

The model uses the following features for prediction:

| Feature | Description |
|---|---|
| `gc_content` | GC content of the surrounding region (0–1) |
| `RNA_type` | RNA type, such as mRNA, lincRNA, lncRNA or pseudogene |
| `RNA_region` | RNA region, such as CDS, intron, 3'UTR or 5'UTR |
| `exon_length` | Exon length |
| `distance_to_junction` | Distance to the nearest splice junction |
| `evolutionary_conservation` | Evolutionary conservation score (0–1) |
| `DNA_5mer` | 5-mer sequence surrounding the candidate adenosine |

## License

This package is released under the MIT License.
