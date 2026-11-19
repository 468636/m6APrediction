## m6APrediction: Random-forest prediction of m6A RNA modification sites
##
## This file contains the package source code only. All testing and
## installation commands are kept outside the package directory.

#' @import randomForest
#' @importFrom stats predict
NULL


#' Predict m6A status for a feature table (multiple sites)
#'
#' Apply the package random forest model to a data frame of m6A-site
#' features, returning the input table augmented with the predicted
#' m6A probability and binary status.
#'
#' @param ml_fit A fitted `randomForest` classification model. The model
#'   shipped in the package is loaded via
#'   `readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))`.
#' @param feature_df A data frame with columns: `gc_content`, `RNA_type`,
#'   `RNA_region`, `exon_length`, `distance_to_junction`,
#'   `evolutionary_conservation`, `DNA_5mer`.
#' @param positive_threshold Numeric threshold in `[0, 1]` on the predicted
#'   positive-class probability above which a site is called `Positive`.
#'   Defaults to `0.5`.
#'
#' @return The input `feature_df` with two added columns:
#'   `predicted_m6A_prob` (numeric, probability of being modified) and
#'   `predicted_m6A_status` (`"Positive"` or `"Negative"`).
#'
#' @export
#'
#' @examples
#' rf_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                                package = "m6APrediction"))
#' example_df <- read.csv(
#'   system.file("extdata", "m6A_input_example.csv",
#'               package = "m6APrediction"),
#'   stringsAsFactors = FALSE
#' )
#' prediction_multiple(rf_fit, example_df)
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {
  required_cols <- c(
    "gc_content", "RNA_type", "RNA_region", "exon_length",
    "distance_to_junction", "evolutionary_conservation", "DNA_5mer"
  )
  stopifnot(all(required_cols %in% colnames(feature_df)))

  pred_df <- cbind(feature_df, dna_encoding(feature_df$DNA_5mer))
  pred_df$RNA_type <- factor(
    pred_df$RNA_type,
    levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")
  )
  pred_df$RNA_region <- factor(
    pred_df$RNA_region,
    levels = c("CDS", "intron", "3'UTR", "5'UTR")
  )

  pred_prob <- predict(ml_fit, newdata = pred_df, type = "prob")[, "Positive"]
  feature_df$predicted_m6A_prob <- pred_prob
  feature_df$predicted_m6A_status <- ifelse(
    pred_prob > positive_threshold, "Positive", "Negative"
  )
  return(feature_df)
}


#' Predict m6A status for a single site from individual features
#'
#' Convenience wrapper around the package function `prediction_multiple()`
#' for one site whose features are supplied as separate arguments rather
#' than a table.
#'
#' @inheritParams prediction_multiple
#' @param gc_content Numeric. GC content of the surrounding region in
#'   `[0, 1]`.
#' @param RNA_type Character. One of `"mRNA"`, `"lincRNA"`, `"lncRNA"`,
#'   `"pseudogene"`.
#' @param RNA_region Character. One of `"CDS"`, `"intron"`, `"3'UTR"`,
#'   `"5'UTR"`.
#' @param exon_length Numeric. Exon length feature.
#' @param distance_to_junction Numeric. Distance to the nearest splice
#'   junction.
#' @param evolutionary_conservation Numeric. Conservation score in
#'   `[0, 1]`.
#' @param DNA_5mer Character. The 5-mer sequence around the candidate
#'   adenosine (five bases, A/T/C/G only).
#'
#' @return A named character vector with two elements:
#'   `predicted_m6A_prob` (probability) and `predicted_m6A_status`
#'   (`"Positive"` or `"Negative"`).
#'
#' @export
#'
#' @examples
#' rf_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                                package = "m6APrediction"))
#' prediction_single(
#'   rf_fit,
#'   gc_content = 0.66,
#'   RNA_type = "lncRNA",
#'   RNA_region = "intron",
#'   exon_length = 0,
#'   distance_to_junction = 10.67,
#'   evolutionary_conservation = 0.016,
#'   DNA_5mer = "GGACC"
#' )
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction,
                              evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5) {
  single_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )
  result_df <- prediction_multiple(
    ml_fit, single_df, positive_threshold
  )
  returned_vector <- c(
    predicted_m6A_prob = result_df$predicted_m6A_prob,
    predicted_m6A_status = result_df$predicted_m6A_status
  )
  return(returned_vector)
}


#' Load the bundled random forest model
#'
#' Read the pre-trained random forest classifier stored in
#' `inst/extdata/rf_fit.rds`.
#'
#' @return A fitted `randomForest` classification object.
#' @export
#'
#' @examples
#' rf_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                                package = "m6APrediction"))
load_rf_model <- function() {
  readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
}
