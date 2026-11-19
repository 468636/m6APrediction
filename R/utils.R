#' Encode a 5-mer DNA sequence into per-position nucleotide factors
#'
#' Split each 5-mer string in a character vector into five single-base
#' columns (`nt_pos1` ... `nt_pos5`), converted to factors with the levels
#' A, T, C, G. This is the feature encoding expected by the random forest
#' model shipped with the package.
#'
#' @param dna_strings A character vector of 5-mer sequences (A/T/C/G).
#'
#' @return A data frame whose columns are `nt_pos1` ... `nt_pos5`, each a
#'   factor with levels `c("A", "T", "C", "G")`.
#'
#' @keywords internal
dna_encoding <- function(dna_strings) {
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(
    unlist(strsplit(dna_strings, "")),
    ncol = nn,
    byrow = TRUE
  )
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m, stringsAsFactors = FALSE)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}
