#' Abundance table for example 2
#'
#' It contains the isotopic abundance for the analysis.
#' This example contains all the isotopes 196 through 204.
#' The isotopic composition of the spikes can be either measured along with the samples,
#' or you can use the certificates values provided by the company when you bought the enriched isotopic solution.
#' If you used isotopic dilution with double incubation to measure the transformations (e.g. IHg to MMHg and MMHg to IHg),
#' you need to add all four of the spikes to your abundance table.
#' Additionally, one isotope needs to remains unchanged, this means that the
#' isotope that you choose to be your natural need to have the natural abundance values.
#'
#' @format A data frame with 7 rows and 6 columns:
#' \describe{
#'   \item{isotope}{isotopic component}
#'   \item{ isotope_199, isotope_200, isotope_201, isotope_202, isotope_204}{isotopic composition of each spike and natural}
#' }
#' @source {Measured in-house to served as an example}
#' @examples
#' abundances_2
#'
"abundances_2"
