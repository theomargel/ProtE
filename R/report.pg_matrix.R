#' DIA-NN Protein Quantification Matrix
#'
#' This dataset contains protein quantification results from DIA-NN or fragPipe analysis.
#' The matrix "report.pg_matrix" includes normalized protein intensities for each sample along with protein annotations like IDs and protein names in the format of protein groups.
#' This output is typically used for downstream analysis such as differential expression analysis, pathway analysis, and clustering.
#'
#' @format A data frame (or matrix) with rows representing proteins / rotein groups and columns representing samples or experimental conditions.
#' The matrix includes both protein annotations and quantification data:
#' \describe{
#'   \item{Protein.Group}{Character. The unique identifier for each protein group (e.g., UniProt ID).}
#'   \item{Protein.Ids}{Character. The same as Protein.Group.}
#'   \item{Sample 1}{Numeric. The quantitative value (e.g., intensity or area) for each protein in Sample 1.}
#'   \item{Sample 2}{Numeric  The quantitative value for each protein in Sample 2 }
#'   \item{Sample 3}{Numeric  The quantitative value for each protein in Sample 3 }
#'   \item{Sample 4}{Numeric  The quantitative value for each protein in Sample 4 }
#'   \item{Sample 5}{Numeric  The quantitative value for each protein in Sample 5 }
#'   \item{Sample 6}{Numeric  The quantitative value for each protein in Sample 6 }
#'   \item{Sample 7}{Numeric  The quantitative value for each protein in Sample 7 }
#'   \item{Sample 8}{Numeric  The quantitative value for each protein in Sample 8 }
#'   \item{Sample 9}{Numeric  The quantitative value for each protein in Sample 9 }
#'   \item{Sample 10}{Numeric  The quantitative value for each protein in Sample 10 }
#'   \item{Sample 11}{Numeric  The quantitative value for each protein in Sample 11 }
#'   \item{Sample 12}{Numeric  The quantitative value for each protein in Sample 12 }
#'   \item{Sample 13}{Numeric  The quantitative value for each protein in Sample 13 }
#'   \item{Sample 14}{Numeric  The quantitative value for each protein in Sample 14 }
#'   \item{Sample 15}{Numeric  The quantitative value for each protein in Sample 15 }
#'   \item{Sample 16}{Numeric  The quantitative value for each protein in Sample 16 }
#'   \item{Sample 17}{Numeric  The quantitative value for each protein in Sample 17 }
#'   \item{Sample 18}{Numeric  The quantitative value for each protein in Sample 18 }
#'   \item{Sample 19}{Numeric  The quantitative value for each protein in Sample 19 }
#'   \item{Sample 20}{Numeric  The quantitative value for each protein in Sample 20 }
#'   \item{Sample 21}{Numeric  The quantitative value for each protein in Sample 21 }
#'   \item{Sample 22}{Numeric  The quantitative value for each protein in Sample 22 }
#'   \item{Sample 23}{Numeric  The quantitative value for each protein in Sample 23 }
#'  \item{Protein.Names}{Character. The name of the proteins}
#'   \item{Genes}{Character. The gene symbol of the protein}
#'   \item{First.Protein.Description}{Empty in the output as of 2024}
#'
#'
#' }
#' @source DIA-NN analysis results, generated from quantitative proteomics experiments.
#' @name report.pg_matrix
#' @docType data
NULL
