#' California proposition 99
#'
#' @docType data
#' @name CA_prop99
#'
#'
#' @format A data frame with 1209 rows and 4 variables:
#' \describe{
#'   \item{State}{US state name, character string}
#'   \item{Year}{Year, integer}
#'   \item{PacksPerCapita}{per-capita cigarette consumption, numeric}
#'   \item{treated}{the treatmed indicator 0: control, 1: treated, numeric}
#' }
#' @source Abadie, Alberto, Alexis Diamond, and Jens Hainmueller.
#'  "Synthetic control methods for comparative case studies: Estimating the effect of California’s tobacco control program."
#'   Journal of the American statistical Association 105, no. 490 (2010): 493-505.
#'
#' @usage data(CA_prop99)
#'
NULL


#' CPS
#'
#' @docType data
#' @name CPS
#'
#' @format A data frame with 2000 rows and 8 variables.
#' \describe{
#'   \item{state}{state}
#'   \item{year}{year}
#'   \item{log_wage}{log_wage}
#'   \item{hours}{hours}
#'   \item{urate}{urate}
#'   \item{min_wage}{min_wage}
#'   \item{open_carry}{open_carry}
#'   \item{abort_ban}{abort_ban}
#' }
#'
#' @usage data(CPS)
#'
NULL

#' CHC
#'
#' @docType data
#' @name CHC
#'
#' @format A data frame with 91770 rows and 6 variables.
#' \describe{
#'   \item{county}{County id}
#'   \item{mort}{Mortality: deaths per 100,000 residents}
#'   \item{year}{Year}
#'   \item{CHC}{Treatment indicator}
#'   \item{urb_share}{Urban share: urban population percentage in 1960}
#'   \item{popwt}{Population of the county in 1960}
#' }
#' @source Bailey, Martha J., and Andrew Goodman-Bacon.
#'   "The War on Poverty's experiment in public medicine: Community health centers and the mortality of older Americans."
#'    American Economic Review 105.3 (2015): 1067-1104.
#'
#' @usage data(CHC)
#'
NULL
