#' The Framingham Heart Study
#'
#' A subset of data from the Framingham Heart Study
#' Syed S Mahmood et al. “The Framingham Heart Study and the epidemiology of cardiovascular disease: a historical perspective”. 
#' In: The lancet 383.9921 (2014), pp. 999–1008.6
#' Data downloaded from SurvSet repository (https://github.com/ErikinBC/SurvSet), 
#' described https://arxiv.org/pdf/2203.03094.pdf
#'
#' @format ## `example_data`
#' A data frame with 4,699 rows and 10 columns:
#' \describe{
#'   \item{event}{Event occurred during follow-up}
#'   \item{time}{Time to event (or censoring)}
#'   \item{sbp, dbp}{Systolic and diastolic blood pressure}
#'   \item{scl}{scl}
#'   \item{age}{Age}
#'   \item{sex}{Sex}
#'   \item{month}{Month}
#'   \item{bmi_cat}{BMI in 3 categories: healthy, overweight and obese}
#' }
#' @source <https://github.com/ErikinBC/SurvSet/blob/main/SurvSet/_datagen/output/Framingham.csv>
"example_data"