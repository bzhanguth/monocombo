#' monocombo: Mono-Combo Seamless Clinical Trial Design, Simulation, and Analysis
#'
#' The \pkg{monocombo} package provides a unified framework for designing, evaluating
#' operating characteristics of, and analyzing multi-arm multi-stage (MAMS) seamless
#' clinical trials that start with monotherapy and roll in combination therapy mid-trial.
#' It supports four design flavors: OR-Platform, Platform, 2-stage MAMS, and 3-stage MAMS.
#'
#' @section Core user-facing functions:
#' \describe{
#'   \item{\code{\link{mc_Design}}}{Compute critical values and sample size.}
#'   \item{\code{\link{mc_OC}}}{Evaluate operating characteristics via simulation.}
#'   \item{\code{\link{mc_Analyze}}}{Stage-wise decision analysis on observed data.}
#' }
#'
#' @section Internal simulation engine:
#' \describe{
#'   \item{\code{\link{mc_pro}}}{Single-scenario simulation engine (called by \code{mc_OC}).}
#'   \item{\code{\link{proportion_test}}}{Unpooled two-sample proportion z-test.}
#' }
#'
#' @importFrom mvtnorm pmvnorm
#' @importFrom MASS fractions
#'
#' @docType package
#' @name monocombo
"_PACKAGE"
