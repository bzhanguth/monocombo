#' Evaluate Operating Characteristics for a Mono-Combo Seamless Clinical Trial Design
#'
#' Evaluates operating characteristics (OC) for a mono-combo seamless trial design across
#' multiple response rate scenarios and, optionally, multiple T2 arrival times. The function
#' wraps the simulation engine \code{mc_pro} and returns rejection probabilities, early stopping
#' probabilities, and average sample sizes for each scenario.
#'
#' @param n Numeric. Per-arm per-stage sample size (consistent with the output of
#'   \code{mc_Design}).
#' @param cv Numeric or numeric vector. Critical value(s) for the efficacy boundary. A single
#'   value is used for all scenarios (\code{sensitive = "method1"}); a vector of length equal to
#'   \code{t_arrive} applies the matching critical value for each arrival time
#'   (\code{sensitive = "method2"}).
#' @param p0 Numeric. Response rate in the control group.
#' @param p1 Numeric vector. Response rates for arm 1 (monotherapy) across scenarios.
#' @param p2 Numeric vector. Response rates for arm 2 (combination therapy) across scenarios.
#'   Must be the same length as \code{p1}.
#' @param t Numeric. Planned T2 addition time as a fraction of per-arm enrollment (e.g.,
#'   \code{t = 0.5} means T2 is added when half of arm 1 enrollment is complete). Default is
#'   \code{0.5}.
#' @param t_arrive Numeric or numeric vector. Actual T2 arrival time(s). If \code{NULL}
#'   (default), \code{t} is used. Supply a vector to evaluate OC across multiple arrival times.
#' @param platform Character. \code{"true"} for OR-Platform or Platform designs;
#'   \code{"false"} for MAMS designs. Default is \code{"true"}.
#' @param method Character. Alpha-spending method: \code{"OBF"} (O'Brien-Fleming) or
#'   \code{"POC"} (Pocock).
#' @param order Character. \code{"true"} for OR-Platform (order restriction \eqn{p_2 \ge p_1});
#'   \code{"false"} for Platform (no order assumption). Controls the number of output columns.
#' @param j Integer or \code{NULL}. Number of stages for MAMS designs
#'   (\code{platform = "false"}). Must be \code{2} or \code{3}. Ignored when
#'   \code{platform = "true"}. If \code{NULL} (default), \code{j} is inferred from
#'   \code{t}: a scalar \code{t} implies \code{j = 2}; a vector \code{t} of length
#'   \code{k} implies \code{j = k}.
#' @param sensitive Character or \code{NULL}. Sensitivity analysis method for critical value
#'   and boundary selection when \code{t_arrive} differs from \code{t}.
#'   \code{"method1"} uses a single fixed \code{cv} with boundaries based on the original
#'   planned \code{t}; \code{"method2"} uses the \code{j}-th element of \code{cv} with
#'   boundaries re-derived from the corresponding \code{t_arrive}. Default is \code{NULL}.
#' @param seed Integer. Random seed for reproducibility. Default is \code{100}.
#' @param nsim Integer. Number of simulation replicates per scenario.
#'
#' @return A matrix with one row per scenario (combination of \code{p1}/\code{p2} and
#'   \code{t_arrive}) and the following columns:
#' \describe{
#'   \item{p0}{Control response rate.}
#'   \item{p1}{Monotherapy response rate for this scenario.}
#'   \item{p2}{Combination therapy response rate for this scenario.}
#'   \item{t_arrive}{Actual T2 arrival time for this scenario.}
#'   \item{reject_1}{Probability of rejecting \eqn{H_0^1} (monotherapy efficacy).}
#'   \item{reject_2not1}{Probability of rejecting \eqn{H_0^2} (combination) but not \eqn{H_0^1}.}
#'   \item{reject_2}{Probability of rejecting \eqn{H_0^2} (combination efficacy).}
#'   \item{reject_at_least_one}{Disjunctive power; equals FWER under the global null.}
#'   \item{reject_1not2}{(Platform only) Probability of rejecting \eqn{H_0^1} but not \eqn{H_0^2}.}
#'   \item{PES_Eff}{Probability of early stopping for efficacy.}
#'   \item{PES_Fut}{Probability of early stopping for futility.}
#'   \item{PES}{Total probability of early stopping (\code{PES_Eff + PES_Fut}).}
#'   \item{ave_N0}{Average sample size for the control arm.}
#'   \item{ave_N1}{Average sample size for the monotherapy arm.}
#'   \item{ave_N2}{Average sample size for the combination arm.}
#'   \item{ave_TotalN}{Average total sample size across all arms.}
#' }
#'
#' @examples
#' # OR-Platform OC: FWER under global null and power under alternative
#' oc_or <- mc_OC(
#'   n        = 18,
#'   cv       = 2.32,
#'   p0       = 0.15,
#'   p1       = c(0.15, 0.4),
#'   p2       = c(0.15, 0.4),
#'   t        = 0.5,
#'   platform = "true",
#'   method   = "POC",
#'   order    = "true",
#'   seed     = 100,
#'   nsim     = 5000
#' )
#'
#' @export
#'
mc_OC = function(n, cv, p0, p1, p2, t = 0.5, t_arrive = NULL, platform = "true",
                 method, j = NULL, order = NULL, sensitive = NULL, seed, nsim) {

  # For MAMS designs, infer j from t when j is not user-specified:
  #   scalar t  -> j = 2  (2-stage design, t is the interim fraction)
  #   vector t  -> j = length(t)  (user provides full stage-fraction vector)
  if (platform == "false") {
    if (is.null(j)) {
      j <- if (length(t) > 1) length(t) else 2
    }
  }

  # MAMS designs (platform="false") always return 12 result values (including
  # reject1not2); platform designs vary by order flag.
  if (platform == "false") {
    ncol = 16
  } else if (order == "true") {
    ncol = 15
  } else {
    ncol = 16
  }

  # For MAMS (platform="false"), t is the stage-fraction vector — not an arrival
  # time grid. Loop only over p1/p2 scenarios (single pass, t_arrive = NULL).
  # For Platform (platform="true"), loop over t_arrive values if provided.
  if (platform == "false") {
    t_vals = NA          # single dummy pass
    total_rows = length(p1)
  } else {
    t_vals = if(is.null(t_arrive)) t else t_arrive
    total_rows = length(p1) * length(t_vals)
  }

 # total_rows = length(p1) * length(t_vals)
  results = matrix(nrow = total_rows, ncol = ncol)

  row_idx = 1

  for (i in 1:length(p1)) {

    for (ta_idx in 1:length(t_vals)) {

      # For MAMS: pass t_arrive=NULL so mc_pro uses its own stage fractions.
      # For Platform: use the current element of t_vals as the arrival time.
      if (platform == "false") {
        current_t_arrive = NULL
      } else {
        current_t_arrive = t_vals[ta_idx]
      }

      # Critical value selection (only meaningful for Platform sensitivity analysis)
      if (!is.null(sensitive) && sensitive == "method1") {
        current_cv = cv
      } else {
        current_cv = if(length(cv) >= ta_idx) cv[ta_idx] else cv[1]
      }

      # Call the simulation engine
      result = mc_pro(n = n, cv = current_cv, method = method, p0 = p0,
                      p1 = p1[i], p2 = p2[i], t = t, t_arrive = current_t_arrive,
                      platform = platform, j = j, sensitive = sensitive,
                      seed = seed, order = order, nsim = nsim)

      # Fill the matrix
      results[row_idx, 1] = p0
      results[row_idx, 2] = p1[i]
      results[row_idx, 3] = p2[i]
      results[row_idx, 4] = if(platform == "false") NA else current_t_arrive
      results[row_idx, 5:ncol] = result

      row_idx = row_idx + 1
    }
  }

  # Set Column Names
  base_cols = c("p0", "p1", "p2", "t_arrive")
  if (platform == "false" || order == "false") {
    colnames(results) = c(base_cols, "reject_1", "reject_2not1", "reject_2",
                          "reject_at_least_one", "reject_1not2", "PES_Eff",
                          "PES_Fut", "PES", "ave_N0", "ave_N1", "ave_N2", "ave_TotalN")
  } else {
    colnames(results) = c(base_cols, "reject_1", "reject_2not1", "reject_2",
                          "reject_at_least_one", "PES_Eff", "PES_Fut", "PES",
                          "ave_N0", "ave_N1", "ave_N2", "ave_TotalN")
  }

  return(results)
}
