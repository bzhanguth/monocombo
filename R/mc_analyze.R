#' Stage-Wise Decision Analysis for a Monotherapy-Combination Therapy Trial
#'
#' Applies the pre-specified hierarchical decision rules at each interim and final analysis
#' of a mono-combo seamless trial involving a monotherapy arm (T1) and a combination therapy
#' arm (T2). Supply only the data available at the current stage; the function determines
#' which stage has been reached from the columns present in \code{data}.
#'
#' @param data A list or data frame containing observed response counts and the total per-arm
#'   sample size. Include only the columns corresponding to stages that have been completed:
#'   \itemize{
#'     \item \code{n}: \strong{Total} per-arm sample size across all planned stages (consistent
#'       with the \code{n} output of \code{mc_Design}). Stage-specific sizes are
#'       derived internally as \code{round(n * t_arrive)} for Stage 1,
#'       \code{round(n * (1 - t_arrive))} for Stage 2, and \code{round(n * t_arrive)} for
#'       Stage 3.
#'     \item \code{y0s1}, \code{y1s1}: Observed responses for control and monotherapy in Stage 1.
#'     \item \code{y0s2}, \code{y1s2}, \code{y2s2}: Observed responses for control, monotherapy,
#'       and combination in Stage 2.
#'     \item \code{y0s3}, \code{y2s3}: Observed responses for control and combination in Stage 3.
#'   }
#' @param u Numeric. Common final-analysis efficacy boundary for T1 and T2;
#'   corresponds to \code{u21 = u32} in the \code{mc_Design} output. First-look
#'   boundaries are derived internally from \code{u}, \code{method}, \code{t}, and
#'   \code{t_arrive}.
#' @param t Numeric. Planned T2 arrival time as a fraction of per-arm enrollment.
#'   Default is \code{0.5}.
#' @param t_arrive Numeric. Realized T2 arrival time. Defaults to \code{t}.
#' @param N_max Character or \code{NULL}. Set to \code{"true"} when analyzing a
#'   delayed-arrival design whose total sample size was held fixed. This preserves
#'   the planned \code{t} when deriving the T2 first-look boundary and Stage-3
#'   sample size. Otherwise, a delayed design is treated as a complete redesign
#'   based on \code{t_arrive}. Default is \code{NULL}.
#' @param method Character. Boundary family: \code{"OBF"} or \code{"POC"}.
#' @param platform Character. Must be \code{"true"}; the current decision tree is
#'   for a rolling platform design.
#' @param order Character. Must be \code{"true"}; the current decision tree uses
#'   the order-restricted Stage-3 combined statistic.
#'
#' @return A named character/numeric vector. Numeric elements are the z-statistics computed
#'   at each look; character elements are the stage decisions:
#'   \itemize{
#'     \item \code{z11}: Z-statistic for T1 vs. control at Stage 1.
#'     \item \code{z21}: Cumulative Z-statistic for T1 vs. control at Stage 2.
#'     \item \code{z22}: Z-statistic for T2 vs. concurrent control at Stage 2.
#'     \item \code{z32}: Cumulative Z-statistic for T2 vs. control at Stage 3.
#'     \item \code{z312}: Pooled isotonic Z-statistic at Stage 3 (OR-Platform only; computed
#'       when \code{z32 < u} to leverage the order assumption \eqn{p_2 \ge p_1}).
#'     \item \code{Stage_1_Decision}: One of "Early Stop Trial - Select T1",
#'       "Proceed T1 to Stage 2 and Add T2", or "Drop T1, Add T2".
#'     \item \code{Stage_2_Decision}: One of "Early Stop Trial - Select T1 and T2",
#'       "Early Stop Trial - T2", "Proceed T2 to Stage 3", or
#'       "Early Terminate Trial due to Futile".
#'     \item \code{Stage_3_Decision}: One of "Select T1 and T2", "Select T2", or "Both Futile".
#'   }
#'   Only z-statistics and decisions relevant to the stages present in \code{data} are returned.
#'
#' @details
#' The lower futility boundary is fixed at \code{l1 = 0} for all stages. The decision
#' tree at each stage is:
#' \enumerate{
#'   \item \strong{Stage 1} (requires \code{y0s1}, \code{y1s1}):
#'     \itemize{
#'       \item \eqn{z_{11} \ge u_{11}}: early stop, select T1.
#'       \item \eqn{0 \le z_{11} < u_{11}}: proceed T1 to Stage 2 and add T2.
#'       \item \eqn{z_{11} < 0}: drop T1, add T2 only.
#'     }
#'   \item \strong{Stage 2} (requires additionally \code{y0s2}, \code{y1s2}, \code{y2s2}):
#'     \itemize{
#'       \item Cumulative T1 (\eqn{z_{21}}) and concurrent T2 (\eqn{z_{22}}) are tested.
#'       \item Various combinations of \eqn{z_{21}} and \eqn{z_{22}} relative to
#'         \code{u11}, \code{u22}, and \code{u} determine early stopping,
#'         continuation, or futility.
#'     }
#'   \item \strong{Stage 3} (requires additionally \code{y0s3}, \code{y2s3}):
#'     \itemize{
#'       \item Cumulative T2 (\eqn{z_{32}}) is tested against \code{u}.
#'       \item If \eqn{z_{32} < u_2}, the pooled isotonic statistic \eqn{z_{312}} is
#'         computed by pooling T1 and T2 data under the order assumption.
#'       \item The Stage-2 control cohort enters the pooled control mean twice and
#'         therefore contributes a factor of \eqn{2^2=4} to its variance.
#'     }
#' }
#'
#' @examples
#' # Stage 1 only: proceed to Stage 2 and add combination arm
#' data_s1 <- list(n = 36, y0s1 = 4, y1s1 = 8)
#' mc_Analyze(data_s1, u = 2.32, t = 0.5)
#'
#' # Stage 2: combination arm borderline, proceed to Stage 3
#' data_s2 <- list(y0s1 = 3, y0s2 = 4, y1s1 = 6, y1s2 = 5, y2s2 = 9, n = 36)
#' mc_Analyze(data_s2, u = 2.32, t = 0.5)
#'
#' # Final analysis (Stage 3): isotonic test selects both arms
#' data_s3 <- list(y0s1 = 3, y0s2 = 4, y0s3 = 5,
#'                 y1s1 = 7, y1s2 = 9,
#'                 y2s2 = 7, y2s3 = 7, n = 36)
#' mc_Analyze(data_s3, u = 2.32, t = 0.5)
#'
#' @export
mc_Analyze <- function(data, u, t = 0.5, t_arrive = t, N_max = NULL,
                       method = "OBF", platform = "true", order = "true") {
  if (!identical(platform, "true")) {
    stop("mc_Analyze() currently supports only rolling platform designs (platform = \"true\").")
  }
  if (!identical(order, "true")) {
    stop("mc_Analyze() currently supports only the order-restricted decision tree (order = \"true\").")
  }
  if (!method %in% c("OBF", "POC")) {
    stop("method must be either \"OBF\" or \"POC\".")
  }
  if (length(u) != 1L || !is.finite(u)) stop("u must be one finite numeric value.")
  if (length(t) != 1L || !is.finite(t) || t <= 0 || t >= 1) {
    stop("t must be one number strictly between 0 and 1.")
  }
  if (length(t_arrive) != 1L || !is.finite(t_arrive) || t_arrive <= 0 || t_arrive >= 1) {
    stop("t_arrive must be one number strictly between 0 and 1.")
  }

  fixed_N <- identical(N_max, "true")
  stage3_fraction <- if (fixed_N) 0.5 * (t + t_arrive) else t_arrive
  if (method == "POC") {
    u11 <- u22 <- u
  } else {
    u11 <- u / sqrt(t_arrive)
    t_for_t2 <- if (fixed_N) t else t_arrive
    u22 <- u * sqrt(((1 - t_arrive) + 0.5 * (t_for_t2 + t_arrive)) /
                    (1 - t_arrive))
  }
  result <- c()  # Initialize an empty vector to store results
  present_columns <- names(data)  # Get the column names from the input data
  n <- data$n  # Get the base sample size
  l1 <- 0  # Lower bound is fixed at 0

  # ---- Stage-specific sample sizes (n = total per-arm, consistent with mc_pro) ----
  n_s1 <- round(n * t_arrive)          # Stage 1 per-arm size
  n_s2 <- round(n * (1 - t_arrive))    # Stage 2 per-arm size
  n_s3 <- round(n * stage3_fraction)   # Stage 3 per-arm size

  # ---- Stage 1 Decision ----
  if ("y0s1" %in% present_columns && "y1s1" %in% present_columns) {
    y10 <- data$y0s1
    y11 <- data$y1s1
    t11 <- round(proportion_test(y11, y10, n_s1, n_s1), 3)  # Proportion test for Stage 1
    result["z11"] <- t11  # Save z-test result for Stage 1

    if (t11 >= u11) {
      result["Stage_1_Decision"] <- "Early Stop Trial - Select T1"

    }

    else if (l1 <= t11 && t11 < u11) {

      # ---- Stage 2 Decision ----
      if ("y0s2" %in% present_columns) {
        y20 <- data$y0s1 + data$y0s2  # Aggregate control group data for Stages 1 and 2
        y21 <- data$y1s1 + data$y1s2
        t21 <- round(proportion_test(y21, y20, n_s1 + n_s2, n_s1 + n_s2), 3)
        result["z21"] <- t21  # Save z-test result for T1 in Stage 2

        y22 <- data$y2s2
        t22 <- round(proportion_test(y22, data$y0s2, n_s2, n_s2), 3)
        result["z22"] <- t22  # Save z-test result for T2 in Stage 2

        if (t21 >= u){

          if (t22 >= u22) {
            result["Stage_2_Decision"] <- "Early Stop Trial - Select T1 and T2"  }
          else if (t22 < u22) {
            # ---- Stage 3 Decision ----
            if ("y2s3" %in% present_columns && "y0s3" %in% present_columns) {
              y30 <- data$y0s2 + data$y0s3  # Aggregate control group data for Stages 2 and 3
              y32 <- data$y2s2 + data$y2s3  # Aggregate T2 data for Stages 2 and 3
              t32 <- round(proportion_test(y32, y30, n_s2 + n_s3, n_s2 + n_s3), 3)  # Proportion test for Stage 3
              result["z32"] <- t32  # Save z-test result for Stage 3

              if (t32 >= u) {
                result["Stage_3_Decision"] <- "Select T1 and T2" }
              else {
                y30 <- data$y0s1 + 2 * data$y0s2 + data$y0s3
                y312 <- data$y1s1 + data$y1s2 + data$y2s2 + data$y2s3
                n_combined <- n_s1 + 2 * n_s2 + n_s3

                # The Stage-2 control count has coefficient 2 in the weighted
                # control mean, so its variance contribution has coefficient 4.
                p_t <- y312 / n_combined
                p_c <- y30 / n_combined
                var_t <- p_t * (1 - p_t) / n_combined
                var_c <- (n_s1 + 4 * n_s2 + n_s3) *
                  p_c * (1 - p_c) / n_combined^2
                se312 <- sqrt(var_t + var_c)
                if (!is.finite(se312) || se312 <= 0) {
                  n_combined_adjusted <- n_combined + 2
                  p_t <- (y312 + 1) / n_combined_adjusted
                  p_c <- (y30 + 1) / n_combined_adjusted
                  var_t <- p_t * (1 - p_t) / n_combined_adjusted
                  var_c <- (n_s1 + 4 * n_s2 + n_s3 + 2) *
                    p_c * (1 - p_c) / n_combined_adjusted^2
                  se312 <- sqrt(var_t + var_c)
                }
                t312 <- if (is.finite(se312) && se312 > 0) {
                  round((p_t - p_c) / se312, 3)
                } else {
                  0
                }
                result["z312"] <- t312
                if (t312 >= u) {result["Stage_3_Decision"] <- "Select T1 and T2" }
                else{result["Stage_3_Decision"] <- "Both Futile"}
              }
            }
            else {result["Stage_2_Decision"] <- "Proceed T2 to Stage 3"}
          }
        }

        if (t21 < u){
          if (t22 >= u22) {
            result["Stage_2_Decision"] <- "Early Stop Trial - T2" }
          else if (t22 < u22 && t22>=l1) {
            # ---- Stage 3 Decision ----
            if ("y2s3" %in% present_columns && "y0s3" %in% present_columns) {
              y30 <- data$y0s2 + data$y0s3  # Aggregate control group data for Stages 2 and 3
              y32 <- data$y2s2 + data$y2s3  # Aggregate T2 data for Stages 2 and 3
              t32 <- round(proportion_test(y32, y30, n_s2 + n_s3, n_s2 + n_s3), 3)  # Proportion test for Stage 3
              result["z32"] <- t32  # Save z-test result for Stage 3

              if (t32 >= u) {
                result["Stage_3_Decision"] <- "Select T2" }
              else { result["Stage_3_Decision"] <- "Both Futile" }

            }

            else {result["Stage_2_Decision"] <- "Proceed T2 to Stage 3"}
          }

          else {result["Stage_2_Decision"] <- "Early Terminate Trial due to Futile"}
        }

      }

      else {result["Stage_1_Decision"] <- "Proceed T1 to Stage 2 and Add T2"}

    }

    else {
      if ("y0s2" %in% present_columns) {
        y22 <- data$y2s2
        t22 <- round(proportion_test(y22, data$y0s2, n_s2, n_s2), 3)
        result["z22"] <- t22  # Save z-test result for T2 in Stage 2

        if (t22 >= u22) {
          result["Stage_2_Decision"] <- "Early Stop Trial - T2" }
        else if (t22 < u22 && t22>=l1) {
          # ---- Stage 3 Decision ----
          if ("y2s3" %in% present_columns && "y0s3" %in% present_columns) {
            y30 <- data$y0s2 + data$y0s3  # Aggregate control group data for Stages 2 and 3
            y32 <- data$y2s2 + data$y2s3  # Aggregate T2 data for Stages 2 and 3
            t32 <- round(proportion_test(y32, y30, n_s2 + n_s3, n_s2 + n_s3), 3)  # Proportion test for Stage 3
            result["z32"] <- t32  # Save z-test result for Stage 3

            if (t32 >= u) {
              result["Stage_3_Decision"] <- "Select T2" }
            else { result["Stage_3_Decision"] <- "Both Futile" }

          }

          else {result["Stage_2_Decision"] <- "Proceed T2 to Stage 3"}
        }

        else {result["Stage_2_Decision"] <- "Early Terminate Trial due to Futile"}
      }
      else {
        result["Stage_1_Decision"] <- "Drop T1, Add T2" }

    }
  }

  return(result)
}
