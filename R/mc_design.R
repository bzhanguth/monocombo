
# Unified functions supporting:
#   OR-Platform  (platform="true",  order="true")
#   Platform     (platform="true",  order="false")
#   2-stage MAMS (platform="false", j=2)
#   3-stage MAMS (platform="false", j=3)
#
# User-facing function: mc_Design()
# Internal helpers:     get_corr_unified(), mc_pro_analytical(), mc_analytic_OC()

# --- 1. Unified Correlation Matrix ---
# platform="true"  : t is a scalar (planned addition time), t_arrive optional
# platform="false" : t is a cumulative-fraction vector, j = length(t)
get_corr_unified <- function(n, t, p0, p1, p2, platform = "false", t_arrive = NULL, j = NULL) {
  v0 <- p0 * (1 - p0); v1 <- p1 * (1 - p1); v2 <- p2 * (1 - p2)

  if (platform == "true") {
    # Existing Platform Logic
    nc1 <- if(is.null(t_arrive)) n*t else n*t_arrive
    nc2 <- if(is.null(t_arrive)) n*(1-t) else n*(1-t_arrive)
    nc3 <- if(is.null(t_arrive)) n*t else n*0.5*(t+t_arrive)

    lambda1 <- lambda2 <- 1
    v10=v20=v30=v0; v11=v21=v1; v22=v32=v2
    n11=lambda1*nc1; n12=lambda1*nc2; n22=lambda2*nc2; n23=lambda2*nc3
    nE1=n11+n12; nE2=n22+n23; nE=nE1+nE2
    w1_12=nc1/(nc1+nc2); w2_12=nc2/(nc1+nc2); w2_23=nc2/(nc2+nc3); w3_23=nc3/(nc2+nc3)

    var11=v11/n11 + v10/nc1
    var21=(nc1^2*(v11/n11+v10/nc1)+nc2^2*(v21/n12+v20/nc2))/(nc1+nc2)^2
    var22=v22/n22+v20/nc2
    var32=(nc2^2*(v22/n22+v20/nc2)+nc3^2*(v32/n23+v30/nc3))/(nc2+nc3)^2

    var_pE1_treat=(nc1^2*v11/n11+nc2^2*v21/n12)/(nc1+nc2)^2
    var_pE2_treat=(nc2^2*v22/n22+nc3^2*v32/n23)/(nc2+nc3)^2
    var_pE=(nE1^2*var_pE1_treat+nE2^2*var_pE2_treat)/nE^2
    var_pCstar=(lambda1^2*nc1*v10+(lambda1+lambda2)^2*nc2*v20+lambda2^2*nc3*v30)/nE^2
    var_ct=var_pE+var_pCstar

    cov11_21=w1_12*var11; cov21_22=w2_12*v20/nc2; cov22_32=w2_23*var22; cov21_32=w2_12*w2_23*v20/nc2
    cov11_ct=(v11+lambda1*v10)/nE; cov22_ct=(v22+2*lambda2*v20)/nE
    cov21_ct=(nc1*v11+lambda1*nc1*v10+nc2*v21+(lambda1+lambda2)*nc2*v20)/(nE*(nc1+nc2))
    cov32_ct=(nc3*v32+lambda2*nc3*v30+nc2*v22+(lambda1+lambda2)*nc2*v20)/(nE*(nc3+nc2))

    Sigma <- diag(1, 5)
    Sigma[1,2] <- Sigma[2,1] <- cov11_21 / sqrt(var11*var21)
    Sigma[2,3] <- Sigma[3,2] <- cov21_22 / sqrt(var21*var22)
    Sigma[3,4] <- Sigma[4,3] <- cov22_32 / sqrt(var22*var32)
    Sigma[2,4] <- Sigma[4,2] <- cov21_32 / sqrt(var21*var32)
    Sigma[1,5] <- Sigma[5,1] <- cov11_ct / sqrt(var11*var_ct)
    Sigma[2,5] <- Sigma[5,2] <- cov21_ct / sqrt(var21*var_ct)
    Sigma[3,5] <- Sigma[5,3] <- cov22_ct / sqrt(var22*var_ct)
    Sigma[4,5] <- Sigma[5,4] <- cov32_ct / sqrt(var32*var_ct)
    return(list(Sigma=Sigma, variances=c(var11, var21, var22, var32, var_ct)))

  } else {
    # Flexible MAMS Correlation Logic
    dim_size <- j * 2
    Sigma <- matrix(0, nrow = dim_size, ncol = dim_size); diag(Sigma) <- 1
    vars <- numeric(dim_size)
    for (k in 1:j) {
      vars[(k*2-1)] <- v1/(n*t[k]) + v0/(n*t[k])
      vars[(k*2)]   <- v2/(n*t[k]) + v0/(n*t[k])
    }
    for (r in 1:dim_size) {
      for (c in 1:dim_size) {
        if (r == c) next
        s_r <- ceiling(r/2); tr_r <- ifelse(r %% 2 == 1, 1, 2)
        s_c <- ceiling(c/2); tr_c <- ifelse(c %% 2 == 1, 1, 2)
        t_max <- t[max(s_r, s_c)]
        v_tr_shared <- if(tr_r == tr_c) (if(tr_r == 1) v1 else v2) else 0
        cov_val <- (v_tr_shared + v0) / (n * t_max)
        Sigma[r, c] <- cov_val / (sqrt(vars[r]) * sqrt(vars[c]))
      }
    }
    return(list(Sigma = Sigma, variances = vars))
  }
}

# --- 2. Unified Probability Calculation ---
mc_pro_analytical <- function(n, method, p0, p1, p2, cv,
                               j = NULL, t = NULL, t_arrive = NULL,
                               platform = "false", order = "true",
                               sensitive = NULL) {

  if (platform == "false") {
    # ---- MAMS Logic ----
    if (is.null(j)) j <- 3
    if (is.null(t)) t <- seq(1/j, 1, length.out = j)
    u <- if (method == "OBF") cv / sqrt(t) else rep(cv, j)
    l <- c(rep(0, j - 1), u[j])

    res <- get_corr_unified(n, t, p0, p1, p2, platform = "false", j = j)
    mu  <- rep(c(p1 - p0, p2 - p0), j) / sqrt(res$variances)

    if (j == 3) {
      regions <- list(
        P1  = list(L = c(u[1], u[1], -Inf, -Inf, -Inf, -Inf), U = c(Inf,  Inf,  Inf,  Inf,  Inf,  Inf)),
        P2  = list(L = c(u[1], -Inf, -Inf, -Inf, -Inf, -Inf), U = c(Inf,  u[1], Inf,  Inf,  Inf,  Inf)),
        P3  = list(L = c(l[1], u[1], -Inf, -Inf, -Inf, -Inf), U = c(u[1], Inf,  Inf,  Inf,  Inf,  Inf)),
        P4  = list(L = c(l[1], l[1], u[2], u[2], -Inf, -Inf), U = c(u[1], u[1], Inf,  Inf,  Inf,  Inf)),
        P5  = list(L = c(l[1], l[1], u[2], -Inf, -Inf, -Inf), U = c(u[1], u[1], Inf,  u[2], Inf,  Inf)),
        P6  = list(L = c(l[1], l[1], l[2], u[2], -Inf, -Inf), U = c(u[1], u[1], u[2], Inf,  Inf,  Inf)),
        P7  = list(L = c(l[1], l[1], l[2], l[2], u[3], u[3]), U = c(u[1], u[1], u[2], u[2], Inf,  Inf)),
        P8  = list(L = c(l[1], l[1], l[2], l[2], u[3], -Inf), U = c(u[1], u[1], u[2], u[2], Inf,  u[3])),
        P9  = list(L = c(l[1], l[1], l[2], l[2], -Inf, u[3]), U = c(u[1], u[1], u[2], u[2], u[3], Inf)),
        P10 = list(L = c(l[1], l[1], l[2], l[2], -Inf, -Inf), U = c(u[1], u[1], u[2], u[2], u[3], u[3])),
        P11 = list(L = c(l[1], l[1], l[2], -Inf, u[3], -Inf), U = c(u[1], u[1], u[2], l[2], Inf,  Inf)),
        P12 = list(L = c(l[1], l[1], l[2], -Inf, -Inf, -Inf), U = c(u[1], u[1], u[2], l[2], u[3], Inf)),
        P13 = list(L = c(l[1], l[1], -Inf, u[2], -Inf, -Inf), U = c(u[1], u[1], l[2], Inf,  Inf,  Inf)),
        P14 = list(L = c(l[1], l[1], -Inf, l[2], -Inf, u[3]), U = c(u[1], u[1], l[2], u[2], Inf,  Inf)),
        P15 = list(L = c(l[1], l[1], -Inf, l[2], -Inf, -Inf), U = c(u[1], u[1], l[2], u[2], Inf,  u[3])),
        P16 = list(L = c(l[1], l[1], -Inf, -Inf, -Inf, -Inf), U = c(u[1], u[1], l[2], l[2], Inf,  Inf)),
        P17 = list(L = c(l[1], -Inf, u[2], -Inf, -Inf, -Inf), U = c(u[1], l[1], Inf,  Inf,  Inf,  Inf)),
        P18 = list(L = c(l[1], -Inf, l[2], -Inf, u[3], -Inf), U = c(u[1], l[1], u[2], Inf,  Inf,  Inf)),
        P19 = list(L = c(l[1], -Inf, l[2], -Inf, -Inf, -Inf), U = c(u[1], l[1], u[2], Inf,  u[3], Inf)),
        P20 = list(L = c(l[1], -Inf, -Inf, -Inf, -Inf, -Inf), U = c(u[1], l[1], l[2], Inf,  Inf,  Inf)),
        P21 = list(L = c(-Inf, u[1], -Inf, -Inf, -Inf, -Inf), U = c(l[1], Inf,  Inf,  Inf,  Inf,  Inf)),
        P22 = list(L = c(-Inf, l[1], -Inf, u[2], -Inf, -Inf), U = c(l[1], u[1], Inf,  Inf,  Inf,  Inf)),
        P23 = list(L = c(-Inf, l[1], -Inf, l[2], -Inf, u[3]), U = c(l[1], u[1], Inf,  u[2], Inf,  Inf)),
        P24 = list(L = c(-Inf, l[1], -Inf, l[2], -Inf, -Inf), U = c(l[1], u[1], Inf,  u[2], Inf,  u[3])),
        P25 = list(L = c(-Inf, l[1], -Inf, -Inf, -Inf, -Inf), U = c(l[1], u[1], Inf,  l[2], Inf,  Inf)),
        P26 = list(L = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), U = c(l[1], l[1], Inf,  Inf,  Inf,  Inf))
      )
      probs <- sapply(regions, function(reg) {
        pmvnorm(lower = reg$L, upper = reg$U, mean = mu, corr = res$Sigma)[1]
      })
      rej_all  <- sum(probs[c("P1","P4","P7")])
      rej1not2 <- sum(probs[c("P2","P5","P8","P11","P17","P18")])
      rej2not1 <- sum(probs[c("P3","P6","P9","P13","P14","P21","P22","P23")])
      result <- c(reject1 = rej_all + rej1not2,
                  reject2 = rej_all + rej2not1,
                  reject2not1 = rej2not1,
                  reject_at_least_one = rej_all + rej1not2 + rej2not1)

    } else {
      # j = 2 MAMS
      regions <- list(
        P1  = list(L = c(u[1], u[1], -Inf, -Inf), U = c(Inf,  Inf,  Inf, Inf)),
        P2  = list(L = c(u[1], l[1], -Inf, -Inf), U = c(Inf,  u[1], Inf, Inf)),
        P3  = list(L = c(l[1], u[1], -Inf, -Inf), U = c(u[1], Inf,  Inf, Inf)),
        P4  = list(L = c(l[1], l[1], u[2], u[2]), U = c(u[1], u[1], Inf, Inf)),
        P5  = list(L = c(l[1], l[1], u[2], -Inf), U = c(u[1], u[1], Inf, u[2])),
        P6  = list(L = c(l[1], l[1], -Inf, u[2]), U = c(u[1], u[1], u[2], Inf)),
        P8  = list(L = c(l[1], -Inf, u[2], -Inf), U = c(u[1], l[1], Inf, Inf)),
        P10 = list(L = c(-Inf, u[1], -Inf, -Inf), U = c(l[1], Inf,  Inf, Inf)),
        P11 = list(L = c(-Inf, l[1], -Inf, u[2]), U = c(l[1], u[1], Inf, Inf))
      )
      probs <- sapply(regions, function(reg) {
        pmvnorm(lower = reg$L, upper = reg$U, mean = mu, corr = res$Sigma)[1]
      })
      rej_all  <- sum(probs[c("P1","P4")])
      rej2not1 <- sum(probs[c("P6","P3","P10","P11")])
      rej1not2 <- sum(probs[c("P2","P8","P5")])
      result <- c(reject1 = rej_all + rej1not2,
                  reject2 = rej_all + rej2not1,
                  reject2not1 = rej2not1,
                  reject_at_least_one = rej_all + rej2not1 + rej1not2)
    }

  } else {
    # ---- Platform Logic ----
    if (is.null(t)) t <- 0.5
    l_plat <- 0

    # Define stage-specific boundaries before building region lists
    if (method == "POC") {
      u11 <- u21 <- u22 <- u32 <- cv
    } else {  # OBF
      # sensitive="method1" forces use of original planned t even when t_arrive differs
      use_t_only <- identical(sensitive, "method1") || is.null(t_arrive)
      if (use_t_only) {
        u11 <- cv
        u21 <- cv * sqrt(t)
        u32 <- cv * sqrt(t)
        u22 <- cv * sqrt(t) / sqrt(1 - t)
      } else {
        u11 <- cv
        u21 <- cv * sqrt(t_arrive)
        u32 <- cv * sqrt(t_arrive)
        u22 <- cv * sqrt(t_arrive) * sqrt(((1 - t_arrive) + 0.5 * (t + t_arrive)) / (1 - t_arrive))
      }
    }

    res     <- get_corr_unified(n, t, p0, p1, p2, platform = "true", t_arrive = t_arrive)
    mu      <- c(p1-p0, p1-p0, p2-p0, p2-p0, ((p1+p2)/2)-p0) / sqrt(res$variances)
    dim_use <- ifelse(order == "true", 5, 4)

    if (order == "false") {
      regions <- list(
        po1  = list(lower = c(u11,    -Inf,    -Inf,    -Inf),   upper = c(Inf,    Inf,  Inf,  Inf)),
        po2  = list(lower = c(l_plat, u21,     u22,     -Inf),   upper = c(u11,    Inf,  Inf,  Inf)),
        po3  = list(lower = c(l_plat, u21,     l_plat,  u32),    upper = c(u11,    Inf,  u22,  Inf)),
        po4  = list(lower = c(l_plat, u21,     l_plat,  -Inf),   upper = c(u11,    Inf,  u22,  u32)),
        po5  = list(lower = c(l_plat, u21,     -Inf,    -Inf),   upper = c(u11,    Inf,  l_plat, Inf)),
        po6  = list(lower = c(l_plat, -Inf,    u22,     -Inf),   upper = c(u11,    u21,  Inf,  Inf)),
        po7  = list(lower = c(l_plat, -Inf,    l_plat,  u32),    upper = c(u11,    u21,  u22,  Inf)),
        po10 = list(lower = c(-Inf,   -Inf,    u22,     -Inf),   upper = c(l_plat, Inf,  Inf,  Inf)),
        po11 = list(lower = c(-Inf,   -Inf,    l_plat,  u32),    upper = c(l_plat, Inf,  u22,  Inf))
      )
    } else {
      regions <- list(
        po1  = list(lower = c(u11,    -Inf,    -Inf,    -Inf,    -Inf), upper = c(Inf,    Inf, Inf, Inf, Inf)),
        po2  = list(lower = c(l_plat, u21,     u22,     -Inf,    -Inf), upper = c(u11,    Inf, Inf, Inf, Inf)),
        po3  = list(lower = c(l_plat, u21,     -Inf,    u32,     -Inf), upper = c(u11,    Inf, u22, Inf, Inf)),
        po4  = list(lower = c(l_plat, u21,     -Inf,    -Inf,    u32),  upper = c(u11,    Inf, u22, u32, Inf)),
        po5  = list(lower = c(l_plat, -Inf,    u22,     -Inf,    -Inf), upper = c(u11,    u21, Inf, Inf, Inf)),
        po6  = list(lower = c(l_plat, -Inf,    l_plat,  u32,     -Inf), upper = c(u11,    u21, u22, Inf, Inf)),
        po9  = list(lower = c(-Inf,   -Inf,    u22,     -Inf,    -Inf), upper = c(l_plat, Inf, Inf, Inf, Inf)),
        po10 = list(lower = c(-Inf,   -Inf,    l_plat,  u32,     -Inf), upper = c(l_plat, Inf, u22, Inf, Inf))
      )
    }

    region_probs <- sapply(regions, function(reg) {
      pmvnorm(lower = reg$lower, upper = reg$upper,
              mean = mu[1:dim_use], corr = res$Sigma[1:dim_use, 1:dim_use])
    })

    if (order == "true") {
      reject_all  <- sum(region_probs[c("po2","po3","po4")])
      reject2not1 <- sum(region_probs[c("po5","po6","po9","po10")])
      result <- c(reject1 = reject_all + region_probs["po1"],
                  reject2 = reject_all + reject2not1,
                  reject2not1 = reject2not1,
                  reject_at_least_one = reject_all + region_probs["po1"] + reject2not1)
    } else {
      reject_all  <- sum(region_probs[c("po2","po3")])
      reject2not1 <- sum(region_probs[c("po6","po7","po10","po11")])
      reject1not2 <- sum(region_probs[c("po1","po4","po5")])
      result <- c(reject1 = reject_all + reject1not2,
                  reject2 = reject_all + reject2not1,
                  reject2not1 = reject2not1,
                  reject_at_least_one = reject_all + reject2not1 + reject1not2)
    }
  }
  return(result)
}

#' Compute Critical Values and Sample Size for a Mono-Combo Seamless Trial
#'
#' Determines the efficacy critical value(s) and per-arm sample size required to control
#' the one-sided family-wise error rate (FWER) at level \code{alpha} and achieve target
#' disjunctive power. Supports four designs selected via \code{platform}, \code{order},
#' and \code{j}.
#'
#' @param alpha Numeric. One-sided FWER target (e.g., \code{0.025}).
#' @param method Character. Alpha-spending method: \code{"OBF"} (O'Brien-Fleming) or
#'   \code{"POC"} (Pocock).
#' @param p0 Numeric. Response rate in the control group under the null.
#' @param p1 Numeric. Target response rate for arm 1 (monotherapy) under the alternative.
#' @param p2 Numeric. Target response rate for arm 2 (combination therapy) under the alternative.
#' @param power Numeric. Target disjunctive power (probability of rejecting at least one null).
#' @param j Integer. Number of stages for MAMS designs (\code{platform = "false"}).
#'   Use \code{j = 2} for 2-stage MAMS and \code{j = 3} for 3-stage MAMS.
#'   Ignored when \code{platform = "true"}. Default is \code{NULL} (treated as \code{3}).
#' @param t Numeric. Planned T2 addition time as a fraction of per-arm enrollment.
#'   For Platform designs this is the fraction at which the combination arm enters
#'   (e.g., \code{t = 0.5}). For MAMS designs, \code{t} is a vector of cumulative stage
#'   fractions (default: equally spaced, \code{seq(1/j, 1, length.out = j)}).
#'   Default is \code{NULL}.
#' @param t_arrive Numeric or \code{NULL}. Actual T2 arrival time (Platform designs only).
#'   \itemize{
#'     \item \code{NULL} (Scenario 1): design is based on the planned \code{t}.
#'     \item Numeric, no \code{N_max} (Scenario 2): \code{t} is forgotten; everything is
#'       re-derived from \code{t_arrive}.
#'     \item Numeric with \code{N_max = "true"} (Scenario 3): total sample size is fixed
#'       from the Scenario 1 baseline plan; only the critical value is re-calibrated.
#'   }
#' @param N_max Character or \code{NULL}. Set to \code{"true"} to invoke Scenario 3 for
#'   Platform designs: the total sample size is held at the value from the baseline plan
#'   (Scenario 1) while the critical value is re-solved for the actual \code{t_arrive}.
#'   Default is \code{NULL}.
#' @param sensitive Character or \code{NULL}. Reserved for internal use; passed through to
#'   \code{mc_pro_analytical}. Default is \code{NULL}.
#' @param platform Character. \code{"true"} for OR-Platform or Platform designs;
#'   \code{"false"} for MAMS designs. Default is \code{"false"}.
#' @param order Character. \code{"true"} for OR-Platform (assumes \eqn{p_2 \ge p_1});
#'   \code{"false"} for Platform (no order assumption). Only relevant when
#'   \code{platform = "true"}. Default is \code{"true"}.
#' @param adjust Character. Set to \code{"true"} to calibrate the critical value
#'   using reproducible cached finite-sample simulations and to verify power after
#'   calibration. This simulation calibration is used for MAMS, planned-entry,
#'   delayed-entry redesign, and fixed-\code{N_max} Platform designs. A fixed-\code{N_max}
#'   design reports whether the power target remains met but does not increase its
#'   sample size.
#'   Default is \code{"false"}.
#' @param calibration_nsim Integer. Number of cached simulations used for each of
#'   the null calibration and alternative power check when simulation calibration
#'   is applied. Default is \code{50000}.
#' @param calibration_seed Integer. Seed used for finite-sample calibration.
#'   Default is \code{271828}.
#'
#' @return
#' \itemize{
#'   \item \strong{Platform designs} (\code{platform = "true"}): a \code{data.frame} with
#'     columns \code{u11}, \code{l1}, \code{u21}, \code{u22}, \code{u32} (stage-specific
#'     critical values), \code{n0} (control arm total), \code{n1} (monotherapy arm total),
#'     \code{n2} (combination arm total), and \code{N_total} (maximum trial sample size).
#'     When \code{adjust = "true"}, these columns contain the final simulation-calibrated
#'     boundaries and sample sizes.
#'   \item \strong{MAMS designs} (\code{platform = "false"}): a one-row
#'     \code{data.frame} with the common stage-specific boundaries \code{u1},
#'     \code{l1}, \code{u2} (and \code{l2}, \code{u3} for a three-stage design),
#'     equal per-arm totals \code{n0}, \code{n1}, \code{n2}, and \code{N_total}.
#'     When \code{adjust = "true"}, these columns contain the final
#'     simulation-calibrated boundaries and sample sizes.
#' }
#'
#' @examples
#' # OR-Platform design (Scenario 1)
#' mc_Design(
#'   alpha = 0.025, method = "POC", p0 = 0.15, p1 = 0.4, p2 = 0.4,
#'   power = 0.8, t = 0.5, platform = "true", order = "true"
#' )
#'
#' # 2-stage MAMS design
#' mc_Design(
#'   alpha = 0.025, method = "POC", p0 = 0.15, p1 = 0.4, p2 = 0.4,
#'   power = 0.8, t = c(0.5, 1), platform = "false", j = 2
#' )
#'
#' # 3-stage MAMS design
#' mc_Design(
#'   alpha = 0.025, method = "POC", p0 = 0.15, p1 = 0.4, p2 = 0.4,
#'   power = 0.8, t = c(1/3, 2/3, 1), platform = "false", j = 3
#' )
#'
#' @importFrom mvtnorm pmvnorm
#' @importFrom MASS fractions
#' @importFrom stats rbinom setNames uniroot
#' @export
# platform="true",  order="true"  -> OR-Platform  (data.frame output)
# platform="true",  order="false" -> Platform     (data.frame output)
# platform="false", j=2           -> 2-stage MAMS (one-row data-frame output)
# platform="false", j=3           -> 3-stage MAMS (one-row data-frame output)
mc_Design <- function(alpha, method, p0, p1, p2, power,
                                   j = NULL, t = NULL, t_arrive = NULL,
                                   N_max = NULL, sensitive = NULL,
                                   platform = "false", order = "true",
                                   adjust = "false",
                                   calibration_nsim = 50000L,
                                   calibration_seed = 271828L) {

  # ---- MAMS path (j=2 or j=3): handled separately with LCM integer constraint ----
  if (platform == "false") {
    stage <- .normalize_mams_t(j, t)
    j <- stage$j
    t <- stage$t

    gcd_fn <- function(a, b) { while (b != 0) { tmp <- b; b <- a %% b; a <- tmp }; a }
    lcm_fn <- function(a, b) { if (a == 0 || b == 0) return(0); abs(a * b) / gcd_fn(a, b) }
    denoms  <- sapply(t, function(x) {
      frac <- MASS::fractions(x)
      as.integer(sub(".*/", "", attr(frac, "fracs")))
    })
    lcm_val <- Reduce(lcm_fn, denoms)

    n_min <- 5; n_max <- 300; finaln <- n_max; prev_n <- 0; c_val <- NA

    while (abs(finaln - prev_n) > 0.1) {
      prev_n <- finaln
      typeI_fn <- function(c1) {
        alpha - .mams_analytical_fixed(
          prev_n, method, p0, p0, p0, c1,
          j = j, t = t, platform = "false", order = order
        )[4]
      }
      c_base <- tryCatch(
        uniroot(typeI_fn, c(1.5, 5), tol = 1e-4)$root,
        error = function(e) uniroot(typeI_fn, c(1, 4), tol = 1e-4)$root
      )
      c_val <- c_base

      power_fn <- function(n_guess) {
        res <- .mams_analytical_fixed(
          n_guess, method, p0, p1, p2, c_val,
          j = j, t = t, platform = "false", order = order
        )
        power - if (p1 == p0) res[3] else res[4]
      }
      n_curr <- tryCatch(
        uniroot(power_fn, c(10, 300), tol = 1e-3)$root,
        error = function(e) finaln
      )
      finaln <- n_curr
    }

    n_final     <- ceiling(finaln / lcm_val) * lcm_val
    power_idx   <- if (p1 == p0) 3L else 4L
    power_check <- .mams_analytical_fixed(
      n_final, method, p0, p1, p2, c_base,
      j = j, t = t, platform = "false"
    )
    if (power_check[power_idx] < power) n_final <- n_final + lcm_val

    c_val <- c_base
    adjust_on <- isTRUE(adjust) || identical(tolower(as.character(adjust)), "true")
    if (adjust_on) {
      repeat {
        cache <- .make_mams_null_cache(
          n = n_final, j = j, t = t, p0 = p0,
          nsim = calibration_nsim, seed = calibration_seed
        )
        calibrated <- .calibrate_mams_cv(
          cache = cache, c_base = c_base, alpha = alpha, method = method
        )
        c_val <- calibrated$cv
        power_check <- .mams_analytical_fixed(
          n_final, method, p0, p1, p2, c_val,
          j = j, t = t, platform = "false", order = order
        )
        power_cache <- .make_mams_cache(
          n = n_final, j = j, t = t,
          p0 = p0, p1 = p1, p2 = p2,
          nsim = calibration_nsim, seed = calibration_seed + 1L
        )
        simulated_power_result <- .mams_rejection_summary_from_cache(
          power_cache, c_val, method
        )
        simulated_power <- if (p1 == p0) {
          simulated_power_result$reject2not1
        } else {
          simulated_power_result$reject_at_least_one
        }
        if (power_check[power_idx] >= power && simulated_power >= power) break
        n_final <- n_final + lcm_val
        if (n_final > n_max) {
          stop("Unable to retain target power after MAMS finite-sample calibration.", call. = FALSE)
        }
      }
    }

    u_final  <- if (method == "OBF") c_val / sqrt(t) else rep(c_val, j)
    l_final  <- c(rep(0, j - 1), u_final[j])

    boundary_values <- unlist(lapply(seq_len(j), function(stage_id) {
      values <- setNames(u_final[stage_id], paste0("u", stage_id))
      if (stage_id < j) {
        values <- c(values, setNames(l_final[stage_id], paste0("l", stage_id)))
      }
      values
    }))
    result <- as.data.frame(as.list(c(
      boundary_values,
      n0 = n_final, n1 = n_final, n2 = n_final,
      N_total = n_final * 3
    )), check.names = FALSE)
    row.names(result) <- NULL
    return(result)
  }

  # ---- Platform path (OR-Platform or Platform): original logic preserved ----

  if (!is.null(N_max) && N_max == "true") {
    t <- .normalize_platform_t(t)
    t_arrive <- .normalize_platform_t(t_arrive)
    plan_res <- mc_Design(alpha, method, p0, p1, p2, power,
                                     t = t, t_arrive = NULL, N_max = NULL,
                                     platform = platform, order = order, adjust = adjust,
                                     calibration_nsim = calibration_nsim,
                                     calibration_seed = calibration_seed)
    N_total <- plan_res$N_total
    n_fixed <- plan_res$n1
    typeI <- function(c1) {
      alpha - .mams_analytical_fixed(
        n_fixed, method, p0, p0, p0, c1,
        j = j, t = t, t_arrive = t_arrive,
        platform = platform, order = order, sensitive = sensitive
      )[4]
    }
    c_base <- tryCatch({
      uniroot(typeI, c(1, 6), tol=1e-4)$root
    }, error = function(e) {
      uniroot(typeI, c(1, 4), tol=1e-4)$root
    })
    c_val <- c_base
    adjust_on <- isTRUE(adjust) || identical(tolower(as.character(adjust)), "true")
    if (adjust_on) {
      null_cache <- .make_platform_null_cache(
        n = n_fixed, t = t, t_arrive = t_arrive, p0 = p0,
        nsim = calibration_nsim, seed = calibration_seed
      )
      calibrated <- .calibrate_platform_cv(
        cache = null_cache, c_base = c_base, alpha = alpha,
        method = method, order = order, sensitive = sensitive
      )
      c_val <- calibrated$cv
      power_idx <- if (p1 == p0) 3L else 4L
      analytical_power_result <- .mams_analytical_fixed(
        n_fixed, method, p0, p1, p2, c_val,
        j = j, t = t, t_arrive = t_arrive,
        platform = platform, order = order, sensitive = sensitive
      )
      power_cache <- .make_platform_cache(
        n = n_fixed, t = t, t_arrive = t_arrive,
        p0 = p0, p1 = p1, p2 = p2,
        nsim = calibration_nsim, seed = calibration_seed + 1L
      )
      simulated_power_result <- .platform_rejection_summary_from_cache(
        power_cache, c_val, method, order, sensitive
      )
      simulated_power <- if (p1 == p0) {
        simulated_power_result$reject2not1
      } else {
        simulated_power_result$reject_at_least_one
      }
    }

    stage_sizes <- .platform_stage_sizes(n_fixed, t, t_arrive)
    if (stage_sizes$N_total != N_total) {
      stop(
        paste0(
          "Integer stage rounding at `t_arrive` gives N_total = ",
          stage_sizes$N_total, " rather than the fixed N_max = ", N_total,
          ". Choose a compatible arrival fraction or revise the integer allocation."
        ),
        call. = FALSE
      )
    }
    upper <- .platform_boundaries(c_val, method, t, t_arrive, sensitive)
    result <- data.frame(
      u11 = upper["u11"], l1 = 0, u21 = upper["u21"],
      u22 = upper["u22"], u32 = upper["u32"],
      n0 = stage_sizes$n0, n1 = stage_sizes$n1, n2 = stage_sizes$n2,
      N_total = N_total,
      row.names = NULL
    )
    return(result)
  }

  else { ### start calculation without same N
    # Scenario 1 (only t): design based on t.
    # Scenario 2 (t + t_arrive, no N_max): forget t, redesign from scratch using t_arrive.
    planned_entry <- is.null(t_arrive)
    if (planned_entry) {
      t <- .normalize_platform_t(t)
    } else {
      t <- t_arrive
    }
    n_min <- 5; n_max <- 300; finaln <- n_max; prev_n <- 0

    while (abs(finaln - prev_n) > 1) {
      prev_n <- finaln
      typeI <- function(c1) {
        alpha - .mams_analytical_fixed(
          prev_n, method, p0, p0, p0, c1,
          j = j, t = t, t_arrive = t_arrive,
          platform = platform, order = order, sensitive = sensitive
        )[4]
      }
      c_base <- tryCatch({
        uniroot(typeI, c(1, 6), tol=1e-4)$root
      }, error = function(e) {
        uniroot(typeI, c(1, 4), tol=1e-4)$root
      })
      c_val <- c_base

      n_lo <- n_min; n_hi <- n_max
      while ((n_hi - n_lo) > 1) {
        mid      <- ceiling((n_lo + n_hi) / 2)
        n1_check <- ceiling(mid * t) + ceiling(mid * (1 - t))
        analytic_result <- .mams_analytical_fixed(
          n1_check, method, p0, p1, p2, c_val,
          j = j, t = t, t_arrive = t_arrive,
          platform = platform, order = order, sensitive = sensitive
        )
        t1 <- if (p1 == p0) analytic_result[3] else analytic_result[4]
        if (t1 >= power) n_hi <- mid else n_lo <- mid
      }
      finaln <- n_hi
    }
  }

  l1 = 0
  n1 <- .next_platform_compatible_n(finaln, t)
  n2 <- n1

  adjust_on <- isTRUE(adjust) || identical(tolower(as.character(adjust)), "true")
  if (adjust_on) {
    power_idx <- if (p1 == p0) 3L else 4L
    repeat {
      null_cache <- .make_platform_null_cache(
        n = n1, t = t, p0 = p0,
        nsim = calibration_nsim, seed = calibration_seed,
        t_arrive = t_arrive
      )
      calibrated <- .calibrate_platform_cv(
        cache = null_cache, c_base = c_base, alpha = alpha,
        method = method, order = order, sensitive = sensitive
      )
      c_val <- calibrated$cv
      analytical_power_result <- .mams_analytical_fixed(
        n1, method, p0, p1, p2, c_val,
        j = j, t = t, t_arrive = t_arrive,
        platform = platform, order = order, sensitive = sensitive
      )
      power_cache <- .make_platform_cache(
        n = n1, t = t, p0 = p0, p1 = p1, p2 = p2,
        nsim = calibration_nsim, seed = calibration_seed + 1L,
        t_arrive = t_arrive
      )
      simulated_power_result <- .platform_rejection_summary_from_cache(
        power_cache, c_val, method, order, sensitive
      )
      simulated_power <- if (p1 == p0) {
        simulated_power_result$reject2not1
      } else {
        simulated_power_result$reject_at_least_one
      }
      if (analytical_power_result[power_idx] >= power && simulated_power >= power) break
      n1 <- .next_platform_compatible_n(n1 + 1L, t)
      n2 <- n1
      if (n1 > n_max) {
        stop("Unable to retain target power after Platform finite-sample calibration.", call. = FALSE)
      }
    }
  }

  upper <- .platform_boundaries(c_val, method, t, t_arrive, sensitive)
  stage_sizes <- .platform_stage_sizes(n1, t, t_arrive)
  result <- data.frame(
    u11 = upper["u11"], l1 = l1, u21 = upper["u21"],
    u22 = upper["u22"], u32 = upper["u32"],
    n0 = stage_sizes$n0, n1 = stage_sizes$n1, n2 = stage_sizes$n2,
    N_total = stage_sizes$N_total,
    row.names = NULL
  )
  return(result)
}

# --- 4. Analytical OC Function ---
mc_analytic_OC <- function(p0, p1, p2, n, cv, t = 0.5, t_arrive = NULL,
                           platform = "true", method = "OBF", order = "true",
                           sensitive = "method1") {

  t_vals <- if(is.null(t_arrive)) t else t_arrive
  total_rows <- length(p1) * length(t_vals)
  results <- matrix(nrow = total_rows, ncol = 8)

  row_idx <- 1
  for (i in 1:length(p1)) {
    for (j in 1:length(t_vals)) {
      current_t_arr <- t_vals[j]

      # Determine CV based on sensitivity method
      if (!is.null(sensitive) && sensitive == "method2") {
        current_cv <- if(length(cv) >= j) cv[j] else cv[1]
      } else {
        current_cv <- cv[1]
      }

      res <- mc_pro_analytical(
        n = n, method = method, p0 = p0, p1 = p1[i], p2 = p2[i],
        cv = current_cv, t = t, t_arrive = current_t_arr,
        platform = platform, order = order, sensitive = sensitive
      )

      results[row_idx, 1] <- p0
      results[row_idx, 2] <- p1[i]
      results[row_idx, 3] <- p2[i]
      results[row_idx, 4] <- current_t_arr
      results[row_idx, 5:8] <- res

      row_idx <- row_idx + 1
    }
  }

  colnames(results) <- c("p0", "p1", "p2", "t_arrive", "reject1", "reject2", "reject2not1", "reject_at_least_one")
  return(as.data.frame(results))
}
