#' Unpooled Two-Sample Proportion Z-Test with Boundary Safeguard
#'
#' Computes the ordinary unpooled z-statistic for comparing two proportions.
#' Only when its plug-in standard error is zero or undefined, one success and
#' one failure are added to each group and the statistic is recomputed.
#' Returns 0 when either original sample size is non-positive or the fallback
#' standard error remains zero or undefined.
#'
#' @param success1 Numeric. Number of successes in group 1 (treatment).
#' @param success2 Numeric. Number of successes in group 2 (control).
#' @param n1 Numeric. Sample size for group 1.
#' @param n2 Numeric. Sample size for group 2.
#'
#' @return Numeric. Ordinarily, the z-statistic
#'   \eqn{(\hat{p}_1 - \hat{p}_2) / SE}, where
#'   \eqn{SE = \sqrt{\hat{p}_1(1-\hat{p}_1)/n_1 +
#'   \hat{p}_2(1-\hat{p}_2)/n_2}}. If this standard error is zero or
#'   undefined, the same formula is evaluated using
#'   \eqn{\tilde{p}_k=(x_k+1)/(n_k+2)} and adjusted sizes \eqn{n_k+2}.
#'   Returns 0 if \code{n1 <= 0}, \code{n2 <= 0}, or \code{SE <= 0}.
#'
#' @examples
#' proportion_test(success1 = 8, success2 = 4, n1 = 18, n2 = 18)
#'
#' @export
proportion_test <- function(success1, success2, n1, n2) {
  if (n1 <= 0 || n2 <= 0) return(0)
  p1 <- success1 / n1
  p2 <- success2 / n2
  se <- sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
  if (is.finite(se) && se > 0) return((p1 - p2) / se)

  # Conditional Agresti--Caffo-style add-one safeguard for SE = 0.
  n1_adjusted <- n1 + 2
  n2_adjusted <- n2 + 2
  p1 <- (success1 + 1) / n1_adjusted
  p2 <- (success2 + 1) / n2_adjusted
  se <- sqrt(
    p1 * (1 - p1) / n1_adjusted +
      p2 * (1 - p2) / n2_adjusted
  )
  if (!is.finite(se) || se <= 0) return(0)
  return((p1 - p2) / se)
}

#' Simulation Engine for Mono-Combo Seamless Trial
#'
#' Runs Monte Carlo simulations for a single scenario of the mono-combo seamless trial and
#' returns rejection probabilities, early stopping probabilities, and average sample sizes.
#' This is the internal engine called by \code{mc_OC}; users typically do not call it directly.
#'
#' @param n Numeric. Per-arm per-stage sample size.
#' @param method Character. Alpha-spending method: \code{"OBF"} (O'Brien-Fleming) or
#'   \code{"POC"} (Pocock).
#' @param p0 Numeric. Response rate in the control group.
#' @param p1 Numeric. Response rate for arm 1 (monotherapy).
#' @param p2 Numeric. Response rate for arm 2 (combination therapy).
#' @param cv Numeric. Critical value for the efficacy boundary.
#' @param j Integer. Number of stages for MAMS designs (\code{platform = "false"}).
#'   Ignored when \code{platform = "true"}. Default is \code{3}.
#' @param t Numeric. Planned T2 addition time as a fraction of per-arm enrollment.
#'   Default is \code{0.5}.
#' @param t_arrive Numeric or \code{NULL}. Actual T2 arrival time. If \code{NULL},
#'   \code{t} is used. Default is \code{NULL}.
#' @param sensitive Character or \code{NULL}. Boundary selection method when
#'   \code{t_arrive} differs from \code{t}: \code{"method1"} uses boundaries based on the
#'   original planned \code{t}; other values use boundaries re-derived from \code{t_arrive}.
#'   Default is \code{NULL}.
#' @param platform Character. \code{"true"} for OR-Platform or Platform designs;
#'   \code{"false"} for MAMS designs. Default is \code{"true"}.
#' @param seed Integer. Random seed for reproducibility. Default is \code{100}.
#' @param order Character. \code{"true"} for OR-Platform (order restriction \eqn{p_2 \ge p_1});
#'   \code{"false"} for Platform. Default is \code{"true"}.
#' @param nsim Integer. Number of simulation replicates. Default is \code{10000}.
#'
#' @return A named numeric vector of operating characteristics. For OR-Platform
#'   (\code{order = "true"}): \code{reject1}, \code{reject_2not1}, \code{reject2},
#'   \code{reject_at_least_one}, \code{PES_Eff}, \code{PES_Fut}, \code{PES},
#'   \code{aveN0}, \code{aveN1}, \code{aveN2}, \code{aveTotalN}.
#'   For Platform (\code{order = "false"}): same columns plus \code{reject1not2}.
#'
#' @seealso \code{\link{mc_OC}} which loops over multiple scenarios and wraps this function.
#'
#' @export
mc_pro <- function(n, method=c("OBF","POC"), p0, p1, p2, cv, j=3, t=0.5,
                   t_arrive=NULL, sensitive=NULL, platform="true",
                   seed=100, order="true", nsim=10000) {

  set.seed(seed)
  method <- match.arg(method)

  # --- 1. Parameter Setup ---
  if (platform == "true") {

    power=reject_all=reject2not1=reject1not2=reject2=reject1=reject_at_least_one=0;
    PES=PES_Fut=PES_Eff=0;
    po1=po3=po2=po4=po5=po6=po7=po8=po9=po10=po11=po12=po13=po14=po15=po16=pt=0
    nosele_futi=0
    aveN0=aveN1=aveN2= aveTotalN=0
    result <- vector(length = 11)

    u1=cv;
    l1=0;

    # Platform designs typically assume 2 stages of decision per arm
    if (is.null(t_arrive)) {
      r10=r11=r12=t; r20=r21=r22=1-t; r30=r32=t
    } else {
      r10=r11=r12=t_arrive; r20=r21=r22=1-t_arrive; r30=r32=0.5*(t+t_arrive)
    }
  }
  else {    # MAMS Design Setup
    stage <- .normalize_mams_t(j, t)
    j <- stage$j
    t_vec <- stage$t

    n_cum <- round(n * t_vec)
    n_inc <- c(n_cum[1], diff(n_cum))
    u_mams <- if(method == "OBF") cv / sqrt(t_vec) else rep(cv, j)
    l_mams <- c(rep(0, j-1), u_mams[j])
  }

  # --- 2. Data Pre-generation ---
  # For MAMS logic
  if (platform == "false") {
    Y0 <- matrix(rbinom(nsim * j, rep(n_inc, each = nsim), p0), nrow=nsim)
    Y1 <- matrix(rbinom(nsim * j, rep(n_inc, each = nsim), p1), nrow=nsim)
    Y2 <- matrix(rbinom(nsim * j, rep(n_inc, each = nsim), p2), nrow=nsim)
  }

  # Initial Counters
  po <- numeric(30) # Path counters
  aveN0 <- aveN1 <- aveN2 <- 0

  # --- 3. Simulation Loop ---
  for (i in 1:nsim) {

    if (platform == "true") {


      # --- PLATFORM LOGIC (Original) ---
      y0s1 <- rbinom(1, round(n*r10), p0); y1s1 <- rbinom(1, round(n*r11), p1); y2s1 <- rbinom(1, round(n*r12), p2)
      t11 <- proportion_test(y1s1, y0s1, round(n*r11), round(n*r10))

      y0s2 <- rbinom(1, round(n*r20), p0); y1s2 <- rbinom(1, round(n*r21), p1); y2s2 <- rbinom(1, round(n*r22), p2)
      t21 <- proportion_test(y1s1+y1s2, y0s1+y0s2, round(n*(r11+r21)), round(n*(r10+r20)))

      # Platform specific T2 comparison (concurrent vs cumulative)
      t22 <- proportion_test(y2s2, y0s2, round(n*r22), round(n*r20))

      y0s3 <- rbinom(1, round(n*r30), p0); y2s3 <- rbinom(1, round(n*r32), p2)
      t32 <- proportion_test(y2s2+y2s3, y0s2+y0s3, round(n*(r22+r32)), round(n*(r20+r30)))

      yc=y0s1+2*y0s2+y0s3
      yt=y1s1+y1s2+y2s2+y2s3

      n_E <- (r11 + r21 + r22 + r32) * n
      n_C <- (r10 + 2*r20 + r30) * n

      pE <- yt / n_E
      pC <- yc / n_C
      varE <- (pE * (1 - pE)) / n_E

      # Control variance (weighted to account for 2*y0stag2)
      n_c1 <- n * r10
      n_c2 <- n * r20
      n_c3 <- n * r30

      varC_weighted <- (
        (n_c1 + 4 * n_c2 + n_c3) * pC * (1 - pC)
      ) / (n_C^2)

      se_ct <- sqrt(varE + varC_weighted)
      if (!is.finite(se_ct) || se_ct <= 0) {
        n_E_adjusted <- n_E + 2
        n_C_adjusted <- n_C + 2
        pE <- (yt + 1) / n_E_adjusted
        pC <- (yc + 1) / n_C_adjusted
        varE <- pE * (1 - pE) / n_E_adjusted
        varC_weighted <- (
          (n_c1 + 4 * n_c2 + n_c3 + 2) * pC * (1 - pC)
        ) / (n_C_adjusted^2)
        se_ct <- sqrt(varE + varC_weighted)
      }
      tct <- if (is.finite(se_ct) && se_ct > 0) (pE - pC) / se_ct else 0

      # 1. Handle the simplest case first: POC method is always constant
      if (method == "POC") {
        u11 <- u21 <- u22 <- u32 <- u2 <- cv

      }
      else if (method == "OBF") {

        # 2. Determine which time parameter to use for the base calculation
        # Use 't' if 'sensitive' is "method1" OR if 't_arrive' is NULL
        use_t_only <- identical(sensitive, "method1") || is.null(t_arrive)

        if (use_t_only) {
          # Standard OBF logic
          u2  <- cv * sqrt(t)
          u11 <- cv
          u21 <- u2
          u22 <- u11 * sqrt(t) / sqrt(1 - t)
          u32 <- u2

        } else {
          # Arrival-time specific OBF logic
          u2  <- cv * sqrt(t_arrive)
          u11 <- cv
          u21 <- u2
          # Simplified the fraction for readability
          u22 <- u2 * sqrt(((1 - t_arrive) + 0.5 * (t + t_arrive)) / (1 - t_arrive))
          u32 <- u2
        }
      }

      if (order=="true"){
        ####  T1 show efficiency at the first stage  ###
        ########### 1################## select t1 and t2 ####
        if(t11>=u11){                                    ##### 1. Early stop due to efficiency   ###
          po1=po1+1/nsim }                               #### select T1


        ##### T1 did not show efficiency at the first stage####
        if (t11< u11 & t11>=l1) { ##### T1 move to stage 2 and T2 is newly added in stage 2

          #### Following 5 conditions based on t22 and t21
          if (t21>=u21) {
            if (t22>=u22) { po2=po2+1/nsim   }            #### 2. Select t1 and t2

            if (t22<u22) {###proceed T2 to stage 3
              if (t32>=u32){po3=po3+1/nsim}         ### 3 ### select t1 and t2
              if(t32<u32){
                if (tct>u32) {po4=po4+1/nsim} ### 4 ###  select t1 and t2
                else {po13=po13+1/nsim}  ##### 13 #### both futile
              }
            }
          }

          if (t21< u21) {
            if ( t22>=u22)  { po5=po5+ 1/nsim}  ####5 Early stop, select t2 only ####

            if ( t22< u22&t22>=l1) {#### Proceed T2 to next stage

              if (t32>=u32){po6=po6+1/nsim} ### 6 ### select t2
              if (t32<u32){po7=po7 +1/nsim} }###7 no selection duet to futility

            if ( t22<l1)  {po8=po8+1/nsim } ###8 early stop due to futility and no selection##3

          }
        }


        ##### When T1 show futility at the first stage####
        if (t11<l1){

          if ( t22>=u22){  po9=po9+1/nsim}    ### 9 Early stop select T2

          if (t22<u22& t22>=l1) { #### proceed T2 to next stage
            if (t32>=u32){po10=po10+1/nsim}     ### 10 ### select t2
            if (t32< u32){ po11=po11+1/nsim }} ####11 No selection due to futility

          if (t22<l1)
          {  po12=po12+1/nsim}  ###12 Early termination and no selection
        }
      }


      if (order=="false"){
        ####  T1 show efficiency at the first stage  ###
        ########### 1################## select t1 and t2 ####
        if(t11>=u11){                                    ##### 1. Early stop due to efficiency   ###
          po1=po1+1/nsim }                               #### select T1


        ##### T1 did not show efficiency at the first stage####
        if (t11< u11& t11>=l1) { ##### T1 move to stage 2 and T2 is newly added in stage 2
          ### Test base on only concurrent control data##

          #### Following 5 conditions based on t22 and t21
          if (t21>=u21) {
            if (t22>=u22) { po2=po2+1/nsim   }                     #### 2. Select t1 and t2

            if (t22>=l1 & t22<u22) {###proceed T2 to stage 3

              if (t32>=u32){po3=po3+1/nsim}         ### 3 ### select t1 and t2
              if (t32<u32){po4=po4+1/nsim}}       ### 4 ###  Select T1

            if (t22< l1){                       #### 5 Early stop trial and select t1}
              po5=po5+ 1/nsim}
          }


          if (t21< u21) {
            if ( t22>=u22)  { po6=po6+ 1/nsim}  ####6 Early stop, select t2 only ####

            if ( t22< u22 &t22>=l1) {#### Proceed T2 to next stage

              if (t32>=u32){po7=po7+1/nsim} ### 7 ### select t2
              if (t32<u32){po8=po8 +1/nsim} }###8 no selection duet to futility

            if ( t22<l1)  {po9=po9+1/nsim } ###9 early stop due to futility and no selection##3

          }
        }

        ##### When T1 show futility at the first stage####
        if (t11<l1){

          ### Drop T1 and add T2
          ### Test base on concurrent control data##
          if ( t22>=u22){  po10=po10+1/nsim}    ### 10 Early stop select T2

          if (t22<u22& t22>=l1) { #### proceed T2 to next stage

            if (t32>=u32){po11=po11+1/nsim}     ### 11 ### select t2
            if (t32< u32){ po12=po12+1/nsim }} ####12 No selection due to futility

          if (t22<l1)
          {  po13=po13+1/nsim}  ###13 Early termination and no selection
        }
      }


    } else {
      # --- MAMS LOGIC (New Optimized) ---
      y0c <- Y0[i, 1]; y1c <- Y1[i, 1]; y2c <- Y2[i, 1]
      z11 <- proportion_test(y1c, y0c, n_cum[1], n_cum[1])
      z12 <- proportion_test(y2c, y0c, n_cum[1], n_cum[1])

      if (j == 3) {

        # Path Mapping for j=3
        if (z11 >= u_mams[1]) {
               if (z12 >= u_mams[1]) po[1] <- po[1] + 1/nsim ##p1 ### select T1 and T2
               else po[2] <- po[2] + 1/nsim } ### select T1


        else if (z11 >= l_mams[1]) {
               if (z12 >= u_mams[1]) po[3] <- po[3] + 1/nsim ### select T2
               else if (z12 >= l_mams[1]) {
                     y0c <- y0c + Y0[i, 2]; y1c <- y1c + Y1[i, 2]; y2c <- y2c + Y2[i, 2]
                     z21 <- proportion_test(y1c, y0c, n_cum[2], n_cum[2])
                      z22 <- proportion_test(y2c, y0c, n_cum[2], n_cum[2])
                     if (z21 >= u_mams[2]) {
                             if (z22 >= u_mams[2]) po[4] <- po[4] + 1/nsim ### select both
                             else po[5] <- po[5] + 1/nsim} ### select T1
                     else if (z21 >= l_mams[2]) {
                           if (z22 >= u_mams[2]) po[6] <- po[6] + 1/nsim ### select T2
                           else if (z22 >= l_mams[2]) {
                                 y0c <- y0c+Y0[i,3]; y1c <- y1c+Y1[i,3]; y2c <- y2c+Y2[i,3]
                                 z31 <- proportion_test(y1c, y0c, n_cum[3], n_cum[3])
                                 z32 <- proportion_test(y2c, y0c, n_cum[3], n_cum[3])
                                    if (z31 >= u_mams[3] && z32 >= u_mams[3]) po[7] <- po[7] + 1/nsim ### select T1, T2
                                    else if (z31 >= u_mams[3]) po[8] <- po[8] + 1/nsim### select T1
                                    else if (z32 >= u_mams[3]) po[9] <- po[9] + 1/nsim ### select  T2
                                          else po[10] <- po[10] + 1/nsim  } ### select none
                           else {   y0c <- y0c+Y0[i,3]; y1c <- y1c+Y1[i,3]
                                  if (proportion_test(y1c, y0c, n_cum[3], n_cum[3]) >= u_mams[3]) po[11] <- po[11] + 1/nsim ### select T1
                                          else po[12] <- po[12] + 1/nsim  }  } ### select none
                  else {### z21<l
                         if (z22 >= u_mams[2]) po[13] <- po[13] + 1/nsim ### select  T2
                         else if (z22 >= l_mams[2]) {
                                  y0c <- y0c+Y0[i,3]; y2c <- y2c+Y2[i,3]
                                   if (proportion_test(y2c, y0c, n_cum[3], n_cum[3]) >= u_mams[3]) po[14] <- po[14] + 1/nsim  ### select  T2
                                      else po[15] <- po[15] + 1/nsim   } ### select none
                         else {po[16] <- po[16] + 1/nsim}  ### select none
                        } ### end of z21<l
                   }### end of l<z12<u1
               else { # T2 drops S1
                    y0c <- y0c+Y0[i,2]; y1c <- y1c+Y1[i,2]
                     if (proportion_test(y1c, y0c, n_cum[2], n_cum[2]) >= u_mams[2]) po[17] <- po[17] + 1/nsim ### select T1
                     else if(proportion_test(y1c, y0c, n_cum[2], n_cum[2]) >= l_mams[2]){
                                   y0c <- y0c + Y0[i,3]; y1c <- y1c + Y1[i,3]
                                   if  (proportion_test(y1c, y0c, n_cum[3], n_cum[3]) >= u_mams[3]) po[18] <- po[18] + 1/nsim ### select T1
                                  else po[19] <- po[19] + 1/nsim  ### select none
                     }
                     else{po[20] <- po[20] + 1/nsim } ### select none

                      }
                   } ### end of l<z11<u1
        else { # T1 drops S1
          if (z12 >= u_mams[1]) po[21] <- po[21] + 1/nsim  ### select T2
          else if (z12 >= l_mams[1]) {
               y0c <- y0c+Y0[i,2]; y2c <- y2c+Y2[i,2]
                z22 <- proportion_test(y2c, y0c, n_cum[2], n_cum[2])
                      if (z22 >= u_mams[2]) po[22] <- po[22] + 1/nsim ### select T2
                      else if (z22 >= l_mams[2]) {
                             y0c <- y0c+Y0[i,3]; y2c <- y2c+Y2[i,3]
                              if (proportion_test(y2c, y0c, n_cum[3], n_cum[3]) >= u_mams[3]) po[23] <- po[23] + 1/nsim ### select T2
                               else po[24] <- po[24] + 1/nsim ### select none
            } else po[25] <- po[25] + 1/nsim  ### select none
          } else po[26] <- po[26] + 1/nsim  ### select none
        }


      }
      else {
        # --- j=2 MAMS Logic ---
        if (z11 >= u_mams[1]) {
          if (z12 >= u_mams[1]) po[1] <- po[1] + 1/nsim else po[2] <- po[2] + 1/nsim
        } else if (z11 >= l_mams[1]) {
          if (z12 >= u_mams[1]) po[3] <- po[3] + 1/nsim
          else if (z12 >= l_mams[1]) {
            y0c <- y0c+Y0[i,2]; y1c <- y1c+Y1[i,2]; y2c <- y2c+Y2[i,2]
            z21 <- proportion_test(y1c, y0c, n_cum[2], n_cum[2])
            z22 <- proportion_test(y2c, y0c, n_cum[2], n_cum[2])
            if (z21 >= u_mams[2] && z22 >= u_mams[2]) po[4] <- po[4] + 1/nsim
            else if (z21 >= u_mams[2]) po[5] <- po[5] + 1/nsim
            else if (z22 >= u_mams[2]) po[6] <- po[6] + 1/nsim else po[7] <- po[7] + 1/nsim
          } else { # T2 drops S1
            y0c <- y0c+Y0[i,2]; y1c <- y1c+Y1[i,2]
            if (proportion_test(y1c, y0c, n_cum[2], n_cum[2]) >= u_mams[2]) po[8] <- po[8] + 1/nsim else po[9] <- po[9] + 1/nsim
          }
        } else {
          if (z12 >= u_mams[1]) po[10] <- po[10] + 1/nsim
          else if (z12 >= l_mams[1]) {
            y0c <- y0c+Y0[i,2]; y2c <- y2c+Y2[i,2]
            if (proportion_test(y2c, y0c, n_cum[2], n_cum[2]) >= u_mams[2]) po[11] <- po[11] + 1/nsim else po[12] <- po[12] + 1/nsim
          } else po[13] <- po[13] + 1/nsim
        }
      }
    }
  }

  # --- 4. Result Aggregation ---
  if (platform == "true") {
    if (order=="true"){
      PES_Eff=po1+po2+po5+po9
      PES_Fut=po8+po12
      PES=PES_Eff+PES_Fut

      reject_all=po2+po3+po4
      reject2not1=po5+po6+po9+po10
      reject2=reject_all+reject2not1
      reject1=reject_all+po1

      reject_at_least_one=reject_all+reject2not1+po1


      nosele_futi=po7+po8+po11+po12+po13
      aveN0=n*r11+  # from z11
        r20*n*(1-po1)+ #from z22
        # (1-t)*n(po2+po3+po4+p13+po5+po6+po7+po8)+# from z21
        r30*n*(po3+po4+po13+po6+po7+po10+po11) # from z32

      aveN1=n*r11+  # from z11
        r21*n*(po2+po3+po4+po13+po5+po6+po7+po8) #from z21

      aveN2= r22*n*(1-po1)+ # from z22
        r32*n*(po3+po4+po13+po6+po7+po10+po11) # from z32

      aveTotalN=2*n*r11+  # from z11
        3*r21*n*(po2+po3+po4+po13+po5+po6+po7+po8)+ #from z21
        2*r22*n*(po9+po10+po11+po12)+# from z22
        2*r32*n*(po3+po4+po13+po6+po7+po10+po11) # from z32

      #aveN0=n*(1+1-po1+po3+po4+po13+po6+po7+po10+po11)
      #aveN1=n*(1+1-po1-po9-po10-po11-po12)
      #aveN2=n*(1-po1+po3+po4+po13+po6+po7+po10+po11)
      #aveTotalN=2*n+ 2*n*(po9+po12)+3*n*(po2+po5+po8 ) + 4*n*(po10+po11)+5*n*(po3+po4+po13+po7+po6)
      result=c(reject1,reject2not1,  reject2,reject_at_least_one,PES_Eff,PES_Fut,PES, aveN0, aveN1, aveN2, aveTotalN)
      result=setNames(result, c ("reject1", "reject_2not1" , "reject2"  ,"reject_at_least_one"," PES_Eff" ,  "PES_Fut",   "PES"  ," aveN0",  "aveN1",  " aveN2", "aveTotalN"))
    }

    if (order=="false"){
      PES_Eff=po1+po2+po5+po6+po10
      PES_Fut=po9+po13
      PES=PES_Eff+PES_Fut

      reject_all=po2+po3
      reject2not1=po6+po7+po10+po11
      reject1not2=po4+po5+po1
      reject_at_least_one=reject_all+reject2not1+reject1not2
      reject2=reject_all+reject2not1
      reject1=reject_all+reject1not2

      nosele_futi=po8+po9+po12+po13
      # totalp=rejctatleastone+nosele_futi
      #aveN0=n*(1+1-po1+po3+po4+po8+po7+po12+po11)
      #aveN1=n*(1+1-po1-po10-po11-po12-po13)
      #aveN2=n*(1-po1+po3+po4+po8+po7+po12+po11)
      #aveTotalN=2*n+ 2*n*(po10+po13)+3*n*(po2+po5+po6+po9 ) + 4*n*(po11+po12)+5*n*(po3+po4+po7+po8)

      aveN0= n*r10+ # from z11
        r20*n*(1-po1)+ #from z22
        n*r30*(po3+po4+po8+po7+po12+po11) # from z32

      aveN1=n*r11+ # from z11
        n*r21*(po2+po3+po4+po8+po7+po5+po6+po9) #from z21

      aveN2=r22*n*(1-po1)+ #from z22
        n*r32*(po3+po4+po8+po7+po12+po11) # from z32

      aveTotalN=2*n*r11+ # from z11
        3* n*r22*(po2+po3+po4+po8+po7+po5+po6+po9)+ #from both z21 and z22
        2*r22*n*(po10+po11+po12+po13)+ #from z22
        2*n*r32*(po3+po4+po8+po7+po12+po11) # from z32

      result=c(reject1,reject2not1, reject2, reject_at_least_one, reject1not2, PES_Eff,PES_Fut,PES, aveN0, aveN1, aveN2,aveTotalN)
      result=setNames(result, c ("reject1", "reject2not1" , "reject2"  ,"reject_at_least_one","reject1not2", " PES_Eff " ,  "PES_Fut  ",   "PES"  ," aveN0",  "aveN1",  " aveN2", "aveTotalN"))
    }
  }

    else {
      # Result Aggregation for MAMS designs
      if (j == 3) {
        reject1 = sum(po[c(1, 2, 4, 5, 7, 8, 11, 17, 18)]); reject2 = sum(po[c(1, 3, 4, 6, 7, 9, 13, 14, 21, 22, 23)])
        reject_all = sum(po[c(1, 4, 7)]); reject1not2 = reject1 - reject_all; reject2not1 = reject2 - reject_all; reject_at_least_one = reject1 + reject2 - reject_all
        PES_Eff = sum(po[c(1:6, 13, 17, 21, 22)]); PES_Fut = sum(po[c(16, 20, 25, 26)]); PES = PES_Eff + PES_Fut
        aveN0 = n/3 + n/3*(1-sum(po[c(1,2,3,21,26)])) + n/3*(sum(po[c(7:12, 14, 15, 18, 19, 23, 24)]))
        aveN1 = n/3 + n/3*(1-sum(po[c(1,2,3,21:26)])) + n/3*(sum(po[c(7:12, 18, 19)]))
        aveN2 = n/3 + n/3*(1-sum(po[c(1,2,3,17:21, 26)])) + n/3*(sum(po[c(7:10, 14, 15, 23, 24)]))
        aveTotalN = aveN0 + aveN1 + aveN2
        res = c(reject1, reject2not1, reject2, reject_at_least_one, reject1not2, PES_Eff, PES_Fut, PES, aveN0, aveN1, aveN2, aveTotalN)
        result=setNames(res, c ("reject1", "reject2not1" , "reject2"  ,"reject_at_least_one","reject1not2", " PES_Eff " ,  "PES_Fut  ",   "PES"  ," aveN0",  "aveN1",  " aveN2", "aveTotalN"))

      } else {
        reject1 = sum(po[c(1, 2, 4, 5, 8)]); reject2 = sum(po[c(1, 3, 4, 6, 10, 11)])
        reject_all = sum(po[c(1, 4)]); reject1not2 = reject1 - reject_all; reject2not1 = reject2 - reject_all; reject_at_least_one = reject1 + reject2 - reject_all
        PES_Eff = po[1]+po[2]+po[3]+po[10]; PES_Fut = po[13]; PES = PES_Eff + PES_Fut
        interim_fraction <- t_vec[1]
        aveN0 = n*interim_fraction + n*(1-interim_fraction)*(1-po[1]-po[2]-po[3]-po[10]-po[13])
        aveN1 = n*interim_fraction + n*(1-interim_fraction)*(po[4]+po[5]+po[6]+po[7]+po[8]+po[9])
        aveN2 = n*interim_fraction + n*(1-interim_fraction)*(po[4]+po[5]+po[6]+po[7]+po[11]+po[12])
        aveTotalN = aveN0 + aveN1 + aveN2
        res = c(reject1, reject2not1, reject2, reject_at_least_one, reject1not2, PES_Eff, PES_Fut, PES, aveN0, aveN1, aveN2, aveTotalN)
        result=setNames(res, c ("reject1", "reject2not1" , "reject2"  ,"reject_at_least_one","reject1not2", " PES_Eff " ,  "PES_Fut  ",   "PES"  ," aveN0",  "aveN1",  " aveN2", "aveTotalN"))
      }
    }



  return(result)
}
