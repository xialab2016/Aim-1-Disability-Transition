# Extracts parameter names and estimates from a fitted msm model
# Assumes 16 transitions and 16 covariate types
f <- function(fit, keep = NULL, ub.sgpv=0) {
  digits <- 5
  qt <- c("12","13","14","15","21","23","24","25",
          "31","32","34","35","41","42","43","45")  # transition labels
  qcov <- c("bq", "bMA", "bOA", "bH", "bL", "bdH", "bdL", "bMS", "bdiagD", 
            "bPS", "bDD", "bTD", "bHM", "bHO", "bLM", "bLO")  # covariate types
  
  qnames <- paste(rep(qcov, each = 16), qt, sep = "")  # full parameter names (e.g., bH12, bdH12, etc.)
  trans <- rep(c("1-2", "1-3", "1-4", "1-5", "2-1", "2-3", "2-4","2-5",
                 "3-1","3-2","3-4","3-5","4-1","4-2","4-3","4-5"), length(qcov))
  beta <- fit$estimates
  beta<-ifelse(beta>=ub.sgpv, beta/exp(beta), beta)#scale 
  beta.se <- sqrt(diag(fit$covmat))
  wald <- (beta / beta.se)^2
  pvalue <- 1 - pchisq(wald, df = 1)
  
  df <- cbind.data.frame(cov = qnames, beta = round(beta, digits))
  return(df)
}
# Computes overall treatment effect (on_high + ddmt_high [+ interaction]) by age group and transition
# Adjusts standard errors using sqrt(2*log(k)/n), where n = transition count from matrix Q

compute_dmt_effects_by_age <- function(coefs, covmat, Q, k = 8,
                                       transitions = c("12", "13", "14", "15", "23", "24","25","34", "35", "45")) {
  results <- data.frame()
  
  for (tr in transitions) {
    from <- as.numeric(substr(tr, 1, 1))
    to   <- as.numeric(substr(tr, 2, 2))
    n <- Q[from, to]
    
    if (n == 0) next
    
    transition_label <- paste0(from, "-", to)  # Converts "12" → "1-2"
    
    # Define terms for each DMT × AgeGroup combination
    params <- list(
      Young_High  = c(paste0("bH", tr), paste0("bdH", tr)),
      Middle_High = c(paste0("bH", tr), paste0("bdH", tr), paste0("bHM", tr)),
      Older_High  = c(paste0("bH", tr), paste0("bdH", tr), paste0("bHO", tr)),
      Young_Low   = c(paste0("bL", tr), paste0("bdL", tr)),
      Middle_Low  = c(paste0("bL", tr), paste0("bdL", tr), paste0("bLM", tr)),
      Older_Low   = c(paste0("bL", tr), paste0("bdL", tr), paste0("bLO", tr))
    )
    
    for (label in names(params)) {
      terms <- params[[label]]
      if (!all(terms %in% names(coefs)) || !all(terms %in% rownames(covmat))) next
      
      est <- sum(coefs[terms])
      
      # Variance calculation including covariances
      var <- sum(sapply(terms, function(x) covmat[x, x]))
      for (i in 1:(length(terms) - 1)) {
        for (j in (i + 1):length(terms)) {
          var <- var + 2 * covmat[terms[i], terms[j]]
        }
      }
      
      se <- sqrt(var) * sqrt(2 * log(k) / n)
      
      results <- rbind(results, data.frame(
        Transition = transition_label,
        AgeGroup = sub("_.+", "", label),
        Treatment = ifelse(grepl("High", label), "Higher-efficacy", "Standard-efficacy"),
        Estimate = round(est, 4),
        StdError = round(se, 4),
        HR = round(exp(est), 4),
        LowerCI = round(exp(est - 1.96 * se), 4),
        UpperCI = round(exp(est + 1.96 * se), 4)
      ))
    }
  }
  
  rownames(results) <- NULL
  return(results)
}

# Computes standard errors, z-scores, p-values and 95% confidence intervals
# Adjusts SE based on transition-specific counts from Q

get_parameter_estimates <- function(coefs, covmat, Q, digits = 4, k = 8) {
  if (is.data.frame(covmat)) covmat <- as.matrix(covmat)
  if (is.null(names(coefs))) stop("Coefficient vector must be named.")
  if (is.null(rownames(covmat))) stop("Covariance matrix must have rownames.")
  
  common_names <- intersect(names(coefs), rownames(covmat))
  results <- data.frame()
  
  for (param in common_names) {
    est <- coefs[param]
    se_raw <- sqrt(covmat[param, param])
    
    # Extract transition index (e.g. bH12 → from=1, to=2)
    tr_match <- regmatches(param, regexec("([0-9])([0-9])$", param))[[1]]
    if (length(tr_match) == 3) {
      from <- as.numeric(tr_match[2])
      to   <- as.numeric(tr_match[3])
      n <- Q[from, to]
    } else {
      n <- 1  # default
    }
    se <- se_raw * sqrt(2 * log(k) / n)
    zval <- est / se
    pval <- 2 * (1 - pnorm(abs(zval)))
    lower <- est - 1.96 * se
    upper <- est + 1.96 * se
    
    results <- rbind(results, data.frame(
      Parameter = param,
      Estimate = round(est, digits),
      StdError = round(se, digits),
      Z = round(zval, digits),
      P = signif(pval, digits),
      LowerCI = round(lower, digits),
      UpperCI = round(upper, digits)
    ))
  }
  
  rownames(results) <- NULL
  return(results)
}


######################################
## Meta-analysis of AMSLS & PROMOTE
######################################

##Function for meta-nalysis
meta.F <- function(est_promote, est_amsls, method = c("fixed", "random")) {
  method <- match.arg(method)
  
  # Merge the data frames by covariate and transition
  merged_data <- merge(
    est_promote, est_amsls,
    by = c("cov", "trans"),
    suffixes = c("_promote", "_amsls")
  )
  
  # Initialize result storage
  results <- data.frame(
    Covariate = merged_data$cov,
    Transition = merged_data$trans,
    HazardRatio_Meta = NA,
    LogHR_Meta = NA,
    LogHR_Meta_SE = NA,
    LogHR_LowerCI = NA,
    LogHR_UpperCI = NA,
    P_Meta = NA,
    Q_Statistic = NA,
    I2_Percent = NA,
    Tau2 = NA,
    PROMOTE = NA,
    AMSLS = NA
  )
  
  # Perform meta-analysis row by row
  for (i in seq_len(nrow(merged_data))) {
    # Extract log HRs and SEs
    b.est <- c(merged_data$beta_promote[i], merged_data$beta_amsls[i])
    se <- c(merged_data$se_promote[i], merged_data$se_amsls[i])
    
    # Fixed-effect weights
    w.fixed <- 1 / se^2
    b.fixed <- sum(w.fixed * b.est) / sum(w.fixed)
    se.fixed <- sqrt(1 / sum(w.fixed))
    
    # Q-statistic and heterogeneity
    q.stat <- sum(w.fixed * (b.est - b.fixed)^2)
    df.q <- length(b.est) - 1
    i2 <- max(0, (q.stat - 10) / q.stat) * 100
    
    # DerSimonian-Laird tau² for random effects
    tau2 <- max(0, (q.stat - df.q) / (sum(w.fixed) - sum(w.fixed^2) / sum(w.fixed)))
    w.random <- 1 / (se^2 + tau2)
    
    # Random effects pooled estimate
    b.random <- sum(w.random * b.est) / sum(w.random)
    se.random <- sqrt(1 / sum(w.random))
    
    # Choose method
    if (method == "fixed") {
      b.meta <- b.fixed
      se.meta <- se.fixed
    } else {
      b.meta <- b.random
      se.meta <- se.random
    }
    
    # CI and p-value
    ci.lower <- b.meta - 1.96 * se.meta
    ci.upper <- b.meta + 1.96 * se.meta
    p.meta <- pchisq((b.meta / se.meta)^2, df = 1, lower.tail = FALSE)
    
    # Hazard ratios
    hr.meta <- exp(b.meta)
    hr.lower <- exp(ci.lower)
    hr.upper <- exp(ci.upper)
    hr_promote <- exp(b.est[1])
    hr_amsls <- exp(b.est[2])
    
    # Store results
    results[i, ] <- c(
      merged_data$cov[i],
      merged_data$trans[i],
      round(hr.meta, 2),
      round(b.meta, 4),
      round(se.meta, 4),
      round(hr.lower, 2),
      round(hr.upper, 2),
      formatC(p.meta, format = "e", digits = 2),
      round(q.stat, 2),
      round(i2, 2),
      round(tau2, 4),
      paste0(round(hr_promote, 2), " (", round(merged_data$se_promote[i], 2), ")"),
      paste0(round(hr_amsls, 2), " (", round(merged_data$se_amsls[i], 2), ")")
    )
  }
  
  return(results)
}


#########################################
## Function to make forest plots
########################################

custom_forest_plot <- function(data, cex.hr, cex.xp, atx, xlim=c(-4.0, 3.2), rectp=3, covp=10, txtcovp_wd=10, stdp=c(1,1), trp=1, hrp=1,  outp=1,alpha=0.2) {
  tryCatch({
    # Ensure necessary columns exist
    required_columns <- c("Covariate", "Transition","PROMOTE", "AMSLS",   "HazardRatio_Meta", "LogHR_Meta", "LogHR_LowerCI", "LogHR_UpperCI", "Outcome", "LogHR_Meta_SE")
    if (!all(required_columns %in% colnames(data))) {
      stop("Data frame must contain required columns: ", paste(required_columns, collapse = ", "))
    }
    
    # Identify groups for each covariate
    covariate_groups <- rle(as.character(data$Covariate))
    group_starts <- cumsum(c(1, head(covariate_groups$lengths, -1)))
    group_lengths <- covariate_groups$lengths
    
    # Prepare display columns
    data$HR_CI_Display <- paste0(
      round(data$HazardRatio_Meta, 2), " [",
      round(data$LogHR_LowerCI, 2), ", ",
      round(data$LogHR_UpperCI, 2), "]"
    )
    data$I2_Outcome <-data$Outcome
    
    # Get row positions
    row_positions <- seq_len(nrow(data))
    
    # Generate the forest plot without the "Study" header
    metafor::forest(
      annotate = FALSE,
      x = data$LogHR_Meta,
      sei = data$LogHR_Meta_SE,
      slab = NA,               # No slab to prevent the "Study" header
      rows = row_positions,
      xlab = "Hazard Ratio",
      refline = 0,
      atransf = exp,
      at = atx,
      digits = c(2, 2),
      cex = cex.hr,
      pch = 15,
      xlim = xlim,
      lwd = 1.2,
      col = "darkgreen",
      header = FALSE           # Suppress the "Study" label
    )
    
    # Add column headers
    text(-covp, max(row_positions) + 2, "Covariate", pos = 4, font = 2, cex = cex.xp)
    text(-stdp[1], max(row_positions) + 2, "PROMOTE", pos = 4, font = 2, cex = cex.xp)
    text(-stdp[2], max(row_positions) + 2, "AMSLS", pos = 4, font = 2, cex=cex.xp)
    text(-trp-0.1, max(row_positions) + 2, "Transition", pos = 4, font = 2, cex=cex.xp)
    text(hrp, max(row_positions) + 2, "HR [95% CI]", pos = 2, font = 2, cex=cex.xp)
    text(outp, max(row_positions) + 2, "State", pos = 4, font = 2, cex=cex.xp)
    
    # Add separators and display values aligned to row positions
    for (i in row_positions) {
      if (i %in% group_starts) {
        lines(x = xlim, y = rep(i - 0.5, 2), lty = "solid", col = "black")
      }
      
      # Display covariate name only for the first row of each group
      #if (i %in% group_starts) {
      # text(-3.0, i, data$Covariate[i], pos = 4, cex = cex.xp, font = 2)
      #}
      # Display covariate name in the middle of the group
      for (g in seq_along(group_starts)) {
        group_center <- group_starts[g] + floor(group_lengths[g] / 2) - 1
        if (i == group_center) {
          wrapped_label <- paste(strwrap(data$Covariate[i], width = txtcovp_wd), collapse = "\n")
          text(-covp, i, wrapped_label, pos = 4, cex = cex.xp, font = 2)
        }
      }
      # Display other columns
      text(-trp, i, data$Transition[i], pos = 4, cex = cex.xp)
      text(-stdp[1] , i, data$PROMOTE[i], pos = 4, cex = cex.xp)
      text(-stdp[2] , i, data$AMSLS[i], pos = 4, cex = cex.xp)
      text(hrp, i, data$HR_CI_Display[i], pos = 2, cex = cex.xp)
      text(outp, i, data$I2_Outcome[i], pos = 4, cex = cex.xp)
    }
    
    # Add background shading for rows with transparency
    for (i in row_positions) {
      if (i %% 2 == 0) {
        rect(
          -rectp, i - 0.5, 3.2, i + 0.5,
          col = adjustcolor("#297B52", alpha.f = alpha),  # 0.5 is 50% transparency
          border = NA
        )
      }
    }
  }, error = function(e) {
    message("Error in forest plot generation: ", e$message)
  })
}

#################################################
##Function to compute HR and treatment for dmt
#effects by group at fixed treatment duration
################################################

compute_HR_CI <- function(age_group, tdelay_value = 1, dmt_duration = 1, dmt_type, coef, vcov_mat) {
  # Initialize coefficient vector
  c_vec <- rep(0, length(coef))
  names(c_vec) <- names(coef)
  
  # Main DMT effects + time on DMT + interaction with age if applicable
  if (dmt_type == "high") {
    c_vec["dmt_high"] <- 1
    c_vec["ddmth"] <- dmt_duration   # Duration effect for high DMT
    if (age_group == "Middle-aged") {
      c_vec["tagecat1:dmt_high"] <- 1
    } else if (age_group == "Older") {
      c_vec["tagecat2:dmt_high"] <- 1
    }
    # No interaction for Young
  } else if (dmt_type == "low") {
    c_vec["dmt_low"] <- 1
    c_vec["ddmtl"] <- dmt_duration   # Duration effect for low DMT
    if (age_group == "Middle-aged") {
      c_vec["tagecat1:dmt_low"] <- 1
    } else if (age_group == "Older") {
      c_vec["tagecat2:dmt_low"] <- 1
    }
    # No interaction for Young
  }
  
  # Treatment delay applies to all age groups
  c_vec["tdelay"] <- tdelay_value
  
  # Compute log(HR), variance, SE
  logHR <- sum(c_vec * coef)
  var_logHR <- as.numeric(t(c_vec) %*% vcov_mat %*% c_vec)
  SE_logHR <- sqrt(var_logHR)
  
  # HR and 95% CI
  HR <- exp(logHR)
  CI_lower <- exp(logHR - 1.96 * SE_logHR)
  CI_upper <- exp(logHR + 1.96 * SE_logHR)
  
  # Create output table
  output <- data.frame(
    Effect = paste(dmt_type, "-", age_group, "- delay:", tdelay_value, "yr - DMT:", dmt_duration, "yr"),
    HR = round(HR, 3),
    SE = round(SE_logHR, 3),
    CI = paste0("[", round(CI_lower, 3), ", ", round(CI_upper, 3), "]"),
    stringsAsFactors = FALSE
  )
  
  return(output)
}

#############################################################
## Function to compute HR and over time for treatment effects
#############################################################

compute_HR <- function(age_group, delay, duration, dmt_type, coef) {
  logHR <- 0
  
  if (dmt_type == "high") {
    logHR <- logHR + coef["dmt_high"] + coef["ddmth"] * duration
    if (age_group == "Middle-aged") {
      logHR <- logHR + coef["tagecat1:dmt_high"]
    } else if (age_group == "Older") {
      logHR <- logHR + coef["tagecat2:dmt_high"]
    }
  } else if (dmt_type == "low") {
    logHR <- logHR + coef["dmt_low"] + coef["ddmtl"] * duration
    if (age_group == "Middle-aged") {
      logHR <- logHR + coef["tagecat1:dmt_low"]
    } else if (age_group == "Older") {
      logHR <- logHR + coef["tagecat2:dmt_low"]
    }
  }
  
  logHR <- logHR + coef["tdelay"] * delay
  
  return(exp(logHR))
}


#####################################################
# Function to extract mapped estimates per transition
#####################################################

var_map <- c(
  bH   = "dmt_high",
  bL   = "dmt_low",
  bdH  = "ddmth",
  bdL  = "ddmtl",
  bHM  = "tagecat1:dmt_high",
  bHO  = "tagecat2:dmt_high",
  bLM  = "tagecat1:dmt_low",
  bLO  = "tagecat2:dmt_low",
  bTD  = "tdelay"
)


extract_coef_list <- function(meta.res, transitions = c("1-4", "1-5", "3-5", "4-5", "2-4", "2-5")) {
  coef_list <- list()
  
  for (tr in transitions) {
    tr_data <- subset(meta.res, Transition == tr)
    
    # Automatically find matching covariate estimates
    coef_vec <- sapply(names(var_map), function(prefix) {
      cov_name <- paste0(prefix, gsub("-", "", tr))  # e.g., bH14 from 1-4
      value <- as.numeric(tr_data$LogHR_Meta[tr_data$Covariate == cov_name])
      if (length(value) == 0) NA_real_ else value
    }, USE.NAMES = TRUE)
    
    # Rename using readable labels
    names(coef_vec) <- unname(var_map)
    
    # Store in list
    coef_list[[tr]] <- coef_vec
  }
  
  return(coef_list)
}


##############################################################################
# Compute transition probabilities/intensities from the fitted msm model
############################################################################

predict_trans_prob_from_data <- function(model,
                                         data,
                                         return_type = c("probability", "intensity"),
                                         covariate_names = NULL,
                                         time_var = "obstime",
                                         state_var = "state",
                                         ci=c("delta", "normal", "bootstrap", "none")) {
  return_type <- match.arg(return_type)
  all_states <- 1:5
  valid_transitions <- list(
    "1" = c(2, 3, 4, 5),
    "2" = c(1, 3, 4, 5),
    "3" = c(1, 2, 4, 5),
    "4" = c(1, 2, 3, 5),
    "5" = c()
  )
  
  # Rename for easier reference
  data <- data %>%
    arrange(id, .data[[time_var]]) %>%
    group_by(id) %>%
    mutate(
      lagstate = lag(.data[[state_var]]),
      lagtime = lag(.data[[time_var]]),
      delta_time = .data[[time_var]] - lagtime,
      rowid = row_number()
    ) %>%
    ungroup() %>%
    filter(!is.na(lagstate) & lagstate != 5 & delta_time > 0)
  
  results <- list()
  
  for (i in 1:nrow(data)) {
    row <- data[i, ]
    from_state <- as.character(row$lagstate)
    to_states <- valid_transitions[[from_state]]
    time <- row$delta_time
    
    covariates <- if (!is.null(covariate_names)) {
      as.list(row[, covariate_names, drop = TRUE])
    } else {
      list()
    }
    
    row_values <- sapply(all_states, function(to) {
      if (to %in% to_states) {
        tryCatch({
          if (return_type == "probability") {
            pmatrix.msm(model, t = time, covariates = covariates)[as.numeric(from_state), to]
          } else {
            qmatrix.msm(model, covariates = covariates, ci=ci)[as.numeric(from_state), to]
          }
        }, error = function(e) NA)
      } else {
        NA
      }
    })
    
    results[[length(results) + 1]] <- cbind(data[i, ], as.list(setNames(row_values, paste0("to_state_", all_states))))
  }
  
  final_df <- bind_rows(results)
  return(final_df)
}

##Usage
#predict_trans_prob_from_data(
#  model = promote.fit,
#  data = msdat1,
#  return_type = "probability",
#  covariate_names = c("dd", "tdelay", "tagecat1", "tagecat2", "on_high", "on_low", "ddmt_high", "ddmt_low", "sex", "tsymp", "progstat"),
#  time_var = "obstime",
#  state_var = "state"
#)


########################################################################
# Function to compute PMR and related metrics for each DMT and age group
########################################################################
compute_pmr_metrics <- function(data, age_breaks = seq(20, 80, by = 10)) {
  # Categorize age
  data <- data %>%
    mutate(
      age_group = cut(rage, breaks = age_breaks, include.lowest = TRUE, right = FALSE)
    )
  
  # Identify forward and reverse transitions (e.g., 2->3 vs 3->2)
  data <- data %>%
    rowwise() %>%
    mutate(
      P_worsen = if (lagstate < state) get(paste0("to_state_", state)) else NA_real_,
      P_improve = if (lagstate > state) get(paste0("to_state_", state)) else NA_real_
    ) %>%
    ungroup()
  
  # Keep only rows with a valid transition probability
  pmr_data <- data %>%
    filter(!is.na(P_worsen) | !is.na(P_improve)) %>%
    mutate(
      dmt_group = case_when(
        on_high == 1 ~ "HE",
        on_low == 1 ~ "LE",
        on_high == 0 & on_low == 0 ~ "UT"
      )
    )
  
  # Summarize by age group and DMT
  summary <- pmr_data %>%
    group_by(age_group, dmt_group) %>%
    summarise(
      P_worsen = mean(P_worsen, na.rm = TRUE),
      P_improve = mean(P_improve, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      PMR = ifelse(P_improve == 0, NA, P_worsen / P_improve),
      CondWorsen = ifelse((P_worsen + P_improve) == 0, NA, P_worsen / (P_worsen + P_improve)),
      TotalTransition = P_worsen + P_improve
    ) %>%
    pivot_wider(
      names_from = dmt_group,
      values_from = c(P_worsen, P_improve, PMR, CondWorsen, TotalTransition, n),
      names_glue = "{.value}_{dmt_group}"
    ) %>%
    mutate(
      RR_HE_LE = P_worsen_HE / P_worsen_LE,
      RR_HE_UT = P_worsen_HE / P_worsen_UT,
      RR_LE_UT = P_worsen_LE / P_worsen_UT,
      RD_HE_LE = P_worsen_HE - P_worsen_LE,
      RD_HE_UT = P_worsen_HE - P_worsen_UT,
      RD_LE_UT = P_worsen_LE - P_worsen_UT,
      PR_HE_LE = 100 * (1 - P_worsen_HE / P_worsen_LE),
      PR_HE_UT = 100 * (1 - P_worsen_HE / P_worsen_UT),
      PR_LE_UT = 100 * (1 - P_worsen_LE / P_worsen_UT),
      RG_HE_over_LE_vs_UT = 100 * (P_worsen_UT - P_worsen_HE) / (P_worsen_UT - P_worsen_LE)
    )
  
  return(summary)
}

#######################################################
## Function to make bubble plot and store PMR metrics
######################################################

bubbleplot <- function(model, data, tp_results=NULL, 
                       return_type = c("probability", "intensity"),
                       covariate_names = NULL,
                       time_var = "obstime",
                       state_var = "state",
                       title = "Transition Metric Bubble Plot",
                       ci=c("delta", "normal", "bootstrap", "none"),
                       saveplotdata=TRUE,
                       metrics = c(
                         "PMR",
                         "PMR_norm",
                         "PMR_percent_reduction",
                         "CondWorsen",
                         "CondWorsen_percent_reduction",
                         "TotalTransition",
                         "TotalWorsenRate",
                         "TotalImproveRate",
                         "PercentReductionVsMax",
                         "PMR_percent_relative"
                       )) {
  
  if (exists("tp_results") && is.data.frame(tp_results)) {
    results <- tp_results
  } else {
    results <- predict_trans_prob_from_data(model, data, return_type,covariate_names, time_var,state_var, ci )
  }
  
  results <- results %>%
    mutate(age_group = cut(tage, breaks = seq(10, 90, by = 10), include.lowest = TRUE)) %>%
    rowwise() %>%
    mutate(
      P_worsen = if (lagstate < state) get(paste0("to_state_", state)) else NA_real_,
      P_improve = if (lagstate > state) get(paste0("to_state_", state)) else NA_real_
    ) %>%
    ungroup()
  
  bubble_data <- results %>%
    filter(!is.na(P_worsen) | !is.na(P_improve), !is.na(dmt_name), !is.na(age_group)) %>%
    group_by(dmt_name, age_group) %>%
    summarise(
      n_patients = n(),
      avg_p_worsen = mean(P_worsen, na.rm = TRUE),
      avg_p_improve = mean(P_improve, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      PMR_norm = ifelse(avg_p_improve > 0, as.numeric(scale(avg_p_worsen/ avg_p_improve)) , NA_real_),
      PMR = ifelse(avg_p_improve > 0, avg_p_worsen/ avg_p_improve , NA_real_),
      CondWorsen = ifelse((avg_p_worsen + avg_p_improve) > 0,
                          avg_p_worsen / (avg_p_worsen + avg_p_improve), NA_real_),
      TotalTransition = ifelse(is.na(avg_p_worsen) & is.na(avg_p_improve),
                               NA_real_, avg_p_worsen + avg_p_improve),
      PMR_percent_reduction = ifelse(avg_p_improve > 0,
                                     100 * (1 - (avg_p_worsen / avg_p_improve)), NA_real_),
      CondWorsen_percent_reduction = ifelse((avg_p_worsen + avg_p_improve) > 0,
                                            100 * (1 - (avg_p_worsen / (avg_p_worsen + avg_p_improve))), NA_real_),
      TotalWorsenRate = avg_p_worsen,
      TotalImproveRate = avg_p_improve
    )
  
  max_pmr <- max(bubble_data$PMR, na.rm = TRUE)
  
  bubble_data <- bubble_data %>%
    mutate(
      PMR_percent_relative = 100 * (1 - PMR / max_pmr)
    )
  
  max_worsen <- max(bubble_data$avg_p_worsen, na.rm = TRUE)
  
  
  bubble_data <- bubble_data %>%
    mutate(
      PercentReductionVsMax = 100 * (1 - avg_p_worsen / max_worsen)
    )
  
  # Filter to selected metrics only
  bubble_long <- bubble_data %>%
    pivot_longer(
      cols = all_of(metrics),
      names_to = "metric_name",
      values_to = "metric_value"
    )
  
  if(saveplotdata==TRUE){save(bubble_long, file="bubblePlotdata.rdata")}
  
  # Plot
  ggplot(bubble_long, aes(x = age_group, y = metric_value, size = n_patients, fill = dmt_name)) +
    geom_point(shape = 21, color = "black", alpha = 0.75) +
    facet_wrap(~ metric_name, scales = "free_y") +
    scale_size_continuous(range = c(2, 12)) +
    labs(
      title = title,
      x = "Age Group",
      y = "Metric Value",
      size = "Number of Patients",
      fill = "DMT Name"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


##Usage
#bubbleplot(
#  model = promote.fit, #fitted msm model
#  data = msdat, #date containing covariates used to fit the model. Names of covariate shoild be same as in model
#  title = "Promote PMR and CondWorsen", #titel of the plots
#  metrics = c("PMR", "CondWorsen", "CondWorsen_percent_reduction", "PMR_percent_relative", "TotalWorsenRate", "TotalImproveRate")  #metrics to be ploted
#)



