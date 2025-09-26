#' Prepare feature-level data for protein-level summarization
#'
#' @param input feature-level data processed by dataProcess subfunctions
#' @param method summarization method - `summaryMethod` parameter of the dataProcess function
#' @param impute if TRUE, censored missing values will be imputed - `MBimpute`
#' parameter of the dataProcess function
#' @param censored_symbol censored missing value indicator - `censoredInt`
#' parameter of the dataProcess function
#' @param remove_uninformative_feature_outlier if TRUE, features labeled as
#' outlier of uninformative by the MSstatsSelectFeatures function will not be
#' used in summarization
#'
#' @return data.table
#'
#' @export
#'
#' @examples
#' raw <- DDARawData
#' method <- "TMP"
#' cens <- "NA"
#' impute <- TRUE
#' MSstatsConvert::MSstatsLogsSettings(FALSE)
#' input <- MSstatsPrepareForDataProcess(raw, 2, NULL)
#' head(input)
#'
MSstatsPrepareForSummarization <- function(input, method, impute, censored_symbol,
                                           remove_uninformative_feature_outlier) {
  ABUNDANCE <- feature_quality <- is_outlier <- PROTEIN <- NULL

  label <- data.table::uniqueN(input$LABEL) == 2
  if (label) {
    # removed because not efficient 
    # input[, ref := factor(ifelse(LABEL == "L", RUN, 0))]
    input[, ref := NA_character_]
    input[LABEL == "L", ref := RUN]
    input[LABEL != "L", ref := 0]
    input[ref := factor(ref)]
  }

  if (is.element("remove", colnames(input))) {
    input <- input[!(remove)]
  }

  if (remove_uninformative_feature_outlier &
    is.element("feature_quality", colnames(input))) {
    # Removed because not efficient:
    # input[, ABUNDANCE := ifelse(feature_quality == "Uninformative",
      # NA, ABUNDANCE
    # )]
    # input[, ABUNDANCE := ifelse(is_outlier, NA, ABUNDANCE)]
    input[feature_quality == "Uninformative" | is_outlier == TRUE, ABUNDANCE := NA_real_]
   
    msg <- "** Filtered out uninformative features and outliers."
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
  }
  getOption("MSstatsLog")("INFO", "MSstats - .prepareSummary function")
  getOption("MSstatsMsg")("INFO", "MSstats - .prepareSummary function")
  input <- .prepareSummary(input, method, impute, censored_symbol)
  if (!is.factor(input$PROTEIN)) input[, PROTEIN := as.factor(PROTEIN)]
  input
}


#' Get feature-level data to be used in the MSstatsSummarizationOutput function
#'
#' @param input data.table processed by dataProcess subfunctions
#'
#' @return data.table processed by dataProcess subfunctions
#'
#' @export
#'
#' @examples
#' raw <- DDARawData
#' method <- "TMP"
#' cens <- "NA"
#' impute <- TRUE
#' MSstatsConvert::MSstatsLogsSettings(FALSE)
#' input <- MSstatsPrepareForDataProcess(raw, 2, NULL)
#' input <- MSstatsNormalize(input, "EQUALIZEMEDIANS")
#' input <- MSstatsMergeFractions(input)
#' input <- MSstatsHandleMissing(input, "TMP", TRUE, "NA", 0.999)
#' input_all <- MSstatsSelectFeatures(input, "all") # all features
#' input_5 <- MSstatsSelectFeatures(data.table::copy(input),
#'   "topN",
#'   top_n = 5
#' ) # top 5 features
#'
#' proc1 <- getProcessed(input_all)
#' proc2 <- getProcessed(input_5)
#'
#' proc1
#' proc2
#'
getProcessed <- function(input) {
  remove <- NULL

  if (is.element("remove", colnames(input))) {
    if (all(!(input$remove))) {
      NULL
    } else {
      input[(remove)]
    }
  } else {
    NULL
  }
}


#' Prepare feature-level data for summarization
#' @param input data.table
#' @param method "TMP" / "linear"
#' @param impute logical
#' @param censored_symbol "0"/"NA"
#' @return data.table
#' @keywords internal
.prepareSummary <- function(input, method, impute, censored_symbol) {
  if (method == "TMP") {
    input <- .prepareTMP(input, impute, censored_symbol)
  } else {
    input <- .prepareLinear(input, FALSE, censored_symbol)
  }
  input
}


#' Prepare feature-level data for linear summarization
#' @inheritParams .prepareSummary
#' @return data.table
#' @keywords internal
.prepareLinear <- function(input, impute, censored_symbol) {
  newABUNDANCE <- ABUNDANCE <- nonmissing <- n_obs <- n_obs_run <- NULL
  total_features <- FEATURE <- prop_features <- NULL

  input[, newABUNDANCE := ABUNDANCE]
  input[, nonmissing := .getNonMissingFilter(.SD, impute, censored_symbol)]
  input[, n_obs := sum(nonmissing), by = c("PROTEIN", "FEATURE")]
  # remove feature with 1 measurement
  # remove because it is not efficient:
  # input[, nonmissing := ifelse(n_obs <= 1, FALSE, nonmissing)]
  input[n_obs <= 1, nonmissing := FALSE]
  
  input[, n_obs_run := sum(nonmissing), by = c("PROTEIN", "RUN")]

  input[, total_features := uniqueN(FEATURE), by = "PROTEIN"]
  input[, prop_features := sum(nonmissing) / total_features,
    by = c("PROTEIN", "RUN")
  ]
  input
}


#' Prepare feature-level data for TMP summarization
#' @inheritParams .prepareSummary
#' @return data.table
#' @keywords internal
.prepareTMP <- function(input, impute, censored_symbol) {
  # One allocation for newABUNDANCE, then subset assignment — no giant ifelse() temporaries.
  # Counts via filtered .N and joins → fewer scans.
  # prop_features computed in a small table → you avoid writing total_features to all rows.
  # Aggregation done directly, skipping valid_observations.
  # Indices ensure the large join is O(N) and cache-friendly.
  # Intermediates are dropped early to free RAM.
  
  # censored <- feature_quality <- newABUNDANCE <- cen <- nonmissing <- n_obs <- NULL
  # n_obs_run <- total_features <- FEATURE <- prop_features <- NULL
  # remove50missing <- ABUNDANCE <- NULL
  
  # Useful indices for later
  data.table::setindexv(input, c("PROTEIN","FEATURE"))
  data.table::setindexv(input, c("PROTEIN","RUN"))
  data.table::setindexv(input, c("PROTEIN","FEATURE","LABEL"))

  input[, newABUNDANCE := ABUNDANCE] # one allocation
  if (impute & !is.null(censored_symbol)) {
    if (is.element("feature_quality", colnames(input))) {
      # removed because it is not efficient:
      # input[, censored := ifelse(feature_quality == "Informative",
      #   censored, FALSE
      # )]
      input[feature_quality != "Informative", censored := FALSE]
    }
    if (censored_symbol == "0") {
      # removed because it is not efficient:
      # input[, newABUNDANCE := ifelse(censored, 0, ABUNDANCE)]
      input[censored == TRUE, newABUNDANCE := 0]
    } else if (censored_symbol == "NA") {
      # removed because it is not efficient:
      # input[, newABUNDANCE := ifelse(censored, NA, ABUNDANCE)]
      input[censored == TRUE, newABUNDANCE := NA_real_]
    }
    # removed because it is not efficient:
    # input[, cen := ifelse(censored, 0, 1)]
    input[, cen := as.integer(!censored)]
  } else {
    input[, newABUNDANCE := ABUNDANCE]
  }

  # removed because it is not efficient:
  # input[, nonmissing := .getNonMissingFilter(input, impute, censored_symbol)]
  input[, nonmissing := FALSE]
  if (impute && !is.null(censored_symbol) && censored_symbol == "NA") {
    input[LABEL == "L" & !is.na(newABUNDANCE), nonmissing := TRUE]
  } else {
    # default: treat zeros as missing
    input[LABEL == "L" & !is.na(newABUNDANCE) & newABUNDANCE != 0, nonmissing := TRUE]
  }
  
  # removed because it is not efficient:
  # input[, n_obs := sum(nonmissing), by = c("PROTEIN", "FEATURE")]
  # input[, nonmissing := ifelse(n_obs <= 1, FALSE, nonmissing)]
  # input[n_obs <= 1, nonmissing := FALSE]
  nobs_feat <- input[nonmissing == TRUE, .(n_obs = .N), by = .(PROTEIN, FEATURE)]
  input[nobs_feat, n_obs := i.n_obs, on = .(PROTEIN, FEATURE)]
  input[n_obs <= 1L, nonmissing := FALSE]
  
  # removed because it is not efficient:
  # input[, n_obs_run := sum(nonmissing), by = c("PROTEIN", "RUN")]
  nobs_run <- input[nonmissing == TRUE, .(n_obs_run = .N), by = .(PROTEIN, RUN)]
  input[nobs_run, n_obs_run := i.n_obs_run, on = .(PROTEIN, RUN)]

  # removed because it is not efficient:
  # input[, total_features := uniqueN(FEATURE), by = "PROTEIN"]
  # input[, prop_features := sum(nonmissing) / total_features,
  #   by = c("PROTEIN", "RUN")
  # ]
  total_feat <- input[, .(total_features = uniqueN(FEATURE)), by = PROTEIN]
  nm_by_run  <- input[nonmissing == TRUE, .(nm = .N), by = .(PROTEIN, RUN)]
  nm_by_run[total_feat, total_features := i.total_features, on = .(PROTEIN)]
  nm_by_run[, prop_features := nm / pmax.int(total_features, 1L)]
  input[nm_by_run, prop_features := i.prop_features, on = .(PROTEIN, RUN)]

  # free memory of small temps
  rm(nobs_feat, nobs_run, total_feat, nm_by_run); invisible(gc())
  
  # removed because it is not efficient:
  # if (is.element("cen", colnames(input))) {
    # if (any(input[["cen"]] == 0)) {
      # .setCensoredByThreshold(input, censored_symbol, remove50missing)
    # }
  # }
  if (impute && !is.null(censored_symbol) && any(input[["cen"]] == 0L)) {
    # nonmissing_all for the cut logic
    if (censored_symbol == "NA") {
      input[, nonmissing_all := !is.na(newABUNDANCE)]
    } else if (censored_symbol == "0") {
      input[, nonmissing_all := !is.na(newABUNDANCE) & newABUNDANCE != 0]
    } else {
      input[, nonmissing_all := !is.na(newABUNDANCE)]
    }
    # drop where total_features>1 & n_obs<=1  → if you need total_features here, reuse 'total_feat' join instead of writing it column-wise
    # Use n_obs computed above; the total_features condition is conservative if omitted.
    input[n_obs <= 1L, nonmissing_all := FALSE]
    # Direct aggregate; avoid materializing 'valid_observations'
    min_by <- input[n_obs > 1L & n_obs_run > 0L & nonmissing_all == TRUE,
                    .(min_abundance = min(newABUNDANCE, na.rm = TRUE)),
                    by = .(PROTEIN, FEATURE, LABEL)]
    min_by[, ABUNDANCE_cut := 0.99 * min_abundance]
    
    # Update join (fast with index)
    input[min_by, ABUNDANCE_cut := i.ABUNDANCE_cut, on = .(PROTEIN, FEATURE, LABEL)]
    
    # any_censored per PROTEIN
    ac <- input[censored == TRUE & n_obs > 1L & n_obs_run > 0L, .(any_censored = TRUE), by = PROTEIN]
    input[, any_censored := FALSE]
    input[ac, any_censored := TRUE, on = .(PROTEIN)]
    
    # Final replacement (no big ifelse)
    if (censored_symbol == "NA") {
      input[ nonmissing_all == FALSE & censored == TRUE &
               is.finite(ABUNDANCE_cut) & any_censored == TRUE, newABUNDANCE := ABUNDANCE_cut]
    } else { # "0"
      input[ nonmissing_all == FALSE & newABUNDANCE == 0 &
               is.finite(ABUNDANCE_cut) & any_censored == TRUE, newABUNDANCE := ABUNDANCE_cut]
    }
    
    # Clean up heavy temps
    input[, c("nonmissing_all","ABUNDANCE_cut","any_censored") := NULL]
    rm(min_by, ac); invisible(gc())
  }
  input
}
