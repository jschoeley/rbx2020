# Specify expected deaths models
#
# Most models take a data frame with training data
# which is used for model estimation, called <df_training>,
# and a data frame used to make prediction given the fitted model,
# called <df_prediction>.
#
# The time series and INLA models take a single data frame with a
# <sample_name> variable indicating training and test data.
#
# All models output a list with two elements called <model>,
# an object containing the fitted model, coefficients etc., and
# <predictions>, which is the input data frame <df_prediction> with
# an added column called <deaths_predicted>.

# Init ------------------------------------------------------------

library(here); library(glue)

wd <- here()

cnst <- list()
cnst <- within(cnst, {
  path_mod_spec = glue('{wd}/tmp/mod_spec.rds')
})

ModSpec <- list()

# Simple averages model -------------------------------------------

# This model and estimates the average mortality rate over
# some years within each week and stratum. The associated
# predict() method multiplies this average mortality with
# given exposures to derive death counts.
ModSpec$ModSimpleAverages <- function (
  df_training, df_prediction,
  week_name, death_name, exposure_name, year_name,
  n_years,
  stat = 'mortality', ...
) {
  
  require(dplyr)
  .strata <- enquos(...)
  .week <- enquo(week_name)
  .deaths <- enquo(death_name)
  .exposures <- enquo(exposure_name)
  .year <- enquo(year_name)

  # last <n_years> number of years before start of test data
  training_years <-
    df_training %>% pull(!!.year) %>% unique() %>%
    sort(decreasing = TRUE) %>% `[`(1:n_years)
    
  df_training_grouped_filtered <-
    df_training %>%
    filter(!!.year %in% training_years) %>%
    group_by(!!!.strata, !!.week)
    
  if (identical(stat, 'mortality')) {
    average <-
      df_training_grouped_filtered %>%
      summarise(avg = mean(!!.deaths / !!.exposures), .groups = 'drop')
    predictions <-
      suppressMessages(left_join(df_prediction, average)) %>%
      mutate(deaths_predicted = avg*!!.exposures) %>%
      select(-avg)
  }
  
  if (identical(stat, 'deaths')) {
    average <-
      df_training_grouped_filtered %>%
      summarise(avg = mean(!!.deaths), .groups = 'drop')
    predictions <-
      suppressMessages(left_join(df_prediction, average)) %>%
      mutate(deaths_predicted = avg) %>%
      select(-avg)
  }
  
  model <- structure(
    list(
      averages = average,
      strata_name = .strata,
      week_name = .week,
      death_name = .deaths,
      exposure_name = .exposures
    ),
    class = 'avg'
  )
  
  return(
    list(model = model, predictions = predictions)
  )
  
}

# Serfling model --------------------------------------------------

ModSpec$ModSerfling <- function (
  df_training, df_prediction,
  formula, family,
  name_week, weeks_for_training = NULL,
  name_year, n_years_for_training = NULL
  
) {

  require(dplyr)
  
  .week <- enquo(name_week)
  .year <- enquo(name_year)
  
  if (!is.null(weeks_for_training)) {
    # only train on these weeks
    df_training <-
      df_training %>% filter(!!.week %in% weeks_for_training)
  }
  if (!is.null(n_years_for_training)) {
    # last <n_years> number of years before start of test data
    training_years <-
      df_training %>% pull(!!.year) %>% unique() %>%
      sort(decreasing = TRUE) %>% `[`(1:n_years_for_training)
    # only train on these years
    df_training <-
      df_training %>% filter(!!.year %in% training_years)
  }
  
  model <- glm(
    formula = formula,
    family = family,
    data = df_training
  )
  
  predictions <-
    df_prediction %>%
    mutate(
      deaths_predicted =
        predict(model, newdata = ., type = 'response')
    )
  
  return(
    list(model = model, predictions = predictions)
  )
  
}

# GAM model -------------------------------------------------------

ModSpec$ModGAM <- function (
  df_training, df_prediction,
  formula, family, method
) {
  
  require(dplyr)
  require(mgcv)
  
  model <- gam(
    formula = formula,
    family = family,
    data = df_training,
    method = method
  )
  
  predictions <-
    df_prediction %>%
    mutate(
      deaths_predicted =
        predict(model, newdata = ., type = 'response')
    )
  
  return(
    list(model = model, predictions = predictions)
  )
  
}

# Time series models ----------------------------------------------

# very rough version of exponential smoothing model
# assumes NAs only at beginning of training series, not within
# does not test for it
ModSpec$ModTS <- function (df, sample_name, deaths_name, date_name, stratum_id_name,
                       lags, model, persistence = NULL,
                       orders = list(ar = c(0), i = c(0), ma = c(0), select = FALSE),
                       distribution = 'default', initial = 'optimal') {
  
  require(smooth); require(dplyr)
  
  .stratum <- enquo(stratum_id_name)
  .sample <- enquo(sample_name)
  .deaths <- enquo(deaths_name)
  .date <- enquo(date_name)
  
  predictions <-
    df %>%
    group_by(!!.stratum) %>%
    group_modify(~{
      
      log1p_deaths_ts_train <-
        .x %>%
        arrange(!!.date) %>%
        drop_na(!!.deaths) %>%
        filter(!!.sample == 'training') %>%
        pull(var = !!.deaths) %>%
        log1p()
      
      # length of prediction data
      n_pred <- .x %>% nrow()
      # length of test data
      n_test <- .x %>% filter(!!.sample == 'test') %>% nrow()
      # length of training data
      n_train <- n_pred - n_test
      # length of no-na time series training data
      n_train_no_na <- length(log1p_deaths_ts_train)
      # length of na time series training data
      n_train_na <- n_train - n_train_no_na
      
      
      es_fit <- adam(
        log1p_deaths_ts_train, model = model, silent = TRUE,
        h = n_test, lags = lags, persistence = persistence,
        orders = orders, initial = initial,
        distribution = distribution
      )
      
      .x %>%
        mutate(
          deaths_predicted = c(rep(NA, n_train_na),
                               expm1(es_fit$fitted),
                               expm1(es_fit$forecast))
        )
      
    }) %>%
    ungroup()
  
  return(
    list(model = NA, predictions = predictions)
  )
  
}

# PCLM model ------------------------------------------------------

#' Complete Incompletely Observed Seasonal Death Counts via PCLM
#'
#' Given a time series of weekly death counts nested within epi-years
#' where the last year is incompletely observed, make predictions for
#' weekly death counts over all the weeks of all the years.
#'
#' t_0: first week of observations
#' t_cutoff: first week without observations
#' t_end: last week which requires estimations
#' w_max: total number of weeks in a year
#' 
#' Idea by Silvia Rizzi, Implementation by Jonas Schöley and Silvia Rizzi
ModSpec$ModPCLM <- function (
  df,
  #deaths_name, epi_year_name, epi_week_name, sample_name, stratum_id_name,
  t_0 = 0, t_cutoff = 26, t_end = 51, w_max = 52,
  lambda = 1e6, diff_deg = 3, cyclic = TRUE
) {
  
  require(data.table)
  
  # names
  .D <- quote(deaths_observed)
  .y <- quote(epi_year_int)
  .w <- quote(epi_week)
  .z <- quote(stratum_id)
  .s <- quote(cv_sample)
  
  # subset data to years with complete information and separate
  # separate training data
  df_all <- as.data.table(df)
  df_complete <- df_all[,.SD[!anyNA(eval(.D))], by = list(eval(.y), eval(.z))]
  df_train <- subset(df_complete, eval(.s) == 'training')
  
  # segments of a year
  
  # partition a year into two segments at these weeks:
  # a (observed part): [t0, tcutoff),
  # b (unobserved part): [tcutoff, tend]
  w_a = t_0:(t_cutoff-1)
  w_b = t_cutoff:t_end
  
  ### Step 1: estimate deaths in different parts of the year
  
  # D_yz: total deaths in year y, stratum z
  # D_ayz: total deaths in segment a of a year y, stratum z
  D <- df_train[,.(
    D_yz = sum(eval(.D)),
    D_ayz = sum(eval(.D)[eval(.w) %in% w_a])),
    by = list(y = eval(.y), z = eval(.z))
  ]
  
  # r_ayz: ratio of total deaths to deaths in segment a, year y stratum z
  D[, r_ayz := D_yz / D_ayz, by = list(y, z)]
  
  # r_az: r_ayz averaged over years
  # this is the estimator of total deaths over an unobserved segment b
  # TODO: ignore incompletely observed years in mean calculation
  # right now I just assume that the last epi-year is not completely
  # observed and thus excluded from forming the mean
  D[, r_az := mean(r_ayz[y != max(y)]), by = z]
  
  # Dhat_yz: estimated total deaths in year y, stratum z
  D[, Dhat_yz := D_ayz*r_az, by = z]
  
  # Dhat_byz: estimated total deaths in segment b, year y, stratum z
  D[, Dhat_byz := Dhat_yz - D_ayz, by = z]
  
  ### Step 2: Estimate weekly number of deaths over periods w_a and w_b
  ###         separately by year and stratum
  
  # univariate pclm source function
  pclm <- function(y, C, X, lambda, deg, cyclic = TRUE){
    # Fit a PCLM (estimate b in ) E(y) = C %*% exp(X %*% b)
    # y = the vector of observed counts of length i
    # C = the composition matrix of dimension IxJ
    # X = the identity matrix of dimension JxJ; or B-spline basis
    # lambda = smoothing parameter
    # deg = order of differences of the components of b
    # show = indicates whether iteration details should be shown
    # cyclic = should the difference between the last and the first observation
    # be penalized?
    
    # Some preparations
    nx <- dim(X)[2]
    D <- diff(diag(nx), diff = deg)
    if (cyclic) {
      D <- rbind(
        # penalize first and second differences between last and first values
        c(1, rep(0, nx-2), -1),
        c(1, rep(0, nx-3), 1, -2),
        c(-2, 1, rep(0, nx-3), 1),
        D
      )
    }
    rowD <- nrow(D)
    bstart <- log(sum(y) / nx);
    b <- rep(bstart, nx);
    
    # Perform the iterations
    for (it in 1:50) {
      b0 <- b
      eta <- X %*% b
      gam <- exp(eta)
      mu <- C %*% gam
      w <- c(1 / mu, rep(lambda, rowD)) 
      Gam <- gam %*% rep(1, nx)
      Q <- C %*% (Gam * X)
      z <- c(y - mu + Q %*% b, rep(0, rowD))
      Fit <- lsfit(rbind(Q, D), z, wt = w, intercept = F) 
      b <- Fit$coef
      
      db <- max(abs(b - b0))
      if (db < 1e-6) break
      
    }
    
    # Regression diagnostics
    R <- t(Q) %*% diag(c(1 / mu)) %*% Q
    H <- solve(R + lambda * t(D) %*% D, R)
    H0 <- solve(R + lambda * t(D) %*% D) # variance-covariance matrix Bayesian approach
    H1 <- H0 %*% R %*% H0 # variance-covariance matrix sandwich estimator
    fit <- list()
    fit$trace <- sum(diag(H))
    ok <- y > 0
    fit$dev <- 2 * sum(y[ok] * log(y[ok] / mu[ok]))
    fit$gamma <- gam
    fit$aic <- fit$dev + 2 * fit$trace
    fit$bic <- fit$dev + log(length(y)) * fit$trace
    fit$mu <- mu
    fit$H0 <- H0
    fit$H1 <- H1
    fit$eta <- eta
    fit
  }
  
  # set up indices
  
  # weeks in each segments
  n_a = length(w_a); n_b = length(w_b)
  # position index over weeks of segments a, b
  i_a = 1:n_a; i_b = tail(i_a, 1)+1:n_b
  # length of ungrouped data
  n_fine = n_a + n_b
  # length of grouped data (number of groups)
  n_coarse = n_a + 1
  # position index over ungrouped data
  index_fine = c(i_a, i_b)
  # position index over grouped data
  index_coarse = c(i_a, rep(i_b[1], n_b))
  
  # make (C)omposition matrix
  C <- matrix(0, n_coarse, n_fine)
  C[cbind(index_coarse, index_fine)] <- 1
  
  # make (B)asis matrix
  B <- diag(n_fine)
  
  # set up data frame for results
  Y <- as.data.table(expand.grid(
    w = c(w_a, w_b),
    y = unique(df_train[,eval(.y)]),
    z = unique(df_train[,eval(.z)]),
    deaths_predicted = as.numeric(NA)
  ))
  Yiter <- Y[w == w_a[1]]
  
  # apply pclm model to each year and stratum in the input data
  for (i in 1:nrow(Yiter)) {
    # subset data
    this_year = Yiter$y[i]; this_stratum = Yiter$z[i]
    #cat(this_year, this_stratum, '\n')
    df_yz <- df_train[eval(.y) == this_year & eval(.z) == this_stratum]
    
    # prepare plcm data
    y_a <- df_yz[eval(.w) %in% w_a, .(D_w = eval(.D)),
                 by = list(eval(.w))][['D_w']]
    y_b <- D[y == this_year & z == this_stratum][['Dhat_byz']]
    y <- c(y_a, y_b)
    
    # fit pclm and predict death counts
    mod <- pclm(y, C, B, lambda = lambda, deg = diff_deg, cyclic = cyclic)
    y_hat <- mod$gamma[c(i_a, i_b)]
    Y[y == this_year & z == this_stratum ,'deaths_predicted'] <- y_hat
  }
  
  # merge back into input data
  predictions <-
    merge.data.table(
      df_all, Y,
      by.x = as.character(c(.y, .w, .z)),
      by.y = c('y', 'w', 'z'),
      all.x = TRUE
    )
  
  return(
    list(model = mod, predictions = predictions)
  )
  
}

# CODA model ------------------------------------------------------

# Idea by James Oeppen; Implementation and smoothing extension by
# Jonas Schöley
ModSpec$ModCODA <- function (
  df_input,
  seasonality_indicator, formula, transform,
  winter_death_weeks, exogenous_weeks,
  gamma = 1, case_weights = 'none', method = 'GCV.Cp',
  stratum_name, week_fac_name, death_name, exposure_name, sample_name, year_name
) {
  
  require(dplyr)
  require(compositions)
  require(mgcv)
  .stratum_id <- enquo(stratum_name)
  .week <- enquo(week_fac_name)
  .deaths <- enquo(death_name)
  .exposures <- enquo(exposure_name)
  .sample <- enquo(sample_name)
  .year <- enquo(year_name)
  
  # because coda can't train on partial information for a year
  # add an indicator stating if a given year contains partial data
  dat_input <-
    df_input %>%
    group_by(!!.stratum_id, !!.year) %>%
    mutate(
      year_contains_partial_data =
        any(!!.sample == 'test') | any(is.na(!!.deaths))
    ) %>%
    ungroup() %>%
    # add row id for later identification
    mutate(row_id = 1:n())
  
  # year and stratum specific summary variables
  dat_yj <-
    dat_input %>%
    # add this because quosure can't be used in log2()
    mutate(.weekly_deaths = !!.deaths) %>%
    group_by(!!.stratum_id, !!.year) %>%
    summarise(
      # for deriving total deaths later on and for seasonality indicators
      deaths_exogenous_yj =
        sum(ifelse(!!.week %in% exogenous_weeks, .weekly_deaths, 0)),
      n_weeks_exogeneous_yj =
        length(exogenous_weeks),
      deaths_exogenous_entropy_yj =
        sum(ifelse(!!.week %in% exogenous_weeks & .weekly_deaths != 0,
                   .weekly_deaths*log2(.weekly_deaths), 0)),
      annual_deaths_yj =
        # annual deaths are unknown for years where only partial
        # training data is available
        sum(ifelse(year_contains_partial_data, NA, .weekly_deaths)),
      deaths_winter_yj =
        sum(ifelse(!!.week %in% winter_death_weeks, .weekly_deaths, 0)),
      exposure_winter_yj =
        sum(ifelse(!!.week %in% winter_death_weeks, !!.exposures, 0)),
      # seasonality indicators
      seasonality_yj =
        switch (seasonality_indicator,
                'deaths' = deaths_winter_yj,
                'mortality' = deaths_winter_yj / exposure_winter_yj,
                'ratio' = deaths_winter_yj/deaths_exogenous_yj,
                'entropy' = (log2(deaths_exogenous_yj)-
                               1/deaths_exogenous_yj*deaths_exogenous_entropy_yj)/
                  log2(n_weeks_exogeneous_yj)
        ),
      .groups = 'drop'
    )
  
  # stratum specific summary variables
  dat_j <-
    dat_yj %>%
    group_by(!!.stratum_id) %>%
    # seasonality indicator average and sd over years by stratum
    summarise(
      avg_seasonality_j =
        mean(seasonality_yj),
      sd_seasonality_j =
        sd(seasonality_yj),
      .groups = 'drop'
    )
  suppressMessages({
    # add multilevel predictors to input data and compute
    # share on annual death by week per stratum
    dat_jyw <-
      dat_input %>%
      left_join(
        dat_yj,
        by = c(as_label(.stratum_id), as_label(.year))
      ) %>%
      left_join(
        dat_j,
        by = c(as_label(.stratum_id))
      ) %>%
      mutate(
        # predictors
        z_seasonality_yj =
          (seasonality_yj - avg_seasonality_j) / sd_seasonality_j,
        # outcome variable
        share_on_annual_deaths_wyj =
          !!.deaths / annual_deaths_yj
      ) %>%
      ungroup()
    
    # apply either ILR or ALR transform to outcome variable
    dat_jyw <-
      dat_jyw %>%
      arrange(!!.stratum_id, !!.year, !!.week) %>%
      group_by(!!.stratum_id, !!.year) %>%
      group_modify(.f = ~{
        x <- .x[['share_on_annual_deaths_wyj']]
        # number of weeks in composition
        # may be 52 or 53
        n_weeks <- length(pull(.x, !!.week))
        # no 0s allowed in CODA...
        fudge <- 1e-6
        if (any(is.na(x))) {
          # if part of the annual share is missing
          # don't bother with a CODA transform
          z <- NA
        } else {
          z <- switch(
            transform,
            ilr = ilr(x+fudge, V = ilrBase(D=n_weeks)),
            alr = alr(x+fudge, ivar = n_weeks) 
          )
          z <- c(z, NA)
        }
        mutate(
          .x, transformed_share_on_annual_deaths = z,
          n_weeks = n_weeks
        )
      }) %>%
      ungroup()
    
    # add interaction effects between seasonal predictor and stratum
    # this allows to specify a stratified "varying coefficients model"
    # via mgcv::gam() where the effect of z_seasonality_yj may vary
    # over time in a smooth fashion
    dat_mod_mat <-
      model.matrix(
        as.formula(paste0('~-1+z_seasonality_yj:', as_label(.stratum_id))),
        data = dat_jyw)
    stratum_names <-
      unique(pull(dat_input, !!.stratum_id))
    dimnames(dat_mod_mat)[[2]] <-
      paste0('z_seasonality_yj',1:length(stratum_names))
    dat_jyw <-
      cbind(dat_jyw, dat_mod_mat)
    
    # remove reference week
    dat_jyw_no_reference_week <-
      dat_jyw %>%
      # within each stratum-year remove the reference week
      group_by(!!.stratum_id, !!.year) %>%
      slice(-n_weeks[1]) %>%
      mutate(across(!!.week, factor)) %>%
      ungroup()
    # prepare training and prediction data
    dat_train <-
      dat_jyw_no_reference_week %>%
      filter(year_contains_partial_data == FALSE)
    dat_pred <- dat_jyw_no_reference_week
    
    # case weights
    if (identical(case_weights, 'none')) {
      dat_train <- mutate(dat_train, case_weight = 1)
    } else {
      dat_train <-
        mutate(dat_train, case_weight = log1p(!!case_weights))
    }
    
    # fit model
    model <-
      gam(
        formula = formula,
        family = gaussian(link = 'identity'),
        data = dat_train,
        weights = case_weight,
        gamma = gamma,
        method = method
      )
    
    # predict from fitted model
    predictions <-
      dat_pred %>%
      mutate(
        # predicted transformed weekly death proportions by
        # epi-year and stratum
        y_hat =
          predict(model, newdata = ., type = 'response')
      ) %>%
      select(row_id, y_hat) %>%
      # merge into original data which does include reference week
      right_join(dat_jyw, by = 'row_id') %>%
      # within each stratum-year, convert predicted transformed
      # weekly death proportions to proportion scale
      group_by(!!.stratum_id, !!.year) %>%
      group_modify(~{
        # predicted transformed share without reference week
        z_hat <- .x$y_hat[-(.x$n_weeks[1])]
        # back-transform from coda space to proportion space
        x_hat <- switch(
          transform,
          ilr = ilrInv(z_hat, V = ilrBase(D=.x$n_weeks[1])),
          alr = alrInv(z_hat)
        )
        mutate(.x, pred_share_on_annual_deaths_yjw = c(x_hat))
      }) %>%
      mutate(
        pred_share_exogenous_deaths_yj =
          sum(ifelse(!!.week %in% exogenous_weeks, pred_share_on_annual_deaths_yjw, 0))
      ) %>%
      select(-n_weeks) %>%
      ungroup() %>%
      # convert to death count scale
      mutate(
        pred_annual_deaths_yj =
          deaths_exogenous_yj/pred_share_exogenous_deaths_yj,
        deaths_predicted =
          pred_share_on_annual_deaths_yjw*pred_annual_deaths_yj
      ) %>%
      select(row_id, deaths_predicted) %>%
      right_join(dat_input, by = 'row_id') %>%
      select(-row_id, -year_contains_partial_data)
  })
  
  return(
    list(model = model, predictions = predictions)
  )
  
}

# INLA ------------------------------------------------------------

# Re-implementation of https://doi.org/10.1038/s41591-020-1112-0
ModSpec$ModINLA <- function (
  df_traintest,
  inla_formula,
  stratum_name,
  sample_name,
  death_name,
  exposure_name,
  tanomaly_name,
  week_name,
  time_name,
  holiday_name,
  threads
) {
  
  require(dplyr)
  require(INLA)
  
  # translate external to internal variables names
  .stratum <- enquo(stratum_name)
  .sample <- enquo(sample_name)
  .death <- enquo(death_name)
  .exposure <- enquo(exposure_name)
  .tanomaly <- enquo(tanomaly_name)
  .week <- enquo(week_name)
  .time <- enquo(time_name)
  .holiday <- enquo(holiday_name)
  
  # prepare data for fit
  ready_for_fit <-
    df_traintest %>%
    # select variables of interest
    select(
      sample = !!.sample,
      stratum = !!.stratum,
      time = !!.time,
      week = !!.week,
      death = !!.death,
      exposure = !!.exposure,
      tanomaly = !!.tanomaly,
      holiday = !!.holiday
    ) %>%
    # if the observation is part of the test-set
    # make the outcome <NA>
    mutate(
      death = ifelse(sample == 'test', NA, death)
    ) %>%
    # add variables needed in model fit
    mutate(
      global_slope = time,
      time_ar = time,
      time_seas = time,
      week_rw = week,
      resid_iid = time,
      logexposure = log(exposure)
    )
  
  model <-
    ready_for_fit %>%
    group_by(stratum) %>%
    group_modify(~{
      
      the_fit <- inla(
        formula = inla_formula,
        offset = logexposure,
        family = 'poisson',
        control.predictor = list(link = 1),
        # we are only interested in predictions, don't
        # compute anything else
        control.compute = list(
          dic = FALSE, config = FALSE, hyperpar = FALSE,
          return.marginals = FALSE, mlik = FALSE
        ),
        num.threads = threads,
        data = .x
      )
      
      mutate(
        .x,
        deaths_predicted = the_fit$summary.fitted.values$mean
      )
      
    })
  
  predictions <-
    df_traintest %>% mutate(
      deaths_predicted = model$deaths_predicted
    )
  
  return(
    list(model = NA, predictions = predictions)
  )
  
}

# Export ----------------------------------------------------------

saveRDS(ModSpec, file = cnst$path_mod_spec)
