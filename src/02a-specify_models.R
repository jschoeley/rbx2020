# Specify expected deaths models
#
# We specify models which, given some training data, predict
# death counts over training and test data. The models all
# expect a data frame <df> as input with a column indicating
# for each row if it is part of the 'training' or the 'test'
# data. The <df> is supposed to feature weekly observations on
# death counts by population strata.
#
# All models output the original input <df> with an added column called
# <deaths_predicted> holding the expected death count and added columns
# <deaths_sim><1:nsim> holding simulated death counts from
# the predicted distribution of counts.

# Init ------------------------------------------------------------

library(here); library(glue)

wd <- here()

cnst <- list()
cnst <- within(cnst, {
  path_mod_spec = glue('{wd}/tmp/mod_spec.rds')
})

ModSpec <- list()

# SerflingGLM -----------------------------------------------------

#' Serfling GLM
#'
#' @param df data frame containing the variables in the model.
#' @param models list of formulas passed to glm().
#' @param family family object passed to glm(). poisson() or quasipoisson().
#' @param col_sample name of column in <df> indicating training or test
#' data. must have values 'training' or 'test'.
#' @param col_stratum name of column in <df> indicating strata.
#' @param weeks_for_training vector of weeks in training data to be used
#' for training the model. default NULL uses all weeks of the year.
#' @param col_week name of column used for <weeks_for_training> selection
#' @param n_years_for_training number of years in training data to be used
#' for training. counts backwards from the last year in training data.
#' default NULL uses all years in training.
#' @param col_year name of column used for <n_years_for_training> selection.
#' @param nsim number of simulated predictions.
#' @param simulate_beta should the simulated predictions contain
#' uncertainty around the parameter estimates of the model? (default = TRUE)
#' @param simulate_y should the simulated predictions contain uncertainty
#' around the sampling distribution of the outcome (default = TRUE)
#'
#' @details
#' A Serfling style GLM is fitted over the training data and expected
#' death counts are predicted over the complete input data frame. The
#' training data is indicated by the column <col_sample> and can further
#' be subset by specifying <weeks_for_training> and <n_years_for_training>.
#' All models in <models> are fit separately by stratum and within
#' each stratum the model with the lowest AIC is chosen for prediction.
#' By default, the input <df> is returned with added expected death
#' counts and <nsim> columns holding simulated death counts from the
#' predicted distribution of counts.
#'
#' @return
#' <df> with added column <deaths_predicted> containing the expected
#' death counts, and columns <deaths_sim><1:nsim> containing simulated
#' expected death counts if simulate_y = FALSE or simulated deaths from
#' the predicted outcome distribution if simulate_y = TRUE.
ModSpec$SerflingGLM <- function (
  df, models, family,
  # column names for training/test split and strata
  col_sample, col_stratum,
  # only fit on part of the year
  weeks_for_training = NULL, col_week = NULL,
  # only fit on part of the available years
  n_years_for_training = NULL, col_year = NULL,
  # simulation parameters
  nsim = 100, simulate_beta = TRUE, simulate_y = TRUE
) {
  
  df['.rowid'] <- 1:nrow(df)
  
  ## subset input data to rows used for fitting the model ##
  
  # index of rows designated as training data
  idx_train <- df[[col_sample]] == 'training'
  # index of rows with weeks suitable for training
  idx_weeks <- TRUE
  if (!is.null(weeks_for_training)) {
    # only train on these weeks
    idx_weeks <- df[[col_week]] %in% weeks_for_training
  }
  # index of rows with years suitable for training
  idx_years <- TRUE
  if (!is.null(n_years_for_training)) {
    # most recent <n_years> in training data
    years_for_training <- sort(unique(df[idx_train,][[col_year]]),
                               decreasing = TRUE)[1:n_years_for_training]
    # only train on these years
    idx_years <- df[[col_year]] %in% years_for_training
  }
  # index of data used for fitting
  idx_fit <- idx_train & idx_years & idx_weeks
  
  # for each stratum, fit model, predict and simulate from model,
  # add results to df
  strata <- unique(df[[col_stratum]])
  for (i in strata) {
    
    ## fit model ##
    
    # stratum subsets of training and prediction data
    df_prd <- df[df[[col_stratum]]==i,]
    df_trn <- df[df[[col_stratum]]==i&idx_fit,]
    
    if (family$family != 'quasipoisson') {
      # fit all the models and extract the one with lowest AIC
      fits <- lapply(models, function (x) {
        glm(formula = x, family = family,
            data = df_trn, method = glm2::glm.fit2)
      })
      k <- which.min(sapply(fits, AIC))
      formula <- models[[k]]
      model <- fits[[k]]
    }
    
    if (family$family == 'quasipoisson') {
      # first fit as poisson to get AIC
      fits <- lapply(models, function (x) {
        glm(formula = x, family = poisson(),
            data = df_trn, method = glm2::glm.fit2)
      })
      # determine model with lowest AIC and refit as quasipoisson
      k <- which.min(sapply(fits, AIC))
      formula <- models[[k]]
      model <- glm(formula, family = family,
                   data = df_trn, method = glm2::glm.fit2)
    }
    
    ## predict from model ##
    
    # create a design matrix for prediction
    # keep input NAs in design matrix
    # matrix has same number of rows as <df>
    formula_rhs <- update(formula, NULL~.)
    XX_prd <- model.frame(formula_rhs, df_prd, na.action = na.pass)
    X_prd <- model.matrix(formula_rhs, XX_prd)
    # ensure that the prediction matrix only includes factor level
    # interactions as seen during fit
    # this avoids errors when a factor has a new level in the prediction
    # data set and there's no coefficient estimated for that model;
    # prediction for novel levels will be the same as prediction for
    # reference level
    admissible_terms <- colnames(model.matrix(model))
    X_prd <- X_prd[,admissible_terms]
    
    # offset will be added to linear predictor, add 0 if no offset
    x_offset <- model.offset(XX_prd)
    if (is.null(x_offset)) x_offset <- 0
    
    # expected death counts
    ILink <- model$family$linkinv
    Ey <- ILink(X_prd %*% coef(model) + x_offset)
    
    ## simulate model predictions ##
    
    # simulated model coefficients
    if (isTRUE(simulate_beta)) {
      beta_sim <- MASS::mvrnorm(nsim, coef(model), vcov(model))
    } else {
      beta_sim <- matrix(rep(coef(model), nsim), nrow = nsim, byrow = TRUE) 
    }
    # simulated expectations of the outcome distribution
    Ey_sim <- apply(beta_sim, 1, FUN = function (b) ILink(X_prd%*%b + x_offset))
    # simulated outcome
    y_sim <- apply(Ey_sim, 2, FUN = function (Ey) {
      # if a simulation algorithm hasn't been defined for a family
      # just return the expectation of the outcome
      y <- mu <- Ey
      # NA's can't be passed to the simulation functions, so keep them out
      idx_na <- is.na(mu); mu_ <- mu[!idx_na]; N <- length(mu_)
      if (model$family$family == 'poisson') {
        y[!idx_na] <- rpois(n = N, lambda = mu_)      
      }
      if (model$family$family == 'quasipoisson') {
        # https://stats.stackexchange.com/q/157575
        # we estimate the rate and dispersion parameter via quasipoisson
        # and then sample from a Negative Binomial distribution with the
        # same rate and dispersion (NBI)
        phi <- summary(model)$dispersion
        # in case of under-dispersion, sample from Poisson
        if (phi < 1) { phi = 1 }
        y[!idx_na] <- rnbinom(n = N, mu = mu_, size = mu_/(phi-1))
      }
      # just return the expectation if outcome simulation is FALSE
      if (!isTRUE(simulate_y)) {
        y <- Ey
      }
      return(y)
    })
    colnames_y_sim <- paste0('deaths_sim', 1:nsim)
    
    # add predictions and simulations to input data
    df[df_prd[['.rowid']], 'deaths_predicted'] <- Ey
    df[df_prd[['.rowid']], colnames_y_sim] <- y_sim
    
  }
  
  df[,'.rowid'] <- NULL
  return(df)
  
}

# CountGAM --------------------------------------------------------

#' Count Prediction with GAM
#'
#' @param df data frame containing the variables in the model.
#' @param formula formula for gam().
#' @param family family object passed to gam(). poisson() or quasipoisson().
#' @param col_sample name of column in <df> indicating training or test
#' data. must have values 'training' or 'test'.
#' @param nsim number of simulated predictions.
#' @param simulate_beta should the simulated predictions contain
#' uncertainty around the parameter estimates of the model? (default = TRUE)
#' @param simulate_y should the simulated predictions contain uncertainty
#' around the sampling distribution of the outcome (default = TRUE)
#'
#' @details
#' A GAM is fitted over the training data and expected
#' death counts are predicted over the complete input data frame. The
#' training data is indicated by the column <col_sample>.
#' By default, the input <df> is returned with added expected death
#' counts and <nsim> columns holding simulated death counts from the
#' predicted distribution of counts.
#'
#' @return
#' <df> with added column <deaths_predicted> containing the expected
#' death counts, and columns <deaths_sim><1:nsim> containing simulated
#' expected death counts if simulate_y = FALSE or simulated deaths from
#' the predicted outcome distribution if simulate_y = TRUE.
ModSpec$CountGAM <- function (
  df, formula, family,
  col_sample,
  nsim = 100, simulate_beta = TRUE, simulate_y = TRUE
) {
  
  require(mgcv)
  
  # training data
  df_trn <- df[df[[col_sample]]=='training',]
  
  ## fit model ##
  
  model <- gam(
    formula = formula,
    family = family,
    data = df_trn
  )
  
  ## predict from model ##
  
  # create a design matrix for prediction
  # keep input NAs in design matrix
  # matrix has same number of rows as <df>
  X_prd <- predict(model, newdata = df, type = 'lpmatrix')
  # estimated coefficients
  beta <- coef(model)
  # linear predictor over prediction data w/o offset
  eta_prd_without_offset <- X_prd %*% beta
  # linear predictor over prediction data with offset included
  eta_prd_with_offset <- matrix(predict(model, newdata = df, type = 'link'), ncol = 1)
  # I know of no easy way with mgcv to extract the offset over "newdata"
  # therefore this rather strange solution
  # offset over prediction data (may be 0)
  offset_prd <- eta_prd_with_offset - eta_prd_without_offset
  # inverse link function
  ILink <- model$family$linkinv
  # expected death counts
  Ey_prd <- ILink(eta_prd_with_offset)
  
  ## simulate model predictions ##
  
  # simulated model coefficients
  if (isTRUE(simulate_beta)) {
    beta_sim <- MASS::mvrnorm(nsim, beta, vcov(model, freq = FALSE, unconditional = TRUE))
  } else {
    beta_sim <- matrix(rep(beta, nsim), nrow = nsim, byrow = TRUE)
  }
  # simulated expectations of the outcome distribution
  Ey_sim <- apply(beta_sim, 1, FUN = function (b) ILink(X_prd%*%b + offset_prd))
  # simulated outcome
  y_sim <- apply(Ey_sim, 2, FUN = function (Ey) {
    # if a simulation algorithm hasn't been defined for a family
    # just return the expectation of the outcome
    y <- mu <- Ey
    # NA's can't be passed to the simulation functions, so keep them out
    idx_na <- is.na(mu); mu_ <- mu[!idx_na]; N <- length(mu_)
    if (model$family$family == 'poisson') {
      y[!idx_na] <- rpois(n = N, lambda = mu_)      
    }
    if (model$family$family == 'quasipoisson') {
      # https://stats.stackexchange.com/q/157575
      # we estimate the rate and dispersion parameter via quasipoisson
      # and then sample from a Negative Binomial distribution with the
      # same rate and dispersion (NBI)
      phi <- summary(model)$dispersion
      # in case of under-dispersion, sample from Poisson
      if (phi < 1) { phi = 1 }
      y[!idx_na] <- rnbinom(n = N, mu = mu_, size = mu_/(phi-1))      
    }
    # just return the expectation if outcome simulation is FALSE
    if (!isTRUE(simulate_y)) {
      y <- Ey
    }
    return(y)
  })
  colnames_y_sim <- paste0('deaths_sim', 1:nsim)
  
  # add predictions and simulations to input data
  
  df[,'deaths_predicted'] <- Ey_prd
  df[,colnames_y_sim] <- y_sim
  
  return(df)
  
}

# KontisLGM -------------------------------------------------------

#' Kontis etal. Count Prediction LGM
#'
#' @param df data frame containing the variables in the model.
#' @param formula1 formula for inla(). terms must only use the
#' following variable names:
#' death, holiday, tanomaly, exposure, global_slope, time_ar, time_seas,
#' week_rw, resid_iid.
#' @param formula2 alternative formula for inla().
#' @param formula2_threshold integer threshold. if the average weekly
#' death counts in the training data within a stratum fall below this
#' number, then use the alternative formula.
#' @param col_stratum name of column in <df> indicating strata.
#' @param col_sample name of column in <df> indicating training or test
#' data. must have values 'training' or 'test'.
#' @param col_death name of column in <df> indicating death counts.
#' @param col_exposure name of column in <df> indicating exposures.
#' @param col_tanomaly name of column in <df> indicating temperature
#' anomalies.
#' @param col_week name of column in <df> indicating week of year.
#' @param col_time name of column in <df> indicating weeks since start
#' of observation.
#' @param col_holiday name of column in <df> indicating public holidays.
#' @param nsim number of simulated predictions.
#' @param weeks_for_training_within_year vector of weeks within a year
#' in training data to be used for training the model.
#' default NULL uses all weeks of the year.
#' @param weeks_for_training_pre_test number of weeks in training data
#' to be used for training. counts backwards from the most recent week
#' in training data. default NULL uses all weeks in training.
#' @param threads number of threads to use.
#'
#' @details
#' Based upon the model in
#' https://doi.org/10.1038/s41591-020-1112-0.
#'
#' @return
#' <df> with added column <deaths_predicted> containing the expected
#' death counts, and columns <deaths_sim><1:nsim> containing simulated
#' death counts from the posterior predictive distribution.
ModSpec$KontisLGM <- function (
  df,
  formula1,
  formula2 = NULL,
  formula2_threshold = NULL,
  # variable names
  col_stratum,
  col_sample,
  col_death,
  col_exposure,
  col_tanomaly,
  col_week,
  col_time,
  col_holiday,
  # number of draws from posterior predictive distribution
  nsim = 100,
  # specify weeks used for training the model
  weeks_for_training_within_year = NULL,#c(15:26, 36:45)
  weeks_for_training_pre_test = NULL,#52*5
  threads = 1
) {
  
  require(dplyr)
  require(INLA)
  
  # translate external to internal variables names
  .stratum <- enquo(col_stratum)
  .sample <- enquo(col_sample)
  .death <- enquo(col_death)
  .exposure <- enquo(col_exposure)
  .tanomaly <- enquo(col_tanomaly)
  .week <- enquo(col_week)
  .time <- enquo(col_time)
  .holiday <- enquo(col_holiday)
  
  # add temporary row id so that we can merge the predictions
  # back into the input data
  df$temporary_row_id <- 1:nrow(df)
  
  # prepare data for fit
  ready_for_fit <-
    df %>%
    # select variables of interest
    select(
      temporary_row_id,
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
      resid_iid = time
    )
  
  # exclude weeks within a year from training
  if (!is.null(weeks_for_training_within_year)) {
    ready_for_fit <-
      ready_for_fit %>%
      mutate(death = ifelse(week %in% weeks_for_training_within_year, death, NA))
  }
  
  # only use a certain number of weeks before start of test to train model
  if (!is.null(weeks_for_training_pre_test)) {
    training_weeks <-
      ready_for_fit %>% filter(sample=='training') %>%
      pull(time) %>% unique() %>%
      sort(decreasing = TRUE) %>% `[`(1:weeks_for_training_pre_test)
    ready_for_fit <-
      ready_for_fit %>%
      mutate(death = ifelse(time %in% training_weeks, death, NA))
  }
  
  SamplePosteriorPredictiveDistribution <- function(inla_fit, n) {
    # draws from posterior distribution of parameters
    draws <- inla.posterior.sample(n, inla_fit)
    # length of data series
    data_length <- attr(draws, '.contents')$length[1]
    
    count_samples <- lapply(draws, function(draw) {
      # posterior sample of lambdas
      lambda <- exp(draw$latent[grep('Predictor', rownames(draw$latent))])
      # posterior sample of Poisson counts
      rpois(data_length, lambda)
    })
    return(count_samples)
  }
  
  model <-
    ready_for_fit %>%
    # fit the model separately by stratum
    group_by(stratum) %>%
    group_modify(~{
      
      # decide which formula to use, main or alternative
      # if the average observed deaths are below a threshold, use
      # the alternative formula
      stratum_specific_formula <- formula1
      if (isTRUE(mean(.x[['death']], na.rm = TRUE) < formula2_threshold)) {
        stratum_specific_formula <- formula2
      }
      
      # fit model
      
      # initial fit to get starting values for full fit
      # aides in convergence
      init_fit <- inla(
        formula = stratum_specific_formula,
        family = 'poisson',
        control.predictor = list(link = 1),
        control.compute = list(dic = FALSE, config = TRUE),
        control.inla = list(
          int.strategy = 'eb', strategy = 'gaussian',
          diagonal = 10000
        ),
        num.threads = threads,
        data = .x
      )
      
      # full fit
      the_fit <- inla(
        formula = stratum_specific_formula,
        family = 'poisson',
        control.predictor = list(link = 1),
        control.compute = list(
          dic = FALSE, config = TRUE),
        control.inla = list(diagonal = 0),
        control.mode = list(result = init_fit, restart = TRUE),
        num.threads = threads,
        data = .x
      )
      
      # add mean prediction
      deaths_predicted <- the_fit$summary.fitted.values$mean
      # add N prediction series sampled from posterior predictive distribution
      deaths_sampled <-
        do.call('cbind',
                SamplePosteriorPredictiveDistribution(the_fit, n = nsim))
      colnames(deaths_sampled) <- paste0('deaths_sim', 1:nsim)
      
      data.frame(
        temporary_row_id = .x$temporary_row_id,
        deaths_predicted = deaths_predicted,
        deaths_sampled
      )
      
    }) %>%
    ungroup() %>%
    select(-stratum)
  
  df <-
    left_join(df, model, by = 'temporary_row_id') %>%
    select(-temporary_row_id)
  
  return(df)
  
}

# Export ----------------------------------------------------------

saveRDS(ModSpec, file = cnst$path_mod_spec)
