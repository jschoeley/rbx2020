# Parameterize expected death models

# Init ------------------------------------------------------------

library(dplyr)
library(here); library(glue)

wd <- here()
setwd(wd)

cnst <- list()
cnst <- within(cnst, {
  path_mod_para = glue('tmp/mod_para.rds')
})

# Specifications of models to test --------------------------------

# Just a big list specifying all the models to be tested,
# and their parametrizations. See specify_models.R for the exact
# implementation of the models.

mod_para <-
  tribble(
    
    ~model_id, ~model_class, ~model_spec,
    
    # ... Average -----------------------------------------------------
    
    'AVGc5', 'avg', list(
      stat = 'deaths', n_years = 5
    ),
    'AVGc3', 'avg', list(
      stat = 'deaths', n_years = 3
    ),
    
    'AVGr5', 'avg', list(
      stat = 'mortality', n_years = 5
    ),
    'AVGr3', 'avg', list(
      stat = 'mortality', n_years = 3
    ),
    
    # ... Serfling ----------------------------------------------------
    
    # https://github.com/EuroMOMOnetwork/MOMO/blob/master/R/excess.R
    'SRFcem', 'serfling', list(
      formula = formula(
        deaths_observed ~
          # log linear long term trend
          origin_weeks*stratum_id +
          # seasonality
          # full year period
          sin(2*pi*epi_week/52)*stratum_id +
          cos(2*pi*epi_week/52)*stratum_id
      ),
      family = quasipoisson(link = 'log'),
      weeks_for_training = c(15:26, 36:45),
      n_years_for_training = 5
    ),
    'SRFc', 'serfling', list(
      formula = formula(
        deaths_observed ~
          # log linear long term trend
          origin_weeks*stratum_id +
          # seasonality
          # full year period
          sin(2*pi*epi_week/52)*stratum_id +
          cos(2*pi*epi_week/52)*stratum_id +
          # half year period
          sin(2*pi*epi_week/26)*stratum_id +
          cos(2*pi*epi_week/26)*stratum_id +
          # adjustment for special weeks
          holiday3*stratum_id
      ),
      family = quasipoisson(link = 'log'),
      weeks_for_training = NULL,
      n_years_for_training = NULL
    ),
    'SRFr', 'serfling', list(
      formula = formula(
        deaths_observed ~
          # log linear long term trend
          origin_weeks*stratum_id +
          # seasonality
          # full year period
          sin(2*pi*epi_week/52)*stratum_id +
          cos(2*pi*epi_week/52)*stratum_id +
          # half year period
          sin(2*pi*epi_week/26)*stratum_id +
          cos(2*pi*epi_week/26)*stratum_id +
          # adjustment for special weeks
          holiday3*stratum_id +
          # exposures
          offset(log(personweeks))
      ),
      family = quasipoisson(link = 'log'),
      weeks_for_training = NULL,
      n_years_for_training = NULL
    ),

    # ... GAM ---------------------------------------------------------
    
    'GAMr', 'gam', list(
      formula = formula(
        deaths_observed ~
          # log linear long term trend
          origin_weeks*stratum_id +
          # penalized cyclic spline for seasonality
          s(epi_week, bs = 'cp', k = 12, by = stratum_id) +
          # adjustment for special weeks
          holiday3*stratum_id +
          # exposures
          offset(log(personweeks))
      ),
      family = quasipoisson(link = 'log'),
      method = 'GCV.Cp'
    ),
    'GAMrt', 'gam', list(
      formula = formula(
        deaths_observed ~
          # log linear long term trend
          origin_weeks*stratum_id +
          # penalized cyclic spline for seasonality
          s(epi_week, bs = 'cp', k = 12, by = stratum_id) +
          # temperature effect
          s(epi_week, bs = 'cp', k = 12, by = temperature_anomaly) +
          # adjustment for special weeks
          holiday3*stratum_id +
          # exposures
          offset(log(personweeks))
      ),
      family = quasipoisson(link = 'log'),
      method = 'GCV.Cp'
    ),
    
    # ... Time series -------------------------------------------------
    
    # 'Seasonal exponential smoothing (AAA, 1,1,52) deaths', 'es-aaa', 'ts', list(
    #   model = 'AAA', lags = c(1,1,52),
    #   persistence = NULL,
    #   distribution = 'dnorm',
    #   initial = 'optimal',
    #   orders = list(ar = c(0), i = c(0), ma = c(0), select = FALSE)
    # ),
    
    # 'Seasonal Arima (1,1,1)(1,1,1)52 deaths', 'sarima', 'ts', list(
    #   model = 'NNN', lags = c(1,52),
    #   persistence = NULL,
    #   distribution = 'dnorm',
    #   initial = 'backcasting',
    #   orders = list(ar=c(1,1),i=c(1,1),ma=c(1,1))
    # ),
    
    # ... INLA ------------------------------------------------------
    
    'LGMr', 'inla', list(
      inla_formula = formula(
        death ~
          1 +
          global_slope +
          holiday +
          f(time_ar,
            model = 'ar1',
            hyper = list(prec = list(prior = 'loggamma', param = c(0.001, 0.001)))
          ) +
          f(time_seas,
            model = 'seasonal', season.length = 52,
            hyper = list(prec = list(prior = 'loggamma', param = c(0.001, 0.001)))
          ) +
          # independent remaining errors
          f(resid_iid,
            model = 'iid',
            hyper = list(prec = list(prior = 'loggamma', param = c(0.001, 0.001)))
          )
      )
    ),
    'LGMrt', 'inla', list(
      inla_formula = formula(
        death ~
          1 +
          global_slope +
          holiday +
          tanomaly +
          f(time_ar,
            model = 'ar1',
            hyper = list(prec = list(prior = 'loggamma', param = c(0.001, 0.001)))
          ) +
          f(time_seas,
            model = 'seasonal', season.length = 52,
            hyper = list(prec = list(prior = 'loggamma', param = c(0.001, 0.001)))
          ) +
          # effect of temperature anomaly varies in a cyclic fashion
          # over the week of a year
          f(week_rw, tanomaly,
            model = 'rw1', cyclic = TRUE,
            hyper = list(prec = list(prior = 'loggamma', param = c(0.001, 0.001)))
          ) +
          # independent remaining errors
          f(resid_iid,
            model = 'iid',
            hyper = list(prec = list(prior = 'loggamma', param = c(0.001, 0.001)))
          )
      )
    ),
    'LGMrt2', 'inla', list(
      inla_formula = formula(
        death ~
          1 +
          global_slope +
          holiday +
          tanomaly +
          f(time_ar,
            model = 'ar', order = 2,
            hyper = list(prec = list(prior = 'loggamma', param = c(0.001, 0.001)))
          ) +
          f(time_seas,
            model = 'seasonal', season.length = 52,
            hyper = list(prec = list(prior = 'loggamma', param = c(0.001, 0.001)))
          ) +
          # effect of temperature anomaly varies in a cyclic fashion
          # over the week of a year
          f(week_rw, tanomaly,
            model = 'rw2', cyclic = TRUE,
            hyper = list(prec = list(prior = 'loggamma', param = c(0.001, 0.001)))
          ) +
          # independent remaining errors
          f(resid_iid,
            model = 'iid',
            hyper = list(prec = list(prior = 'loggamma', param = c(0.001, 0.001)))
          )
      )
    )#,
    
  )

# Export ----------------------------------------------------------

saveRDS(mod_para, file = cnst$path_mod_para)
