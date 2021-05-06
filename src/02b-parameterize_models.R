# Parameterize expected death models

# A list specifying all the models to be tested, and their
# parametrizations. See specify_models.R for the exact
# implementation of the models.

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

mod_para <-
  tribble(
    
    ~model_id, ~model_class, ~model_para,
    
    # ... AVG ---------------------------------------------------------
    
    'AVRc5', 'glm', list(
      models = list(formula(
        deaths_observed ~
          # single coefficient for every weeks
          as.factor(iso_week)
      )),
      family = quasipoisson(link = 'log'),
      n_years_for_training = 5,
      weeks_for_training = NULL
    ),
    
    'AVRr5', 'glm', list(
      models = list(formula(
        deaths_observed ~
          # single coefficient for every weeks
          as.factor(iso_week) +
          offset(log(personweeks))
      )),
      family = quasipoisson(link = 'log'),
      n_years_for_training = 5,
      weeks_for_training = NULL
    ),

    # ... SRF ---------------------------------------------------------
    
    # Euromomo style Serfling
    # https://github.com/EuroMOMOnetwork/MOMO/blob/master/R/excess.R
    # AIC selection of seasonality
    'SRFcem', 'glm', list(
      models = list(
        formula(
          deaths_observed ~
            # log linear long term trend
            origin_weeks +
            # seasonality
            # full year period
            sin(2*pi/52*iso_week) +
            cos(2*pi/52*iso_week) +
            # half year period
            sin(2*pi/26*iso_week) +
            cos(2*pi/26*iso_week) +
            # adjustment for special weeks
            holiday3
        ),
        formula(
          deaths_observed ~
            # log linear long term trend
            origin_weeks +
            # seasonality
            # full year period
            sin(2*pi/52*iso_week) +
            cos(2*pi/52*iso_week) +
            # adjustment for special weeks
            holiday3
        ),
        formula(
          deaths_observed ~
            # log linear long term trend
            origin_weeks +
            # adjustment for special weeks
            holiday3
        )
      ),
      family = quasipoisson(link = 'log'),
      weeks_for_training = c(15:26, 36:45),
      n_years_for_training = NULL
    ),
    # Forecasting Serfling without exposures
    # AIC selection of seasonality
    'SRFc', 'glm', list(
      models = list(
        formula(
          deaths_observed ~
            # log linear long term trend
            origin_weeks +
            # seasonality
            # full year period
            sin(2*pi/52*iso_week) +
            cos(2*pi/52*iso_week) +
            # half year period
            sin(2*pi/26*iso_week) +
            cos(2*pi/26*iso_week) +
            # adjustment for special weeks
            holiday3
        ),
        formula(
          deaths_observed ~
            # log linear long term trend
            origin_weeks +
            # seasonality
            # full year period
            sin(2*pi/52*iso_week) +
            cos(2*pi/52*iso_week) +
            # adjustment for special weeks
            holiday3
        ),
        formula(
          deaths_observed ~
            # log linear long term trend
            origin_weeks +
            # adjustment for special weeks
            holiday3
        )
      ),
      family = quasipoisson(link = 'log'),
      weeks_for_training = NULL,
      n_years_for_training = NULL
    ),
    # Forecasting Serfling with exposures
    # AIC selection of seasonality
    'SRFr', 'glm', list(
      models = list(
        formula(
          deaths_observed ~
            # log linear long term trend
            origin_weeks +
            # seasonality
            # full year period
            sin(2*pi/52*iso_week) +
            cos(2*pi/52*iso_week) +
            # half year period
            sin(2*pi/26*iso_week) +
            cos(2*pi/26*iso_week) +
            # adjustment for special weeks
            holiday3 +
            # exposures
            offset(log(personweeks))
        ),
        formula(
          deaths_observed ~
            # log linear long term trend
            origin_weeks +
            # seasonality
            # full year period
            sin(2*pi/52*iso_week) +
            cos(2*pi/52*iso_week) +
            # adjustment for special weeks
            holiday3 +
            # exposures
            offset(log(personweeks))
        ),
        formula(
          deaths_observed ~
            # log linear long term trend
            origin_weeks +
            # adjustment for special weeks
            holiday3 +
            # exposures
            offset(log(personweeks))
        )
      ),
      family = quasipoisson(link = 'log'),
      weeks_for_training = NULL,
      n_years_for_training = NULL
    ),
    
    # ... GAM ---------------------------------------------------------
    
    # Gam without temperature
    'GAMr', 'gam', list(
      formula = formula(
        deaths_observed ~
          # log linear long term trend
          origin_weeks*stratum_id +
          # penalized cyclic spline for seasonality
          s(iso_week, bs = 'cp', k = 12, by = stratum_id) +
          # adjustment for special weeks
          holiday3*stratum_id +
          # exposures
          offset(log(personweeks))
      ),
      family = quasipoisson(link = 'log')
    ),
    # Gam with temperature
    'GAMrt', 'gam', list(
      formula = formula(
        deaths_observed ~
          # log linear long term trend
          origin_weeks*stratum_id +
          # penalized cyclic spline for seasonality
          s(iso_week, bs = 'cp', k = 12, by = stratum_id) +
          # temperature effect
          s(iso_week, bs = 'cp', k = 12, by = temperature_anomaly) +
          # adjustment for special weeks
          holiday3*stratum_id +
          # exposures
          offset(log(personweeks))
      ),
      family = quasipoisson(link = 'log')
    ),
    
    # ... LGM ---------------------------------------------------------
    
    # Kontis style LGM without temperature effect
    'LGMr', 'lgm', list(
      formula = formula(
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
          f(resid_iid,
            model = 'iid',
            hyper = list(prec = list(prior = 'loggamma', param = c(0.001, 0.001)))
          )
      )
    ),
    # Kontis style LGM with temperature effect
    'LGMrt', 'lgm', list(
      formula = formula(
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
          f(week_rw, tanomaly,
            model = 'rw1', cyclic = TRUE,
            hyper = list(prec = list(prior = 'loggamma', param = c(0.001, 0.001)))
          ) +
          f(resid_iid,
            model = 'iid',
            hyper = list(prec = list(prior = 'loggamma', param = c(0.001, 0.001)))
          )
      )
    ),
    # Kontis style LGM with temperature effect and different order residuals
    'LGMrt2', 'lgm', list(
      formula = formula(
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
    )
    
  )

# Export ----------------------------------------------------------

saveRDS(mod_para, file = cnst$path_mod_para)
