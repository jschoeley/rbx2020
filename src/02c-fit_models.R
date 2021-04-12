# Fit expected deaths models

# Init ------------------------------------------------------------

set.seed(1987)

library(here); library(glue)
library(tidyverse)
library(foreach)
library(doParallel)

# Constants -------------------------------------------------------

wd <- here()
setwd(wd)

# paths
paths <- list(
  dat = 'dat',
  out = 'out',
  tmp = 'tmp',
  mod_spec = 'tmp/mod_spec.rds',
  mod_para = 'tmp/mod_para.rds',
  country_metadata = 'src/region_metadata.csv',
  glob = 'src/00-glob.R',
  log = 'tmp/fitted_models_log.txt',
  mocy_cv = 'tmp/mocy_cv.rds'
)

# global functions and constants
source(paths$glob)

# model specifications
ModSpec <- readRDS(paths$mod_spec)
# model parametrizations
mod_para <- readRDS(paths$mod_para)

# info on regions
region_metadata <- read_csv(paths$country_metadata)

# constants
cnst <- list()
cnst <- within(cnst,{
  # which countries to analyze?
  countries =
    region_metadata %>%
    filter(europe_analysis == 1) %>%
    pull(region_code)
  # how many threads used to fit models?
  cpu_nodes = 8
})

dat <- list()
fig <- list()

# setup parallel computation
cnst$cl <- makeCluster(cnst$cpu_nodes, outfile = paths$log)
registerDoParallel(cnst$cl)

# Load data -------------------------------------------------------

# load data for cross validation
mocy_cv <- readRDS(paths$mocy_cv)

# Prepare data for fit --------------------------------------------

dat$cv_data <-
  mocy_cv %>%
  # filter to countries of interest and add variables of interest
  filter(region_iso %in% cnst$countries) %>%
  ungroup()

# Fit models and predict ------------------------------------------

# merge CV data with model specs
dat$fit_data <-
  dat$cv_data %>%
  nest(
    # those are the dimensions to loop over, e.g.
    # separate fits by country and cv_id
    data = c(-region_iso, -cv_id)
  ) %>%
  expand_grid(
    nest(
      mod_para,
      model_para =
        c(-model_id, -model_class)
    )
  )

# fit models
dat$fitted_models <-
  foreach(
    x = iter(dat$fit_data, by = 'row'),
    .combine = bind_rows,
    .packages = c('dplyr', 'tidyr', 'mgcv', 'compositions', 'INLA')
  ) %dopar% {suppressPackageStartupMessages({
    
    cat(Sys.time(), ' Fit ', x$region_iso,
        ' on CV set ', x$cv_id,
        ' for ', x$model_id, '\n', sep = '')
    
    # prepare training and prediction data
    input_dat <- x[,'data'][[1]][[1]]
    dat_train <- filter(input_dat, cv_sample == 'training')
    dat_pred <- input_dat
    
    # fit models and capture errors
    result <- tryCatch({
      
      # ... Average mortality models ----------------------------------------
      
      if (x$model_class == 'avg') {
        
        # fit model and predict from it
        model_fit_and_predictions <- ModSpec$ModSimpleAverages(
          df_training = dat_train, df_prediction = dat_pred,
          week_name = iso_week, death_name = deaths_observed,
          exposure_name = personweeks, year_name = iso_year,
          n_years = x$model_para[[1]][[1]][[1]][['n_years']],
          stat = x$model_para[[1]][[1]][[1]][['stat']],
          sex, age_group
        )
        
      }

      # ... Time series models ----------------------------------------------

      if (x$model_class == 'ts') {
        
        # fit model and predict from it
        model_fit_and_predictions <- ModSpec$ModTS(
          df = dat_pred,
          sample_name = cv_sample, deaths_name = deaths_observed,
          date_name = date, stratum_id_name = stratum_id,
          model = x$model_para[[1]][[1]][[1]][['model']],
          lags = x$model_para[[1]][[1]][[1]][['lags']],
          persistence = x$model_para[[1]][[1]][[1]][['persistence']],
          distribution = x$model_para[[1]][[1]][[1]][['distribution']],
          initial = x$model_para[[1]][[1]][[1]][['initial']],
          orders = x$model_para[[1]][[1]][[1]][['orders']]
        )
        
      }
      
      # ... Serfling models -------------------------------------------------
      
      if (x$model_class == 'serfling') {
        
        model_fit_and_predictions <- ModSpec$ModSerfling(
          df_training = dat_train, df_prediction = dat_pred,
          formula = x$model_para[[1]][[1]][[1]][['formula']],
          family = x$model_para[[1]][[1]][[1]][['family']],
          name_week = iso_week,
          weeks_for_training = x$model_para[[1]][[1]][[1]][['weeks_for_training']],
          name_year = iso_year,
          n_years_for_training = x$model_para[[1]][[1]][[1]][['n_years_for_training']]
        )
        
      }
      
      # ... GAM models ------------------------------------------------------
      
      if (x$model_class == 'gam') {
        
        model_fit_and_predictions <- ModSpec$ModGAM(
          df_training = dat_train, df_prediction = dat_pred,
          formula = x$model_para[[1]][[1]][[1]][['formula']],
          family = x$model_para[[1]][[1]][[1]][['family']],
          method = x$model_para[[1]][[1]][[1]][['method']]
        )
        
      }

      # ... INLA models -------------------------------------------------

      if (x$model_class == 'inla') {
        
        model_fit_and_predictions <- ModSpec$ModINLA(
          df_traintest = dat_pred,
          inla_formula = x$model_para[[1]][[1]][[1]][['inla_formula']],
          stratum_name = stratum_id, sample_name = cv_sample,
          death_name = deaths_observed, exposure_name = personweeks,
          tanomaly_name = temperature_anomaly, week_name = epi_week,
          time_name = origin_weeks, holiday_name = holiday3,
          # because the iteration is already across multiple cores
          # make INLA only use a single core, this actually speeds-up
          # the whole fitting procedure
          threads = 1
        )
        
      }
      
      # ... Assemble results ------------------------------------------------
      
      result_if_no_error <- bind_cols(
        x,
        tibble(
          #fitted_model = list(model_fit_and_predictions$model),
          predictions = list(model_fit_and_predictions$predictions),
          error_while_fit = FALSE,
          error_message = NA
        )
      )
      return(result_if_no_error)
      
    },
    
    # if the fitting produces and error
    error = function(e) {
      cat(Sys.time(), ' Error', x$region_iso, ' on CV set ', x$cv_id,
          ' for ', x$model_id, ': ', geterrmessage(), '\n')
      result_if_error <- bind_cols(x, tibble(
        #fitted_model = list(NA),
        predictions = list(dat_pred %>% mutate(deaths_predicted = NA)),
        error_while_fit = TRUE,
        error_message = geterrmessage()
      ))
      return(result_if_error)
    })
    
    return(result)
    
  })}

stopCluster(cnst$cl)

# Plot observed vs. fitted ----------------------------------------

dat$fitted_models %>%
  #filter(model_id == '5-year avg. weekly mortality') %>%
  group_by(region_iso, model_id) %>%
  group_walk(~{
    
    fig_dat <-
      # predictions for single country and model
      .x %>%
      unnest(predictions) %>%
      group_by(cv_id, date, cv_sample) %>%
      summarise(
        deaths_observed = sum(deaths_observed),
        deaths_predicted = sum(deaths_predicted)
      )

    fig[[paste(.y[[1]], .y[[2]])]] <<-
      fig_dat %>%
      ggplot(aes(x = date)) +
      geom_point(aes(color = cv_sample, y = deaths_observed),
                 size = 0.3) +
      geom_line(aes(y = deaths_predicted, alpha = cv_sample),
                color = 'red') +
      scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
      scale_alpha_manual(values = c(training = 0.3, test = 1)) +
      scale_color_manual(values = figspec$colors$sample) +
      facet_grid(cv_id~'') +
      guides(color = 'none', alpha = 'none') +
      figspec$MyGGplotTheme(grid = 'xy') +
      labs(
        x = NULL, y = 'Weekly deaths',
        title = paste(.y[[1]], .y[[2]])
      )
  })

# Exports ---------------------------------------------------------

saveRDS(dat$fitted_models, file = 'tmp/fitted_models.rds', compress = 'xz')

ggsave(
  filename = 'fitted_vs_observed.pdf',
  path = paths$tmp,
  plot = gridExtra::marrangeGrob(fig, nrow=1, ncol=1), 
  width = 15, height = 9
)
