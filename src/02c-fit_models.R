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
  # number of draws from predicted distribution of deaths
  nsim = 500
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

# filter to countries of interest and add variables of interest
dat$cv_data <-
  mocy_cv %>%
  filter(region_iso %in% cnst$countries) %>%
  ungroup()

# merge data with model specs
dat$fit_data <-
  dat$cv_data %>%
  nest(data = c(-region_iso, -cv_id)) %>%
  expand_grid(mod_para)

# Fit models and predict ------------------------------------------

# fit models
dat$fitted_models <-
  foreach(
    x = iter(dat$fit_data, by = 'row'),
    .combine = bind_rows,
    .packages = c('dplyr', 'tidyr', 'mgcv', 'INLA')
  ) %dopar% {suppressPackageStartupMessages({
    
    cat(Sys.time(), ' Fit ', x$region_iso,
        ' on CV set ', x$cv_id,
        ' for ', x$model_id, '\n', sep = '')
    
    # extract input data
    input_dat <- x[,'data'][[1]][[1]]
    # model parameterization
    model_para <- x$model_para[[1]]
    
    # fit models and capture errors
    result <- tryCatch({
      
      # GLM
      if (x$model_class == 'glm') {
        
        predictions <- ModSpec$SerflingGLM(
          df = input_dat,
          models = model_para$models,
          family = model_para$family,
          col_sample = 'cv_sample',
          col_stratum = 'stratum_id',
          weeks_for_training = model_para$weeks_for_training,
          col_week = 'iso_week',
          n_years_for_training = model_para$n_years_for_training,
          col_year = 'iso_year',
          nsim = cnst$nsim, simulate_beta = TRUE, simulate_y = TRUE
        )
        
      }
      
      # GAM
      if (x$model_class == 'gam') {
        
        predictions <- ModSpec$CountGAM(
          df = input_dat, formula = model_para$formula,
          family = model_para$family,
          col_sample = 'cv_sample',
          nsim = cnst$nsim, simulate_beta = TRUE, simulate_y = TRUE
        )
        
      }

      # LGM
      if (x$model_class == 'lgm') {
        
        predictions <- ModSpec$KontisLGM(
          df = input_dat, formula1 = model_para$formula,
          formula2 = NULL, formula2_threshold = NULL,
          col_stratum = 'stratum_id', col_sample = 'cv_sample',
          col_death = 'deaths_observed', col_exposure = 'personweeks',
          col_tanomaly = 'temperature_anomaly', col_week = 'iso_week',
          col_time = 'origin_weeks', col_holiday = 'holiday3',
          nsim = cnst$nsim,
          weeks_for_training_within_year = NULL,
          weeks_for_training_pre_test = NULL,
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
          predictions = list(predictions),
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
      # return same object as fitted model, but with NA predictions
      input_dat[,c('deaths_predicted',paste0('deaths_sim', 1:100))] <- NA
      result_if_error <- bind_cols(x, tibble(
        predictions = list(input_dat),
        error_while_fit = TRUE,
        error_message = geterrmessage()
      ))
      return(result_if_error)

    }) # end of tryCatch()
    
    return(result)
    
  })}

stopCluster(cnst$cl)

# remove superfluous data column from fitted_models
# input data is part of predictions
dat$fitted_models <-
  dat$fitted_models %>% select(-data)

# Plot observed vs. fitted ----------------------------------------

dat$fitted_models %>%
  filter(!error_while_fit) %>%
  group_by(region_iso, model_id) %>%
  group_walk(~{
    
    predictions <- unnest(.x, predictions)
    expected_deaths <-
      predictions %>%
      group_by(cv_id, date, cv_sample) %>%
      summarise(
        deaths_observed = sum(deaths_observed),
        deaths_predicted = sum(deaths_predicted)
      ) %>%
      ungroup()
    simulated_deaths <-
      predictions %>%
      pivot_longer(cols = starts_with('deaths_sim'),
                   names_to = 'sim_id', values_to = 'deaths_sim') %>%
      group_by(cv_id, date, cv_sample, sim_id) %>%
      summarise(
        deaths_sim = sum(deaths_sim)
      ) %>%
      group_by(cv_id, date, cv_sample) %>%
      summarise(
        q05 = quantile(deaths_sim, 0.05, na.rm = TRUE),
        q95 = quantile(deaths_sim, 0.95, na.rm = TRUE)
      ) %>%
      ungroup()
    
    fig_dat <- left_join(expected_deaths, simulated_deaths)
    
    fig[[paste(.y[[1]], .y[[2]])]] <<-
      fig_dat %>%
      ggplot(aes(x = date)) +
      geom_ribbon(aes(ymin = q05, ymax = q95),
                  fill = 'grey70', color = NA) +
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
