# Determine error and bias of various models of expected deaths

# Init ------------------------------------------------------------

library(here); library(glue)
library(tidyverse)
library(patchwork)

wd <- here(); setwd(wd)

# paths
path <- list(
  dat = 'dat',
  out = 'out',
  glob = 'src/00-glob.R',
  fitted_models = 'tmp/fitted_models.rds'
)

cnst <- list()
cnst <- within(cnst,{
  model_metadata = read_csv('src/model_metadata.csv')
  models_to_include = model_metadata[model_metadata$include==1,][['code']]
})

# global functions and constants
source(path$glob)

dat <- list()
fig <- list()

# Load data -------------------------------------------------------

# load data for cross validation
dat$fitted_models <- readRDS(path$fitted_models)
# fitted cross-validation series, only test data
dat$cv_test <-
  dat$fitted_models %>%
  filter(
    cv_id != 0,
    model_id %in% cnst$models_to_include,
    error_while_fit == FALSE
  ) %>%
  unnest(predictions) %>%
  filter(cv_sample == 'test') %>%
  select(-contains('deaths_sim'))

# Calculate errors on multiple aggregation levels -----------------

# residuals by week, age and sex
dat$residual_death_week_age_sex <-
  DeathsResiduals(
    dat$cv_test,
    weeks_since_test_start, deaths_observed, deaths_predicted, cv_id,
    model_id,
    region_iso, sex, age_group
  )

# residuals by week
dat$residual_death_week <-
  DeathsResiduals(
    dat$cv_test,
    weeks_since_test_start, deaths_observed, deaths_predicted, cv_id,
    model_id,
    region_iso
  )

# Figure errorbias ------------------------------------------------

PlotErrorsAndBias <- function (df_errors, error_measure, bias_measure, xlab) {
  
  require(tidyverse)
  
  Format <- function (x) { formatC(x, format = 'f', digits = 1) }
  
  observed_errors <-
    df_errors %>%
    select(
      model = model_id,
      error_measure = {{error_measure}},
      bias_measure = {{bias_measure}}
    )

  # summarise errors over countries and possibly strata and weeks
  summarised_errors <-
    observed_errors %>%
    group_by(model) %>%
    summarise(
      error_qlo = quantile(error_measure, p = 0.25, na.rm = TRUE),
      error_qmd = quantile(error_measure, p = 0.5, na.rm = TRUE),
      error_qhi = quantile(error_measure, p = 0.75, na.rm = TRUE),
      bias_qlo = quantile(bias_measure, p = 0.25, na.rm = TRUE),
      bias_qmd = quantile(bias_measure, p = 0.5, na.rm = TRUE),
      bias_qhi = quantile(bias_measure, p = 0.75, na.rm = TRUE)
    ) %>%
    left_join(cnst$model_metadata, c('model' = 'code')) %>%
    mutate(
      model =
        fct_reorder(model, order_1)
    )
  
  ynudge <- 0.2
  sizelarge <- 1
  sizesmall <- 0.6
  textsize <- 2
  sizeribbon <- 7
  
  fig <-
    summarised_errors %>%
    ggplot(aes(y = model, yend = model)) +
    # indicate rows
    geom_segment(
      aes(color = highlight),
      x = -Inf, xend = Inf, size = sizeribbon
    ) +
    geom_vline(xintercept = 0, color = 'grey50', size = 1.5) +
    # plot errors
    geom_segment(
      aes(x = error_qlo, xend = error_qhi),
      position = position_nudge(y = ynudge),
      size = sizesmall
    ) +
    geom_label(
      aes(x = error_qmd, label = Format(error_qmd)),
      position = position_nudge(y = ynudge),
      label.r = unit(0, 'pt'), size = textsize, fontface = 'italic',
      label.padding = unit(1, 'pt')
    ) +
    geom_text(
      aes(x = error_qlo-0.3, label = Format(error_qlo)),
      position = position_nudge(y = ynudge),
      hjust = 'right', size = textsize, fontface = 'italic'
    ) +
    geom_text(
      aes(x = error_qhi+0.3, label = Format(error_qhi)),
      position = position_nudge(y = ynudge),
      hjust = 'left', size = textsize, fontface = 'italic'
    ) +
    # plot bias
    geom_segment(
      aes(x = bias_qlo, xend = bias_qhi),
      position = position_nudge(y = -ynudge),
      size = sizelarge
    ) +
    geom_label(
      aes(x = bias_qmd, label = Format(bias_qmd)),
      position = position_nudge(y = -ynudge),
      label.r = unit(0, 'pt'), size = textsize, fontface = 'bold',
      label.padding = unit(1, 'pt')
    ) +
    geom_text(
      aes(x = bias_qlo-0.3, label = Format(bias_qlo)),
      position = position_nudge(y = -ynudge),
      hjust = 'right', size = textsize
    ) +
    geom_text(
      aes(x = bias_qhi+0.3, label = Format(bias_qhi)),
      position = position_nudge(y = -ynudge),
      hjust = 'left', size = textsize
    ) +
    # misc
    figspec$MyGGplotTheme(show_legend = FALSE) +
    scale_x_continuous(breaks = seq(-10, 10, 5)) +
    scale_color_brewer(type = 'qual', palette = 5) +
    labs(y = NULL, x = xlab) +
    coord_cartesian(clip = 'off', xlim = c(-10,10))
  
  list(errors = summarised_errors, fig = fig)
  
}

# error & bias, total, country level
fig$errorbias_a <-
  dat$residual_death_week$residual_summary %>%
  filter(weeks_since_test_start == 45) %>%
  PlotErrorsAndBias(
    error_measure = mape_cumdeath,
    bias_measure = mpe_cumdeath,
    xlab = NULL
  )

# error & bias, total, stratum level
fig$errorbias_b <-
  dat$residual_death_week_age_sex$residual_summary %>%
  filter(weeks_since_test_start == 45) %>%
  PlotErrorsAndBias(
    error_measure = mape_cumdeath,
    bias_measure = mpe_cumdeath,
    xlab = NULL
  )

# error & bias, weekly, country level
fig$errorbias_c <-
  dat$residual_death_week$residual_summary %>%
  PlotErrorsAndBias(
    error_measure = mape_cumdeath,
    bias_measure = mpe_cumdeath,
    xlab = NULL
  )

# error & bias, weekly, stratum level
fig$errorbias_d <-
  dat$residual_death_week_age_sex$residual_summary %>%
  PlotErrorsAndBias(
    error_measure = mape_cumdeath,
    bias_measure = mpe_cumdeath,
    xlab = 'MPE/MAPE'
  )

# assemble multi-panel figure 3
fig$errorbias <-
  fig$errorbias_a$fig +
  labs(subtitle = 'a. total deaths by country', y = 'Model') +
  fig$errorbias_b$fig +
  labs(subtitle = 'b. total deaths by country and stratum') +
  fig$errorbias_c$fig +
  labs(subtitle = 'c. weekly deaths by country') +
  fig$errorbias_d$fig +
  labs(subtitle = 'd. weekly deaths by country and stratum') +
  plot_layout(ncol = 2, byrow = TRUE)
fig$errorbias

ExportFigure(
  fig$errorbias, path = path$out, filename = 'errorbias',
  add_date = FALSE,
  device = 'pdf',
  width = figspec$fig_dims$width,
  height = figspec$fig_dims$width
)

# Per-timestep prediction error -----------------------------------

dat$weeklyerrors <-
  dat$residual_death_week$residual_summary %>%
  group_by(model_id, weeks_since_test_start) %>%
  summarise(mape = quantile(mape_death, na.rm = TRUE, p = 0.5)) %>%
  ungroup() %>%
  left_join(cnst$model_metadata, c('model_id' = 'code')) %>%
  mutate(model_id = fct_reorder(model_id, order_1))

dat$weeklyerrors_background <-
  dat$weeklyerrors %>%
  select(weeks_since_test_start, mape) %>%
  expand_grid(model_id = unique(dat$weeklyerrors$model_id))

fig$weeklyerrors <-
  dat$weeklyerrors %>%
  ggplot(aes(x = weeks_since_test_start + glob$start_of_test_iso_week, y = mape)) +
  geom_point(
    color = 'grey', size = 0.5,
    data = dat$weeklyerrors_background
  ) +
  geom_line(size = 1) +
  facet_wrap(~model_id) +
  figspec$MyGGplotTheme() +
  labs(x = 'Week of year', y = 'MAPE')

ExportFigure(
  fig$weeklyerrors, path = path$out, filename = 'weeklyerrors',
  add_date = FALSE,
  device = 'pdf',
  width = figspec$fig_dims$width,
  height = figspec$fig_dims$width*0.6
)
