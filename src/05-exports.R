# Format fitted models for CPOP usage

# Init ------------------------------------------------------------

library(here); library(glue)
library(tidyverse)
library(openxlsx)

wd <- here()
setwd(wd)

# paths
paths <- list(
  dat = 'dat',
  out = 'out',
  init = 'src/00-init.R',
  fitted_models = 'out/fitted_models.rds'
)

dat <- list()

# Load fitted models ----------------------------------------------

dat$fitted_models <- readRDS(paths$fitted_models)

dat$expected_deaths <-
  dat$fitted_models %>%
  filter(
    model_id %in% c(1:10, 14:17), cv_id == 0
  ) %>%
  unnest(predictions) %>%
  mutate(deaths_predicted = round(deaths_predicted, 2)) %>%
  filter(iso_year == 2020, iso_week %in% 10:26) %>%
  select(
    model_name_short, obs_id,
    sex, age_start, age_width, iso_year, iso_week,
    region_name, region_iso,
    personweeks, deaths_observed, deaths_predicted
  ) %>%
  pivot_wider(names_from = model_name_short, values_from = deaths_predicted)

# Export ----------------------------------------------------------

saveRDS(dat$expected_deaths, file = glue('{paths$out}/expected_deaths_w10thru26_2020.rds'))

write_csv(dat$expected_deaths, file = glue('{paths$out}/expected_deaths_w10thru26_2020.csv'))

write.xlsx(dat$expected_deaths, glue('{paths$out}/expected_deaths_w10thru26_2020.xlsx'),
           keepNA = TRUE, na.string = '.',
           firstRow = TRUE, firstCol = TRUE)

