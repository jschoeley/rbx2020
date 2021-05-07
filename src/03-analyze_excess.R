# Calculate excess deaths under various models

# Init ------------------------------------------------------------

library(here); library(glue)
library(tidyverse)
library(patchwork)

wd <- here(); setwd(wd)

# paths
path <- list(
  dat = 'dat',
  out = 'out',
  tmp = 'tmp',
  glob = 'src/00-glob.R',
  fitted_models = 'tmp/fitted_models.rds'
)

# global functions and constants
source(path$glob)

# constants
cnst <- list()
cnst <- within(cnst, {
  # analyze this many weeks into test-set
  n_weeks_into_test = glob$forecast_n_weeks
  # number of Poisson replications of weekly predicted death counts
  sim_rep = 100
  model_metadata = read_csv('src/model_metadata.csv')
  region_metadata = read_csv('src/region_metadata.csv')
  models_to_include = model_metadata[model_metadata$include==1,][['code']]
  # model to compare other model predictions to
  reference_model_name = quo(`AVRc5`)
})

dat <- list()
fig <- list()
tab <- list()

# Load data -------------------------------------------------------

# load data on predicted death counts
dat$fitted_models <- readRDS(path$fitted_models)

# Subset to predictions over COVID period -------------------------

# subset to total data series with test set starting early 2020
# only include models of interest
dat$predictions_covid <-
  dat$fitted_models %>%
  select(-model_para) %>%
  unnest(predictions) %>%
  filter(
    model_id %in% cnst$models_to_include,
    cv_id == 0,
    weeks_since_test_start >= 0,
    weeks_since_test_start < cnst$n_weeks_into_test
  )

# save memory
dat$fitted_models <- NULL

# Calculate excess deaths -----------------------------------------

# excess deaths by model and region
dat$excess_deaths <-
  dat$predictions_covid_replicate %>%
  ExcessDeaths(
    date = date,
    observed_deaths = deaths_observed,
    expected_deaths = deaths_predicted,
    simulated_deaths = deaths_simulated,
    sim_id = sim_id,
    model_id,
    region_iso
  )

# cumulative excess deaths by model and region, age and sex
dat$excess_deaths_age_sex <-
  dat$predictions_covid_replicate %>%
  ExcessDeaths(
    date = date,
    observed_deaths = deaths_observed,
    expected_deaths = deaths_predicted,
    simulated_deaths = deaths_simulated,
    sim_id = sim_id,
    model_id,
    region_iso, sex, age_group
  )

# save calculated excess
saveRDS(dat$excess_deaths, file = glue('{path$out}/excess_deaths.rds'))
saveRDS(dat$excess_deaths_age_sex, file = glue('{path$out}/excess_deaths_age_sex.rds'))

# save memory
dat$predictions_covid <- NULL
dat$predictions_covid_replicate <- NULL

#dat$excess_deaths <- readRDS(glue('{path$out}/excess_deaths.rds'))
#dat$excess_deaths_age_sex <- readRDS(glue('{path$out}/excess_deaths_age_sex.rds'))

# Countries where models disagree ---------------------------------

dat$excess_deaths %>%
  filter(date == max(date)) %>%
  group_by(region_iso) %>%
  summarise(
    # does any model give significant excess deaths?
    any_excess = any(cumpscore_q05 > 0),
    # do all models give significant excess deaths?
    all_excess = all(cumpscore_q05 > 0),
    # which models give significant excess deaths?
    excess_which = paste0(model_id[cumpscore_q05 > 0], collapse = ' '),
    # which models do not give significant excess deaths?
    no_excess_which = paste0(model_id[cumpscore_q05 <= 0], collapse = ' ')
  ) %>%
  filter(all_excess == FALSE)


# Examples --------------------------------------------------------

# Belgium end of year excess deaths
dat$excess_deaths %>%
  filter(region_iso == 'BE', date == max(date)) %>%
  select(model_id, cumexcess_expected, cumexcess_q05, cumexcess_q95)

# Figure excess ---------------------------------------------------

dat$fig_excess$excess_deaths_all <-
  dat$excess_deaths %>%
  left_join(cnst$region_metadata, c('region_iso' = 'region_code')) %>%
  select(region_iso, region_name, model_id, date,
         pscore_expected, pscore_q05, pscore_q95,
         cumpscore_expected, cumpscore_q05, cumpscore_q95)
dat$fig_excess$excess_deaths_last <-
  dat$fig_excess$excess_deaths_all %>% 
  filter(date == max(date)) %>%
  group_by(region_name, date) %>%
  summarise(
    min = min(cumpscore_expected),
    avg = mean(cumpscore_expected),
    max = max(cumpscore_expected),
    which_min = model_id[which.min(cumpscore_expected)],
    which_max = model_id[which.max(cumpscore_expected)]
  )
dat$fig_excess$label_colors <- c(positive = '#471227', negative = '#002B36')
dat$fig_excess$background_colors  <- c(positive = '#fd9476', negative = '#4bc4fc')

fig$excess <-
  dat$fig_excess$excess_deaths_all %>%
  ggplot(aes(x = date, group = model_id)) +
  # background
  geom_rect(
    aes(xmin = min(date), xmax = max(date), ymin = 0, ymax = Inf),
    fill = dat$fig_excess$background_colors['positive']
  ) +
  geom_rect(
    aes(xmin = min(date), xmax = max(date), ymin = 0, ymax = -Inf),
    fill = dat$fig_excess$background_colors['negative']
  ) +
  geom_hline(yintercept = 0, color = 'black') +
  # weekly excess
  geom_point(
    aes(y = pscore_expected*1e2),
    size = 0.1, alpha = 0.3
  ) +
  # cumulative excess
  geom_line(aes(y = cumpscore_expected*1e2), alpha = 0.5) +
  # end of year excess labels
  geom_label(
    aes(
      x = date, y = min*1e2-2,
      label =
        paste(
          formatC(min*1e2, format = 'd', flag = '+'),
          which_min
        ),
      fill = ifelse(sign(min)>0,
                    dat$fig_excess$label_colors['positive'],
                    dat$fig_excess$label_colors['negative'])
    ),
    family = 'robotocondensed',
    hjust = 1, size = 3, vjust = 1,
    label.padding = unit(1, 'pt'), label.size = 0, label.r = unit(0, 'pt'),
    data = dat$fig_excess$excess_deaths_last,
    inherit.aes = FALSE,
    color = 'white', alpha = 0.7
  ) +
  geom_label(
    aes(
      x = date, y = max*1e2+2,
      label =
        paste(
          formatC(max*1e2, format = 'd', flag = '+'),
          which_max
        ),
      fill = ifelse(sign(max)>0,
                     dat$fig_excess$label_colors['positive'],
                     dat$fig_excess$label_colors['negative'])
    ),
    family = 'robotocondensed',
    hjust = 1, size = 3, vjust = 0,
    label.padding = unit(1, 'pt'), label.size = 0, label.r = unit(0, 'pt'),
    data = dat$fig_excess$excess_deaths_last,
    inherit.aes = FALSE,
    color = 'white', alpha = 0.7
  ) +
  # across countries
  facet_wrap(~region_name, ncol = 4) +
  # misc
  scale_x_date(date_breaks = '8 weeks', date_labels = '%b') +
  scale_y_continuous(limits = c(-20, 60)) +
  scale_color_identity() +
  scale_fill_identity() +
  figspec$MyGGplotTheme(grid = '', panel_border = TRUE, axis = '', family = 'roboto') +
  theme(panel.ontop = TRUE) +
  labs(
    x = NULL, y = '(Cumulative) weekly excess death percentage'
  ) +
  coord_cartesian(expand = FALSE)
fig$excess

ExportFigure(
  fig$excess, path = path$out, filename = 'excess',
  device = 'pdf',
  add_date = FALSE,
  width = figspec$fig_dims$width,
  height = figspec$fig_dims$width*1.3
)

# Figure rank -----------------------------------------------------

dat$fig_rank$main <-
  dat$excess_deaths %>%
  left_join(cnst$model_metadata, by = c('model_id' = 'code')) %>%
  mutate(model_id = fct_reorder(model_id, order_1)) %>%
  filter(date == max(date)) %>%
  select(region_iso, model_id, cumpscore_expected) %>%
  arrange(model_id) %>%
  group_by(model_id) %>%
  mutate(rank = rank(cumpscore_expected)) %>%
  ungroup()
dat$fig_rank$labels <-
  dat$fig_rank$main %>%
  filter(model_id == levels(dat$fig_rank$main$model_id)[1])
dat$fig_rank$colors <- rep(c('#131917', '#2D424A', '#5F6468', '#B4786B', '#612A11'), 4)
names(dat$fig_rank$colors) <-
  dat$fig_rank$labels %>% arrange(rank) %>% pull(region_iso)

fig$rank <-
  dat$fig_rank$main %>%
  ggplot(aes(
    x = model_id, y = rank,
    group = region_iso, color = region_iso
  )) +
  geom_ribbon(
    aes(xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 5.5),
    fill = 'grey90', color = NA
  ) +
  geom_ribbon(
    aes(xmin = -Inf, xmax = Inf, ymin = 10.5, ymax = 15.5),
    fill = 'grey90', color = NA
  ) +
  geom_line(size = 2) +
  geom_point(
    aes(color = region_iso),
    fill = 'white', size = 3, shape = 21
  ) +
  geom_text(
    aes(label = region_iso, x = 0.9),
    size = 2.3,
    hjust = 1,
    family = 'robotocondensed', fontface = 'bold',
    data = dat$fig_rank$labels
  ) +
  scale_x_discrete(expand = c(0.07, 0.01)) +
  scale_y_continuous(breaks = 1:100, expand = c(0.01, 0)) +
  scale_color_manual(
    values = dat$fig_rank$colors
  ) +
  labs(y = 'Country rank of P-score', x = 'Model') +
  figspec$MyGGplotTheme(show_legend = FALSE, grid = '', family = 'roboto') +
  coord_cartesian(clip = 'off')
fig$rank

ExportFigure(
  fig$rank, path = path$out, filename = 'rank',
  device = 'pdf',
  width = figspec$fig_dims$width,
  height = figspec$fig_dims$width*0.8
)

# Figure baselinediff ---------------------------------------------

PlotLogratios <- function (
  logratios,
  xlab = 'Percent difference'
) {
  
  require(tidyverse)
    
  Format <- function (x) { formatC(x, format = 'f', digits = 1) }

  sizelarge <- 1
  textsize <- 2
  sizeribbon <- 8
  
  fig <-
    logratios %>%
    ggplot(aes(y = model_id, yend = model_id)) +
    # indicate rows
    geom_segment(
      aes(color = highlight),
      x = -Inf, xend = Inf, size = sizeribbon
    ) +
    geom_vline(xintercept = 0, color = 'grey50', size = 1.5) +
    # plot errors
    geom_segment(
      aes(x = qlo, xend = qhi),
      size = sizelarge
    ) +
    geom_label(
      aes(x = qmd, label = Format(qmd)),
      label.r = unit(0, 'pt'), size = textsize,
      label.padding = unit(1, 'pt')
    ) +
    geom_text(
      aes(x = qlo-0.05, label = Format(qlo)),
      hjust = 'right', size = textsize
    ) +
    geom_text(
      aes(x = qhi+0.05, label = Format(qhi)),
      hjust = 'left', size = textsize
    ) +
    # misc
    scale_x_continuous(breaks = seq(-10, 10, 5)) +
    coord_cartesian(clip = 'off', xlim = c(-10,10)) +
    figspec$MyGGplotTheme(show_legend = FALSE) +
    scale_color_brewer(type = 'qual', palette = 5)
  
  list(dat = logratios, fig = fig)
  
}

PartialSpread <- function (df, name, value, strata, spread) {
  df %>%
    select(name, value, strata) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    pivot_longer(c(-all_of(strata), -all_of(spread)), names_to = name)
}

dat$fig_baselinediff_a <-
  dat$excess_deaths %>%
  filter(date == max(date)) %>%
  PartialSpread('model_id','cumbaseline_expected',c('region_iso','date'),
                rlang::as_name(cnst$reference_model_name)) %>%
  mutate(r = (value-!!cnst$reference_model_name)/!!cnst$reference_model_name*100) %>%
  group_by(model_id) %>%
  summarise(
    qlo = quantile(r, p = 0.25, na.rm = TRUE),
    qmd = quantile(r, p = 0.5, na.rm = TRUE),
    qhi = quantile(r, p = 0.75, na.rm = TRUE)
  ) %>%
  left_join(cnst$model_metadata, by = c('model_id' = 'code')) %>%
  mutate(model_id = fct_reorder(model_id, order_1))

fig$baselinediff_a <-
  PlotLogratios(dat$fig_baselinediff_a)

dat$fig_baselinediff_b <-
  dat$excess_deaths %>%
  PartialSpread('model_id','baseline_expected',c('region_iso','date'),
                rlang::as_name(cnst$reference_model_name)) %>%
  mutate(r = (value-!!cnst$reference_model_name)/!!cnst$reference_model_name*100) %>%
  group_by(model_id) %>%
  summarise(
    qlo = quantile(r, p = 0.25, na.rm = TRUE),
    qmd = quantile(r, p = 0.5, na.rm = TRUE),
    qhi = quantile(r, p = 0.75, na.rm = TRUE)
  ) %>%
  left_join(cnst$model_metadata, by = c('model_id' = 'code')) %>%
  mutate(model_id = fct_reorder(model_id, order_1))

fig$baselinediff_b <- PlotLogratios(dat$fig_baselinediff_b)

dat$fig_baselinediff_c <-
  dat$excess_deaths_age_sex %>%
  filter(date == max(date)) %>%
  PartialSpread('model_id','baseline_expected',c('region_iso','sex','age_group','date'),
                rlang::as_name(cnst$reference_model_name)) %>%
  mutate(r = (value-!!cnst$reference_model_name)/!!cnst$reference_model_name*100) %>%
  group_by(model_id) %>%
  summarise(
    qlo = quantile(r, p = 0.25, na.rm = TRUE),
    qmd = quantile(r, p = 0.5, na.rm = TRUE),
    qhi = quantile(r, p = 0.75, na.rm = TRUE)
  ) %>%
  left_join(cnst$model_metadata, by = c('model_id' = 'code')) %>%
  mutate(model_id = fct_reorder(model_id, order_1))

fig$baselinediff_c <- PlotLogratios(dat$fig_baselinediff_c)

dat$fig_baselinediff_d <-
  dat$excess_deaths_age_sex %>%
  PartialSpread('model_id','baseline_expected',c('region_iso','sex','age_group','date'),
                rlang::as_name(cnst$reference_model_name)) %>%
  mutate(r = (value-!!cnst$reference_model_name)/!!cnst$reference_model_name*100) %>%
  group_by(model_id) %>%
  summarise(
    qlo = quantile(r, p = 0.25, na.rm = TRUE),
    qmd = quantile(r, p = 0.5, na.rm = TRUE),
    qhi = quantile(r, p = 0.75, na.rm = TRUE)
  ) %>%
  left_join(cnst$model_metadata, by = c('model_id' = 'code')) %>%
  mutate(model_id = fct_reorder(model_id, order_1))

fig$baselinediff_d <- PlotLogratios(dat$fig_baselinediff_d)

fig$baselinediff <-
  fig$baselinediff_a$fig +
  labs(subtitle = 'a. total deaths by country', y = 'Model', x = NULL) +
  fig$baselinediff_b$fig +
  labs(subtitle = 'b. total deaths by country and stratum', y = NULL, x = NULL) +
  fig$baselinediff_c$fig +
  labs(subtitle = 'c. weekly deaths by country', y = NULL, x = NULL) +
  fig$baselinediff_d$fig +
  labs(subtitle = 'd. weekly deaths by country and stratum', y = NULL,
       x = 'Percent difference of predicted deaths against AVRc5') +
  plot_layout(ncol = 2, byrow = TRUE)
fig$baselinediff

ExportFigure(
  fig$baselinediff, path = path$out, filename = 'baselinediff',
  device = 'pdf',
  add_date = FALSE,
  width = figspec$fig_dims$width,
  height = figspec$fig_dims$width
)

# Figure rank over strata -----------------------------------------

dat$fig_rank_strata$main <-
  dat$excess_deaths_age_sex %>%
  left_join(cnst$model_metadata, by = c('model_id' = 'code')) %>%
  mutate(model_id = fct_reorder(model_id, order_1)) %>%
  filter(date == max(date)) %>%
  select(region_iso, sex, age_group, model_id, cumpscore_expected) %>%
  arrange(model_id) %>%
  group_by(model_id, sex, age_group) %>%
  mutate(rank = rank(cumpscore_expected)) %>%
  ungroup()
dat$fig_rank_strata$labels <-
  dat$fig_rank_strata$main %>%
  filter(model_id == levels(dat$fig_rank_strata$main$model_id)[1])
dat$fig_rank_strata$main <-
  left_join(
    dat$fig_rank_strata$main,
    dat$fig_rank_strata$labels %>%
      select(-cumpscore_expected, -model_id, initial_rank = rank)
  )
dat$fig_rank_strata$colors <- rep(c('#131917', '#2D424A', '#5F6468', '#B4786B', '#612A11'), 4)

fig$rank_strata <-
  dat$fig_rank_strata$main %>%
  ggplot(aes(
    x = model_id, y = rank,
    group = region_iso, color = as.factor(initial_rank)
  )) +
  geom_ribbon(
    aes(xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 5.5),
    fill = 'grey90', color = NA
  ) +
  geom_ribbon(
    aes(xmin = -Inf, xmax = Inf, ymin = 10.5, ymax = 15.5),
    fill = 'grey90', color = NA
  ) +
  geom_line(size = 1.5) +
  geom_point(
    fill = 'white', size = 2, shape = 21
  ) +
  geom_text(
    aes(x = 0.8, label = region_iso, color = as.factor(rank)),
    size = 2.3,
    hjust = 1,
    family = 'robotocondensed', fontface = 'bold',
    data = dat$fig_rank_strata$labels
  ) +
  facet_grid(age_group~sex) +
  scale_y_continuous(breaks = 1:100) +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_color_manual(
    values = dat$fig_rank_strata$colors
  ) +
  labs(y = 'Country rank of P-score', x = 'Model') +
  figspec$MyGGplotTheme(show_legend = FALSE, grid = '') +
  coord_cartesian(clip = 'off')
fig$rank_strata

ExportFigure(
  fig$rank_strata, path = path$out, filename = 'rank_strata',
  device = 'pdf',
  width = figspec$fig_dims$width,
  height = figspec$fig_dims$width*1.4, scale = 1.2
)

# Figure excess over strata ---------------------------------------

fig$excess_strata <-
  dat$excess_deaths_age_sex %>%
  group_by(age_group, sex) %>%
  group_map(~{

    excess_deaths_all <-
      .x %>%
      left_join(cnst$region_metadata, c('region_iso' = 'region_code')) %>%
      select(region_iso, region_name, model_id, date,
             pscore_expected, pscore_q05, pscore_q95,
             cumpscore_expected, cumpscore_q05, cumpscore_q95)
    excess_deaths_last <-
      excess_deaths_all %>% 
      filter(date == max(date)) %>%
      group_by(region_name, date) %>%
      summarise(
        min = min(cumpscore_expected, na.rm = TRUE),
        avg = mean(cumpscore_expected, na.rm = TRUE),
        max = max(cumpscore_expected, na.rm = TRUE),
        which_min = model_id[which.min(cumpscore_expected)],
        which_max = model_id[which.max(cumpscore_expected)]
      )
    label_colors <- c(positive = '#471227', negative = '#002B36')
    background_colors  <- c(positive = '#fd9476', negative = '#4bc4fc')
    
    fig <-
      excess_deaths_all %>%
      ggplot(aes(x = date, group = model_id)) +
      geom_rect(
        aes(xmin = min(date), xmax = max(date), ymin = 0, ymax = Inf),
        fill = background_colors['positive']
      ) +
      geom_rect(
        aes(xmin = min(date), xmax = max(date), ymin = 0, ymax = -Inf),
        fill = background_colors['negative']
      ) +
      geom_hline(yintercept = 0, color = 'black') +
      geom_point(
        aes(y = pscore_expected*1e2),
        size = 0.1
      ) +
      geom_line(aes(y = cumpscore_expected*1e2), alpha = 0.5) +
      # end of year excess labels
      geom_label(
        aes(
          x = date, y = min*1e2-2,
          label =
            paste(
              formatC(min*1e2, format = 'd', flag = '+'),
              which_min
            ),
          fill = ifelse(sign(min)>0, label_colors['positive'], label_colors['negative'])
        ),
        family = 'robotocondensed',
        hjust = 1, size = 3, vjust = 1,
        label.padding = unit(1, 'pt'), label.size = 0, label.r = unit(0, 'pt'),
        data = excess_deaths_last,
        inherit.aes = FALSE,
        color = 'white', alpha = 0.7
      ) +
      geom_label(
        aes(
          x = date, y = max*1e2+2,
          label =
            paste(
              formatC(max*1e2, format = 'd', flag = '+'),
              which_max
            ),
          fill = ifelse(sign(max)>0, dat$fig_excess$label_colors['positive'], dat$fig_excess$label_colors['negative'])
        ),
        family = 'robotocondensed',
        hjust = 1, size = 3, vjust = 0,
        label.padding = unit(1, 'pt'), label.size = 0, label.r = unit(0, 'pt'),
        data = excess_deaths_last,
        inherit.aes = FALSE,
        color = 'white', alpha = 0.7
      ) +
      scale_x_date(date_breaks = '8 weeks', date_labels = '%b') +
      scale_color_identity() +
      scale_fill_identity() +
      figspec$MyGGplotTheme(grid = '', panel_border = TRUE, axis = '') +
      labs(
        x = NULL, y = '(Cumulative) weekly excess death percentage'
      ) +
      facet_wrap(~region_name, ncol = 4) +
      coord_cartesian(expand = FALSE) +
      theme(panel.ontop = TRUE)

    list(fig = fig, name = glue('excess_{.y$sex}_{.y$age_group}'))
    
  })

names(fig$excess_strata) <- unlist(map(fig$excess_strata, ~.x[['name']]))
fig$excess_strata <- map(fig$excess_strata, ~.x[['fig']])
ExportFiguresFromList(
  fig$excess_strata,
  path = path$out,
  device = 'pdf',
  add_date = FALSE,
  width = figspec$fig_dims$width,
  height = figspec$fig_dims$width*1.3
)
