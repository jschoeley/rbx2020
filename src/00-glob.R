# Global constants and functions

# Global constants ------------------------------------------------

glob <- list()
glob <- within(glob, {

  # iso-week which starts epi-year
  start_of_epi_year_iso_week = 27
  # iso-week which starts test period
  start_of_test_iso_week = 8
  # weeks to forecast
  forecast_n_weeks = 52-start_of_test_iso_week+1 # forecast until end of year
  
})

# Figure specification --------------------------------------------

# fonts
library(showtext)
font_add_google('Roboto', 'roboto')
font_add_google('Roboto Condensed', 'robotocondensed')
showtext_auto()

figspec <- list()
figspec <- within(figspec, {
  
  # color coding
  colors = list(
    sample =
      c(
        training = "grey30",
        test = "red"
      ),
    sex =
      c(
        `Male` = "#004B87",
        `Female` = "#c60c30"
      )
  )
  
  # figure dimensions in mm
  fig_dims = list(width = 180)
  
  # ggplot theme
  # ggplot theme by Jonas Schöley
  MyGGplotTheme <-
    function (
      size = 8,
      family = 'roboto',
      scaler = 1,
      axis = 'x',
      panel_border = FALSE,
      grid = 'y',
      minor_grid = '',
      show_legend = TRUE,
      ar = NA,
      axis_title_just = 'rt',
      axis_ticks = TRUE
    ) {
      
      size_med = size*scaler
      size_sml = round(size*0.7)*scaler
      base_linesize = 0.3*scaler
      
      # justification of axis titles
      xj <- switch(tolower(substr(axis_title_just, 1, 1)), b = 0, 
                   l = 0, m = 0.5, c = 0.5, r = 1, t = 1)
      yj <- switch(tolower(substr(axis_title_just, 2, 2)), b = 0, 
                   l = 0, m = 0.5, c = 0.5, r = 1, t = 1)
      
      list(
        theme_minimal(base_size = size_med, base_family = family),
        theme(
          # basic
          text = element_text(color = 'black'),
          line = element_line(size = base_linesize, lineend = 'square'),
          # axis
          axis.title = element_text(size = size_med, face = 'bold'),
          axis.title.x = element_text(hjust = xj),
          axis.title.y = element_text(hjust = yj),
          axis.title.y.right = element_text(hjust = yj, angle = 90),
          axis.text = element_text(size = size_med, color = 'black'),
          # strips
          strip.text = element_text(color = 'black', size = size_med),
          strip.background = element_blank(),
          # plot
          title = element_text(face = 'bold'),
          plot.subtitle = element_text(color = 'black', size = size_med, face = 'bold'),
          plot.caption = element_text(color = 'black', size = size_sml, face = 'plain'),
          plot.background = element_blank(),
          panel.background = element_blank(),
          #plot.margin = unit(c(1, 0.1, 0.5, 0.5), units = 'mm'),
          # grid
          panel.grid = element_blank()
        ),
        if (isTRUE(axis_ticks)) {
          theme(axis.ticks = element_line(size = rel(0.5), color = 'black'))
        },
        if (identical(grid, 'y')) {
          theme(panel.grid.major.y =
                  element_line(size = base_linesize, linetype = 3, color = 'grey80'))
        },
        if (identical(grid, 'x')) {
          theme(panel.grid.major.x =
                  element_line(size = base_linesize, linetype = 3, color = 'grey80'))
        },
        if (identical(grid, 'xy') | identical(grid, 'yx')) {
          theme(panel.grid.major.y =
                  element_line(size = base_linesize, linetype = 3, color = 'grey80'),
                panel.grid.major.x =
                  element_line(size = base_linesize, linetype = 3, color = 'grey80'))
        },
        if (identical(minor_grid, 'y')) {
          theme(panel.grid.minor.y =
                  element_line(size = base_linesize, linetype = 3, color = 'grey80'))
        },
        if (identical(minor_grid, 'x')) {
          theme(panel.grid.minor.x =
                  element_line(size = base_linesize, linetype = 3, color = 'grey80'))
        },
        if (identical(minor_grid, 'xy') | identical(grid, 'yx')) {
          theme(panel.grid.minor.y =
                  element_line(size = base_linesize, linetype = 3, color = 'grey80'),
                panel.grid.minor.x =
                  element_line(size = base_linesize, linetype = 3, color = 'grey80'))
        },
        if (isTRUE(panel_border)) {
          theme(
            panel.border =
              element_rect(fill = NA)
          )
        },
        if (!isTRUE(show_legend)) {
          theme(legend.position = 'none')
        },
        if (axis == 'x') {
          theme(
            axis.line.x = element_line(linetype = 1, color = 'black')
          )
        },
        if (axis == 'y') {
          theme(
            axis.line.y = element_line(linetype = 1, color = 'black')
          )
        },
        if (axis == 'xy') {
          theme(
            axis.line = element_line(linetype = 1, color = 'black')
          )
        },
        if (!is.na(ar)) {
          theme(
            aspect.ratio = ar
          )
        }
      )
    }
  
})


# Global functions figures ----------------------------------------


#' Export ggplot
#' 
#' @author Jonas Schöley
ExportFigure <-
  function(figure,
           path,
           filename,
           width = 170,
           height = 100,
           scale = 1,
           device = 'png',
           dpi = 300,
           add_date = FALSE) {
    require(ggplot2)
    
    if (missing(filename)) {
      filename <- tolower(gsub('\\.', '_', make.names(deparse(substitute(figure)))))
    }
    if (isTRUE(add_date)) {
      filename <- paste0(Sys.Date(), '-', filename)
    }
    
    arguments <-
      list(
        filename = paste0(filename, '.', device),
        plot = figure,
        path = path,
        width = width,
        height = height,
        units = "mm",
        scale = scale,
        dpi = dpi,
        device = device
      )
    if (device == 'pdf') {
      arguments$useDingbats <- FALSE 
    }
    
    do.call(ggsave, arguments)
  }

#' Export ggplots Stored in List
#' 
#' @author Jonas Schöley
ExportFiguresFromList <- function(lst, path, ...) {
  figure_names <- tolower(gsub('\\.+', '_', make.names(names(lst))))
  Fun <- function (figure, filename, ...) {
    ExportFigure(figure = figure, filename = filename, ...)
  }
  purrr::pwalk(
    list(lst, figure_names),
    Fun, path = path, ...
  )
}

# Global functions date -------------------------------------------

#' Create Unique Row ID
#'
#' @param region_iso iso-3166-1 alpha 2 country code with optional
#'   iso-3166-2 region code added, separated by a hyphen.
#' @param sex 'Male' or 'Female'
#' @param age_start Positive Integer.
#' @param year Positive Integer.
#' @param week Positive Integer.
#'
#' @return
#' String with fixed length row ID constructed from input.
#'
#' @examples
#' GenerateRowID('DE-BW', 'Male', 0, 2020, 10)
GenerateRowID <- function(region_iso, sex, age_start, year, week) {
  region_id <- sapply(region_iso, function (x) {
    expanded_region <- '------'
    substr(expanded_region, 1, nchar(x)) <- x
    return(expanded_region)
  })
  sex_id <- as.character(factor(sex, c('Male', 'Female'), c('M', 'F')))
  age_id <- sprintf('%02d', age_start)
  year_id <- sprintf('%04d', year)
  week_id <- sprintf('%02d', week)
  
  row_id <- paste0(region_id, sex_id, age_id, year_id, week_id)
  
  return(row_id)
}

#' Calculate Weeks Since Some Origin Date
#'
#' @param date Date string.
#' @param origin_date Date string.
#' @param week_format Either 'integer' for completed weeks or
#' 'fractional' for completed fractional weeks.
#'
#' @return Time difference in weeks.
#'
#' @author Jonas Schöley
#'
#' @examples
#' # My age in completed weeks
#' WeeksSinceOrigin(Sys.Date(), '1987-07-03')
WeeksSinceOrigin <-
  function(date, origin_date, week_format = "integer") {
    require(ISOweek)
    fractional_weeks_since_origin <-
      as.double(difftime(
        as.Date(date),
        as.Date(origin_date),
        units = "weeks"
      ))
    switch(
      week_format,
      fractional = fractional_weeks_since_origin,
      integer = as.integer(fractional_weeks_since_origin)
    )
  }

#' Convert Week of Year to Date
#'
#' @param year Year integer.
#' @param week Week of year integer (1 to 53).
#' @param weekday Weekday integer (1, Monday to 7, Sunday).
#' @param offset Integer offset added to `week` before date calculation.
#'
#' @return A date object.
#' 
#' @source https://en.wikipedia.org/wiki/ISO_8601
#'
#' @author Jonas Schöley
#'
#' @examples
#' # the first Week of 2020 actually starts Monday, December 30th 2019
#' ISOWeekDate2Date(2020, 1, 1)
ISOWeekDate2Date <- function (year, week, weekday = 1, offset = 0) {
  require(ISOweek)
  isoweek_string <-
    paste0(
      year, '-W',
      formatC(
        week+offset,
        flag = '0',
        format = 'd',
        digits = 1
      ),
      '-', weekday
    )
  ISOweek2date(isoweek_string)
}

#' Convert ISO-weeks to Epi-weeks
#' 
#' @param iso_week ISO-week vector.
#' @param w_start ISO-week for which Epi-week is 0.
#' @param w_53 Assume a 53 week year (default = FALSE)? 
#' 
#' @return A vector of Epi-weeks corresponding to the ISO-week input.
#' 
#' @author Jonas Schöley
#' 
#' @examples
#' # Epi-week 0 aligned with ISO-week 1
#' IsoWeekToEpiWeek(1:52)
#' # Epi-week 0 aligned with ISO-week 10
#' IsoWeekToEpiWeek(1:52, w_start = 10)
#' # Epi-week 0 aligned with ISO-week 10; 53 week year
#' IsoWeekToEpiWeek(1:52, w_start = 10, w_53 = TRUE)
ISOWeekToEpiWeek <- function (iso_week, w_start = 1, w_53 = FALSE) {
  a <- iso_week - w_start
  ifelse(a < 0, 51 + a + 1 + w_53, a)
}

#' Convert Epi-weeks to ISO-weeks
#' 
#' @param epi_week Epi-week vector.
#' @param w_start ISO-week for which Epi-week is 0.
#' @param w_53 Assume a 53 week year (default = FALSE)? 
#' 
#' @return A vector of ISO-weeks corresponding to the Epi-week input.
#' 
#' @author Jonas Schöley
#' 
#' @examples
#' epi_weeks = 0:51
#' # convert to iso week
#' iso_weeks <- EpiWeekToIsoWeek(epi_weeks, w_start = 10)
#' # convert back to epi week
#' IsoWeekToEpiWeek(iso_weeks, w_start = 10)
EpiWeekToIsoWeek <- function (epi_week, w_start = 1, w_53 = FALSE) {
  a <- epi_week + w_start
  ifelse(a > (52 + w_53), a - (52 + w_53), a)
}

# Global functions excess deaths ----------------------------------

#' Calculate Excess Deaths
#'
#' Calculate excess deaths and cumulative excess deaths with prediction
#' intervals from weekly observed death counts and simulated
#' replications of predicted weekly death counts.
#' 
#' @param df A data frame.
#' @param date Name of date variable.
#' @param observed_deaths Name of observed weekly death count variable.
#' @param expected_deaths Name of expected weekly death count variable.
#' @param simulated_deaths Name of predicted/simulated weekly death
#'                         count variable based on expected deaths.
#' @param sim_id Name of variable identifying simulation series.
#' @param ... Names of stratum variables. Death counts in strata not
#'            mentioned here will be summed across.
#'
#' @return A data frame with cumulative excess deaths over time with
#'         prediction intervals.
ExcessDeaths <-
  function (df, date, observed_deaths,
            expected_deaths, simulated_deaths, sim_id, ...) {
    
    require(dplyr)
    
    .date <- enquo(date)
    .observed_deaths <- enquo(observed_deaths)
    .expected_deaths <- enquo(expected_deaths)
    .simulated_deaths <- enquo(simulated_deaths)
    .sim_id <- enquo(sim_id)
    .strata <- enquos(...)
    
    # aggregate to time series of deaths to desired strata
    death_replicates_by_stratum <-
      df %>%
      group_by(!!.sim_id, !!!.strata, !!.date) %>%
      summarise(
        deaths_observed = sum(!!.observed_deaths),
        deaths_expected = sum(!!.expected_deaths),
        deaths_simulated = sum(!!.simulated_deaths)
      ) %>%
      ungroup()
    
    Excess <- function (observed, baseline) observed - baseline
    Pscore <- function (excess, baseline) excess / baseline
    
    observations_expectations_by_stratum <-
      death_replicates_by_stratum %>%
      filter(!!.sim_id == first(!!.sim_id)) %>%
      group_by(!!!.strata) %>% arrange(!!.date) %>%
      transmute(
        !!.date := !!.date,
        
        deaths_observed = deaths_observed,
        cumdeaths_observed = cumsum(deaths_observed),
        
        baseline_expected = deaths_expected,
        cumbaseline_expected = cumsum(baseline_expected),
        excess_expected = Excess(deaths_observed, baseline_expected),
        cumexcess_expected = cumsum(excess_expected),
        pscore_expected = Pscore(excess_expected, baseline_expected),
        cumpscore_expected = Pscore(cumexcess_expected, cumbaseline_expected)
      ) %>%
      ungroup()
    
    Q05 <- function (x) quantile(x, probs = 0.05, na.rm = TRUE)
    Q95 <- function (x) quantile(x, probs = 0.95, na.rm = TRUE)
    Avg <- function (x) mean(x, na.rm = TRUE)
    
    quantiles_by_stratum <-
      death_replicates_by_stratum %>%
      group_by(!!.sim_id, !!!.strata) %>%
      arrange(!!.date) %>%
      transmute(
        !!.date := !!.date,
        baseline = deaths_simulated,
        cumbaseline = cumsum(baseline),
        excess = Excess(deaths_observed, baseline),
        cumexcess = cumsum(excess),
        pscore = Pscore(excess, baseline),
        cumpscore = Pscore(cumexcess, cumbaseline)
      ) %>%
      # summarise across sim draws
      group_by(!!!.strata, !!.date) %>%
      summarise(
        across(
          c(baseline, cumbaseline, excess, cumexcess, pscore, cumpscore),
          list(q05 = Q05, q95 = Q95, avg = Avg),
          .names = '{.col}_{.fn}'
        )
      ) %>%
      ungroup()
    
    join_vars <- c(unlist(purrr::map(.strata, rlang::as_name)), rlang::as_name(.date))
    left_join(observations_expectations_by_stratum, quantiles_by_stratum, by = join_vars) %>%
      arrange(!!!.strata, !!.date) %>%
      select(all_of(join_vars), deaths_observed, cumdeaths_observed,
             starts_with('baseline'), starts_with('cumbaseline'),
             starts_with('excess'), starts_with('cumexcess'),
             starts_with('pscore'), starts_with('cumpscore')
      )
    
  }

#' Calculate Excess Deaths Residuals and Residual Summaries
#'
#' Calculate excess deaths and cumulative excess deaths with prediction
#' intervals from weekly observed death counts and simulated
#' replications of predicted weekly death counts.
#' 
#' @param df A data frame.
#' @param date Name of date variable.
#' @param observed_deaths Name of observed weekly death count variable.
#' @param predicted_deaths Name of predicted weekly death count variable.
#' @param cv_id Name of variable identifying cross-validation series.
#' @param ... Names of stratum variables. Death counts in strata not
#'            mentioned here will be summed across.
#'
#' @return A data frame with cumulative excess deaths over time with
#'         prediction intervals.
DeathsResiduals <-
  function (df, date, observed_deaths, predicted_deaths, cv_id, ...) {
    
    .date <- enquo(date)
    .observed_deaths <- enquo(observed_deaths)
    .predicted_deaths <- enquo(predicted_deaths)
    .cv_id <- enquo(cv_id)
    .strata <- enquos(...)
    
    # calculate raw residuals of weekly deaths and cumulative deaths
    # after aggregation of death counts over into specified strata
    residual_raw <-
      df %>%
      # aggregate to desired population
      group_by(!!.cv_id, !!!.strata, !!.date) %>%
      summarise(
        observed_deaths = sum(!!.observed_deaths),
        predicted_deaths = sum(!!.predicted_deaths)
      ) %>%
      # calculate residuals in weekly deaths and cumulative
      # deaths by cv_id and stratum
      group_by(!!.cv_id, !!!.strata) %>%
      arrange(!!.date) %>%
      mutate(
        resid_e_death = observed_deaths - predicted_deaths,
        resid_pe_death = resid_e_death / observed_deaths * 100,
        observed_cumdeath = cumsum(observed_deaths),
        predicted_cumdeath = cumsum(predicted_deaths),
        resid_e_cumdeath = observed_cumdeath - predicted_cumdeath,
        resid_pe_cumdeath = resid_e_cumdeath / observed_cumdeath * 100
      )
    
    # summarise weekly and cumulative residuals into measures of
    # error, bias, and variance
    residual_summary <-
      residual_raw %>%
      # summarise residuals across cv_id
      group_by(!!!.strata, !!.date) %>%
      summarise(
        ### CUMULATIVE DEATHS
        # raw error
        me_cumdeath = mean(resid_e_cumdeath),
        se_cumdeath = sd(resid_e_cumdeath),
        # percentage error
        mpe_cumdeath = mean(resid_pe_cumdeath),
        spe_cumdeath = sd(resid_pe_cumdeath),
        # absolute error
        mae_cumdeath = mean(abs(resid_e_cumdeath)),
        sae_cumdeath = sd(abs(resid_e_cumdeath)),
        # absolute percentage error
        mape_cumdeath = mean(abs(resid_pe_cumdeath)),
        sape_cumdeath = sd(abs(resid_pe_cumdeath)),
        
        ### WEEKLY DEATHS
        # raw error
        me_death = mean(resid_e_death),
        se_death = sd(resid_e_death),
        # percentage error
        mpe_death = mean(resid_pe_death),
        spe_death = sd(resid_pe_death),
        # absolute error
        mae_death = mean(abs(resid_e_death)),
        sae_death = sd(abs(resid_e_death)),
        # absolute percentage error
        mape_death = mean(abs(resid_pe_death)),
        sape_death = sd(abs(resid_pe_death))
      ) %>%
      ungroup()
    
    return(list(residual_raw = residual_raw, residual_summary = residual_summary))
    
  }
