install.packages(
  'INLA',
  repos = c(getOption('repos'), INLA = 'https://inla.r-inla-download.org/R/stable'),
  dep = TRUE
)
