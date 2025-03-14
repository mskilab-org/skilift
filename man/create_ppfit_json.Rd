% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R, R/segment-width-distribution.R
\name{create_ppfit_json}
\alias{create_ppfit_json}
\alias{get_segstats}
\title{create_ppfit_json}
\usage{
create_ppfit_json(
  balanced_gg,
  cov_path = NULL,
  path_obj = NULL,
  settings_json,
  ref = NULL,
  max_na = NULL,
  out_file = NULL,
  write_json = TRUE,
  overwrite = FALSE,
  return_table = FALSE,
  cores = 1
)

get_segstats(
  balanced_jabba_gg,
  tumor_coverage,
  coverage_field = "foreground",
  settings = Skilift:::default_settings_path,
  ref = "hg19",
  max_na = 0.9,
  cores = 1
)
}
\arguments{
\item{balanced_gg}{the balanced_gg ggraph or the path to it}

\item{path_obj}{path to flow directory, if null, must supply coverage path}

\item{settings_json}{path to settings.json (from PGV)}

\item{ref}{reference name}

\item{max_na}{max.na used with JaBbA}

\item{out_file}{location to write json}

\item{write_json}{TRUE/FALSE whether to write the json}

\item{overwrite}{TRUE/FALSE whether to overwrite existing file}

\item{return_table}{TRUE/FALSE whether to return the data}

\item{cores}{cores for JaBbA:::segstats}

\item{balanced_jabba_gg}{the balanced_gg ggraph or the path to it}

\item{tumor_coverage}{path to associated coverage, if null, will try to pull path from pairs table}

\item{coverage_field}{field in coverage file to use for segstats}

\item{settings}{path to settings.json (from PGV)}

\item{cov}{path to associated coverage, if null, will try to pull path from pairs table}
}
\value{
NULL or segstats table

segstats table
}
\description{
function to create segmentation plots in case reports

function to create segmentation plots in case reports
}
\author{
Stanley Clarke, Tanubrata Dey
drcln_cov = readRDS(thisp$decomposed_cov[i])

Stanley Clarke, Tanubrata Dey, Shihab Dider
drcln_cov = readRDS(thisp$decomposed_cov[i])
}
