% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metadata.R
\name{add_coverage_parameters}
\alias{add_coverage_parameters}
\title{Add Coverage Parameters}
\usage{
add_coverage_parameters(metadata, tumor_coverage, field = "foreground")
}
\arguments{
\item{metadata}{A data.table containing metadata}

\item{tumor_coverage}{Coverage data for the tumor}

\item{field}{The field to use for coverage calculation}
}
\value{
Updated metadata with coverage parameters added
}
\description{
Adds coverage parameters to the metadata based on tumor coverage
}
