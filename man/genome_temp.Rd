% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{genome_temp}
\alias{genome_temp}
\title{genome_temp}
\usage{
genome_temp(
  patient_id = NA,
  order = NA,
  x = list(NA),
  ref = NA,
  source = "genome.json",
  type = "genome",
  visible = TRUE,
  title = NA,
  max.cn = NULL,
  annotation = list(c("bfb", "chromoplexy", "chromothripsis", "del", "dm", "cpxdm",
    "dup", "pyrgo", "rigma", "simple", "tic", "tyfonas")),
  overwrite = FALSE
)
}
\arguments{
\item{patient_id}{patient.id to be added to pgvdb}

\item{order}{optional entry if you order plots with a column order}

\item{x}{path to a JaBbA ggraph object or object itself as a list}

\item{ref}{reference to use for pgvdb}

\item{type}{genome, can be changed to allelic to render allelic graphs}

\item{visible}{TRUE/FALSE whether the plot is hidden or showing in pgv}

\item{title}{title of the plot in pgvdb}

\item{max.cn}{override max cn of 100}

\item{annotation}{default is list of SVs, make null if no annotations present in object}

\item{overwrite}{TRUE/FALSE to overwrite an existing genome json}
}
\description{
function to create data.table for genome graphs for pgvdb
}
\author{
Stanley Clarke
}
