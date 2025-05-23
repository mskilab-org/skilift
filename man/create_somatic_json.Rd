% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{create_somatic_json}
\alias{create_somatic_json}
\alias{dt2json_mut}
\title{create_somatic_json}
\usage{
create_somatic_json(plot_metadata, datadir, settings = internal_settings_path)

dt2json_mut(
  dt,
  patient.id,
  ref,
  settings,
  file_name = paste(getwd(), "test.json", sep = "/"),
  meta_data = NULL,
  y_col = NULL
)
}
\arguments{
\item{plot_metadata}{data.table with Plot metadata, columns = (patient.id, source, x (contains path to data file or reference to data file object itself), ref, overwrite).}

\item{datadir}{Path to data directory with patient directories.}

\item{settings}{Path to settings.json file.}

\item{dt}{data.table with seqnames,start and end}

\item{patient.id}{patient id to be added to pgvdb or case reports}

\item{ref}{reference for pgv or case reports}

\item{file_name}{the file the json should be written to}

\item{meta_data}{extra information to add to the json}

\item{y_col}{column to be used for y axis values}

\item{somatic_snv_cn}{somatic copy number data table}

\item{out_file}{where the json should be written}

\item{pair}{patient.id to be used for pgvdb or case reports}

\item{pgv_settings}{the settings file used for pgvdb or case reports}

\item{return_table_pgv}{TRUE/FALSE whether to return the data.table to add to pgvdb}

\item{meta_keep}{meta data to write in the json}
}
\value{
NULL.

data.table to add to pgv or NULL, depends on return_table_pgv
}
\description{
Create mutations gGraph JSON file for SOMATIC MUTATIONS.

function to create a mutation json, used in create_somatic_json
}
\author{
Stanley Clarke
}
