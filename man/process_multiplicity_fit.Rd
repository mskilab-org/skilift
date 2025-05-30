% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segment-width-distribution.R
\name{process_multiplicity_fit}
\alias{process_multiplicity_fit}
\title{process_multiplicity_fit}
\usage{
process_multiplicity_fit(
  variants,
  field = "altered_copies",
  mask = TRUE,
  mask_gr = system.file("extdata", "data", "maskA_re.rds", package = "Skilift"),
  save_data = TRUE,
  save_png = TRUE,
  save_html = FALSE,
  bins = 1e+06,
  out_file
)
}
\arguments{
\item{variants}{GRanges object containing the variants}

\item{field}{field to use from the input data; default is altered_copies. other often used paramaters are ["total_snv_copies","altered_copies", "major_snv_copies", "minor_snv_copies"]}

\item{mask}{logical value to mask the data or not; default is TRUE}

\item{mask_gr}{GRanges object containing the mask; default is maskA as provided in the package}

\item{bins}{number of bins for histogram; should specify for lower limit to avoid performance issues; default = 100000}

\item{out_file}{output file path}
}
\description{
a method to lift the mutations as a histogram to JSON for viewing in gOS
}
