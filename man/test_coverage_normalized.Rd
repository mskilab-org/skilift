% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{test_coverage_normalized}
\alias{test_coverage_normalized}
\title{Test coverage values}
\usage{
test_coverage_normalized(
  coverage_values,
  tolerance = 0.1,
  fraction_near_1 = 0.9
)
}
\description{
Test for mean normalized coverage values
}
\details{
Dryclean, fragcounter, cobalt, or some other
binned coverage counter may emit coverage values in  
units of reads or mean normalized, gc- and/or mappability-bias
corrected coverage units.
}
