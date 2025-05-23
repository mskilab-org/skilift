% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filtered-events.R
\name{process_gencode}
\alias{process_gencode}
\title{process_gencode}
\usage{
process_gencode(gencode = NULL, assembly = NULL, seqlevelsstyle = "NCBI")
}
\arguments{
\item{gencode}{path to gencode file. Gencode file must be either rds or some format accepted by rtracklayer::import (e.g. GTF) with just a single entry for each gene (so gencode entries for each gene are collapse to a single range). The input could be .gtf or .rds with GRanges object, or a GRanges object i.e. resulting from importing the (appropriate) GENCODE .gtf via rtracklayer, note: this input is only used in CNA to gene mapping.}
}
\value{
gencode_gr GRanges
}
\description{
Helper script to process gencode parameter
}
\author{
Marcin Imielinski
}
