% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_entropy.R
\name{plot_entropy}
\alias{plot_entropy}
\title{Entropy plot}
\usage{
plot_entropy(
  df,
  host = 1,
  protein_order = "",
  kmer_size = 9,
  ymax = 10,
  line_dot_size = 2,
  base_size = 8,
  all = TRUE,
  highlight_zero_entropy = TRUE
)
}
\arguments{
\item{df}{DiMA JSON converted csv file data}

\item{host}{number of host (1/2)}

\item{protein_order}{order of proteins displayed in plot}

\item{kmer_size}{size of the k-mer window}

\item{ymax}{maximum y-axis}

\item{line_dot_size}{size of the line and dot in plot}

\item{base_size}{word size in plot}

\item{all}{plot both the entropy and total variants (pass FALSE in to plot only the entropy)}

\item{highlight_zero_entropy}{highlight region with zero entropy (default: TRUE)}
}
\value{
A plot
}
\description{
This function plot entropy (black) and total variant (red) incidence of each
k-mer position across the studied proteins and highlight region(s) with zero entropy in yellow.
k-mer position with low support is marked with a red triangle underneath the x-axis line.
}
\examples{
plot_entropy(proteins_1host)
plot_entropy(protein_2hosts, host = 2)
}
