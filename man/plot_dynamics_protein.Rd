% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_dynamics_protein.R
\name{plot_dynamics_protein}
\alias{plot_dynamics_protein}
\title{Dynamics of Diversity Motifs (Protein) Plot}
\usage{
plot_dynamics_protein(
  df,
  host = 1,
  protein_order = "",
  base_size = 8,
  alpha = 1/3,
  line_dot_size = 3,
  bw = "nrd0",
  adjust = 1
)
}
\arguments{
\item{df}{DiMA JSON converted csv file data}

\item{host}{number of host (1/2)}

\item{protein_order}{order of proteins displayed in plot}

\item{base_size}{base font size in plot}

\item{alpha}{any number from 0 (transparent) to 1 (opaque)}

\item{line_dot_size}{dot size in scatter plot}

\item{bw}{smoothing bandwidth of violin plot (default: nrd0)}

\item{adjust}{adjust the width of violin plot (default: 1)}
}
\value{
A plot
}
\description{
This function compactly display the dynamics of diversity motifs (index and its variants: major, minor and unique)
in the form of dot plot(s) as well as violin plots for all the provided individual protein(s).
}
\examples{
plot_dynamics_protein(proteins_1host)
}
