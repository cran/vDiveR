% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_dynamics_proteome.R
\name{plot_dynamics_proteome}
\alias{plot_dynamics_proteome}
\title{Dynamics of Diversity Motifs (Proteome) Plot}
\usage{
plot_dynamics_proteome(df, host = 1, dot_size = 2, word_size = 15, alpha = 1/3)
}
\arguments{
\item{df}{DiMA JSON converted csv file data}

\item{host}{number of host (1/2)}

\item{dot_size}{size of dot in plot}

\item{word_size}{word size in plot}

\item{alpha}{any number from 0 (transparent) to 1 (opaque)}
}
\value{
A plot
}
\description{
This function compactly display the dynamics of diversity motifs (index and its variants: major, minor and unique)
in the form of dot plot as well as violin plot for all the provided proteins at proteome level.
}
\examples{
plot_dynamics_proteome(proteins_1host)
}
