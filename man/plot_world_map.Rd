% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_world_map.R
\name{plot_world_map}
\alias{plot_world_map}
\title{Geographical Distribution of Sequences Plot}
\usage{
plot_world_map(metadata, base_size = 8)
}
\arguments{
\item{metadata}{a dataframe with 3 columns, 'ID', 'region', and 'date'}

\item{base_size}{word size in plot}
}
\value{
A list with 2 elements (a plot followed by a dataframe)
}
\description{
This function plots a world map and color the affected geographical region(s)
from light (lower) to dark (higher), depends on the cumulative number of sequences.
Aside from the plot, this function also returns a dataframe with 2 columns: 'Region' and 'Number of Sequences'.
The input dataframe of this function is obtainable from metadata_extraction(), with NCBI
Protein / GISAID (EpiFlu/EpiCoV/EpiPox/EpiArbo) FASTA file as input.
}
\examples{
geographical_plot <- plot_world_map(metadata)$plot
geographical_df <- plot_world_map(metadata)$df
}
