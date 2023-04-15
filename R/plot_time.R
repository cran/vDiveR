#' Time Distribution of Sequences Plot
#'
#' This function plots the time distribution of provided sequences in the form of bar plot with 'Month' as x-axis and 
#' 'sequence recorded' as y-axis. The input dataframe of this function is obtainable from metadata_extraction(), with NCBI 
#' Protein / GISAID EpiCoV FASTA file as input.
#'
#' @param meta a dataframe with 3 columns, 'ID', 'country', and 'date'
#' @return  A plot
#' 
#' @examples plot_time(metadata)
#' @importFrom ggplot2 scale_x_date
#' @importFrom ggplot2 geom_bar
#' @importFrom scales date_format
#' @export
plot_time <- function(meta){
    Month <- NULL
    meta$Month <- as.Date(cut(as.Date(meta$date, format = "%Y-%m-%d"), breaks = "month"))
    p <- ggplot(data = meta, aes(x = Month)) + geom_bar() + ylab('Number of protein sequence records') +
        scale_x_date(date_breaks = "2 month", labels = date_format("%Y-%b"))  + theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = .5))
    return(p)
}
