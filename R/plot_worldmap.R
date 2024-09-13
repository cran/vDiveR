#' Geographical Distribution of Sequences Plot
#'
#' This function plots a worldmap and color the affected geographical region(s)
#' from light (lower) to dark (higher), depends on the cumulative number of sequences.
#' Aside from the plot, this function also returns a dataframe with 2 columns: 'Region' and 'Number of Sequences'.
#' The input dataframe of this function is obtainable from metadata_extraction(), with NCBI
#' Protein / GISAID EpiCoV FASTA file as input.

#' @param meta a dataframe with 3 columns, 'ID', 'region', and 'date'
#' @param base_size word size in plot
#'
#' @return A list with 2 elements (a plot followed by a dataframe)
#' @examples geographical_plot <- plot_worldmap(metadata)$plot
#' @examples geographical_df <- plot_worldmap(metadata)$df
#' @importFrom ggplot2 geom_polygon scale_fill_gradient map_data
#' @importFrom dplyr left_join %>% group_by ungroup slice 
#' @importFrom stringr str_to_title
#' @export
plot_worldmap <- function(meta, base_size=8){
    long <- lat <- group <- count <- city_ascii <- region.y <- region <- ID <- date <- NULL
    if (nrow(meta) < 1){
        error_msg <- paste("No records found in the metadata dataframe.")
        return(list(plot = NULL, df = error_msg))
    }

    meta$region <- stringr::str_to_title(meta$region)
    region_list <- data.frame(table(meta$region))
    colnames(region_list) <- c('region','count')
  
    #==================== data preparation section ========================#
    build_in_path <- system.file("extdata", "city_mapper.csv", package = "vDiveR")
    city2region <- utils::read.csv(build_in_path, stringsAsFactors = FALSE) 
    city2region_unique <- city2region %>%
        dplyr::group_by(city_ascii) %>%
        dplyr::slice(1) %>%  # if more than 1 match, keep the first occurrence
        dplyr::ungroup()

    # Perform the left join to match 'region' with 'city_ascii' from city_mapper
    meta <- meta %>%
        left_join(city2region_unique, by = c("region" = "city_ascii"))

    # Replace the 'region' column in meta with the matched 'region' from city_mapper
    meta <- meta %>%
        mutate(region = ifelse(!is.na(region.y), region.y, region)) %>%
        select(ID, region,date)

    #================= plotting and tabulating section =====================#
    world_map <- ggplot2::map_data("world")
    missing_region <- region_list$region[! region_list$region %in% world_map$region]

    if(length(missing_region) > 0){
    warning_info <- paste(c(as.character(missing_region), ' do not exist in the ggplot2 world map. Please refer ggplot2 world map region list for the names.'), collapse = ", ")
    warning(warning_info)
    }

    p <- ggplot(world_map, aes(x = long, y = lat, group = group)) +
        geom_polygon(fill="lightgray", colour = "#888888")
    pathogens.map <- dplyr::left_join(region_list, world_map, by = "region")

    p <- p + geom_polygon(data = pathogens.map, aes(fill = count), color = "#888888") +
        scale_fill_gradient(low = "#FFFFFF", high = "#E63F00", name = 'Number of Sequences') +
        theme(plot.background = element_rect(fill = "transparent", colour = NA),
              panel.border = element_blank(), panel.grid = element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
              legend.position = "right") +
        theme_classic(base_size=base_size)

    colnames(region_list) <- c('Region','Number of Sequences')

    return(list(plot=p, df=region_list))

}


