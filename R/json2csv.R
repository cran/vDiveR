#' JSON2CSV
#'
#' This function converts DiMA (v5.0.9) JSON output file to a dataframe with 17
#' predefined columns which further acts as the input for other functions provided in this vDiveR package.
#'
#' @param json_data DiMA JSON output dataframe
#' @param host_name name of the host species
#' @param protein_name name of the protein
#' @return A dataframe which acts as input for the other functions in vDiveR package
#' @examples inputdf<-json2csv(JSON_sample)
#' @importFrom stats aggregate
#' @importFrom dplyr mutate_if right_join distinct
#' @importFrom tidyr replace_na
#' @export
json2csv <-function(json_data, host_name="unknown host", protein_name="unknown protein"){
    # Declare global variables to avoid R CMD check warnings
    Group.2 <- x <- results.position <- results.diversity_motifs <- motif_short <- NULL
    position <- entropy <- support <- low_support <- distinct_variants_incidence <- NULL
    total_variants_incidence <- Ma <- Mi <- U <- I <- proteinName <- NULL
    results.support <- results.low_support <- results.entropy <- NULL
    results.total_variants_incidence <- results.distinct_variants_incidence <- NULL
    average_entropy <- highest_entropy.position <- highest_entropy.entropy <- host <- NULL
    sequence <- multiIndex <- NULL

    data_flatten <- as.data.frame(json_data) %>%
        tidyr::unnest(cols = c(results.diversity_motifs))

    # data transformation
    motifs_incidence <- aggregate(
        data_flatten$incidence,
        list(data_flatten$results.position,data_flatten$motif_short),
        FUN=sum
    ) %>%
        tidyr::spread(Group.2,x) %>%                           # transpose the rows (motif-long & incidence) to columns (index, major, minor, unique)
        dplyr::mutate_if(is.numeric, ~(tidyr::replace_na(., 0)))      # replace NAN with 0
    #rename column 'Group.1' to 'results.position'


    colnames(motifs_incidence)[colnames(motifs_incidence) == 'Group.1'] <- "results.position"

    #sum the number of Index motif found in each position (if > 1 => multiIndex == TRUE)
    multiIndex<-data_flatten%>%
        dplyr::group_by(results.position) %>%
        dplyr::summarize(multiIndex = sum(motif_short=="I"), .groups = "drop") %>%
        as.data.frame()

    #merge multiIndex to motifs_incidence df
    motifs_incidence <-dplyr::left_join(motifs_incidence,multiIndex, by='results.position')%>%
        distinct()

    #replace multiIndex with boolean: x > 1 = TRUE and vice versa
    motifs_incidence$multiIndex <- ifelse(
        tidyr::replace_na(motifs_incidence$multiIndex, 0) > 1,
        TRUE,
        FALSE
    )

    # Ensure expected incidence columns exist even if absent in data
    for (col in c("I", "Ma", "Mi", "U")) {
        if (!col %in% colnames(motifs_incidence)) {
            motifs_incidence[[col]] <- 0
        }
    }

    # ----- Extract index sequence only (for positions that have Index motif)
    # extract the first encountered index kmer sequence if > 1 index is encountered for each position (rarely happens)
    index_sequences<-subset(data_flatten, motif_short == "I") %>% #extract rows that are of index
        dplyr::group_by(results.position) %>% #group them based on position
        dplyr::slice(1) %>% #take the first index encountered per position
        dplyr::ungroup() %>%
        dplyr::select(results.position, sequence) %>%
        as.data.frame() #return the data in dataframe

    # ------ Extract position-level fields from results (not motif-specific)
    results_df <- as.data.frame(json_data$results)

    results_keep <- results_df %>%
        dplyr::transmute(
            results.position = position,
            results.entropy = entropy,
            results.support = support,
            results.low_support = low_support,
            results.distinct_variants_incidence = distinct_variants_incidence,
            results.total_variants_incidence = total_variants_incidence
        )
    
    #replace low_support of NA to FALSE
    results_keep$results.low_support[is.na(results_keep$results.low_support)] <- FALSE

    # ---- Global highest/average entropy (same for all rows) ----
    highest_pos <- if (!is.null(json_data$highest_entropy$position)) json_data$highest_entropy$position else NA
    highest_ent <- if (!is.null(json_data$highest_entropy$entropy)) json_data$highest_entropy$entropy else NA
    avg_ent <- if (!is.null(json_data$average_entropy)) json_data$average_entropy else NA

    results_keep$highest_entropy.position <- highest_pos
    results_keep$highest_entropy.entropy <- highest_ent
    results_keep$average_entropy <- avg_ent

    # ---- all positions in results
    motifs <- dplyr::left_join(results_keep, motifs_incidence, by = "results.position")
    # If any positions had no motifs rows at all (rare), fill incidence NAs with 0
    motifs <- motifs %>%
        dplyr::mutate(
            I = tidyr::replace_na(I, 0),
            Ma = tidyr::replace_na(Ma, 0),
            Mi = tidyr::replace_na(Mi, 0),
            U = tidyr::replace_na(U, 0),
            multiIndex = tidyr::replace_na(multiIndex, FALSE)
        )

    # Join index sequence (will be NA for positions without Index motif)
    if (nrow(index_sequences) > 0) {
        motifs <- dplyr::left_join(
            motifs,
            index_sequences,
            by = "results.position"
        )
    } else {
        # Create empty column if no Index exists anywhere (edge-case)
        motifs$sequence <- NA
    }

    #assign host + protein name
    motifs['host'] <- host_name
    motifs['proteinName'] <- protein_name

    # ---- final rename
    motifs <- motifs %>%
        dplyr::transmute(
            proteinName = proteinName,
            position = results.position,
            count = results.support,              # keeping your original naming
            lowSupport = results.low_support,
            entropy = results.entropy,
            indexSequence = sequence,
            index.incidence = I,
            major.incidence = Ma,
            minor.incidence = Mi,
            unique.incidence = U,
            totalVariants.incidence = results.total_variants_incidence,
            distinctVariant.incidence = results.distinct_variants_incidence,
            multiIndex = multiIndex,
            averageEntropy = average_entropy,
            highestEntropy.position = highest_entropy.position,
            highestEntropy = highest_entropy.entropy,
            host = host
        )

    motifs

}

