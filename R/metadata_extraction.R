#' Metadata Extraction from NCBI/GISAID EpiCoV FASTA file
#'
#' This function retrieves metadata (ID, country, date) from the NCBI/GISAID EpiCoV FASTA file (default FASTA header expected).
#'
#' @param file_path path of fasta file
#' @param source the source of fasta file, either "ncbi" or "GISAID"
#' @return  A dataframe that has three columns consisting ID, collected country and collected date
#' @examples filepath <- system.file('extdata','GISAID_EpiCoV.faa', package = 'vDiveR')
#' @examples meta_gisaid <- metadata_extraction(filepath, 'GISAID')
#' @export
metadata_extraction <- function(file_path, source){
    if(source == 'ncbi'){
        meta <- extract_from_NCBI(file_path)
    }
    if(source == 'GISAID'){
        meta <- extract_from_GISAID(file_path)
    }
    return(meta)
}

#' Extract metadata via fasta file from GISAID
#'
#' This function get the metadata from each header of GISAID fasta file
#' @param file_path path of fasta file
#' @importFrom stringr str_extract
extract_from_GISAID <- function(file_path){
    heads <- c()
    lines <- readLines(file_path, warn=FALSE)
    for(line in lines){
        if(grepl('>', line)){
            line <- substr(line, 2, nchar(line))
            heads <- c(heads, line)
        }
    }
    IDs <- c(); countries <- c(); dates <- c()
    for(head in heads){
        ID <- str_extract(head, "EPI[^|]*"); IDs <- c(IDs, ID)

        country <- tryCatch(strsplit(head, '/', fixed=T)[[1]][2], error = function(e) "NA")
        countries <- c(countries, country)

        date <- str_extract(head, "[0-9]{4}-[0-9]{2}-[0-9]{2}"); dates <- c(dates, date)
    }
    return(data.frame('ID' = IDs, 'country' = countries, 'date' = dates))
}

#' Extract metadata via fasta file from ncbi
#'
#' This function get the metadata from each head of fasta file
#' @param file_path path of fasta file
#' @importFrom rentrez entrez_search entrez_fetch
extract_from_NCBI <- function(file_path){
    IDs <- c()
    lines <- readLines(file_path, warn=FALSE)
    for(line in lines){
        if(grepl('>', line)){
            ID <- strsplit(line, ' ')[[1]][1]
            ID <- substr(ID,2,nchar(ID))
            IDs <- c(IDs, ID)
        }
    }
    countries <- c(); dates <- c(); dropsample <- c(); keepsample <- c()
    for(ID in IDs){
        if(substr(ID,1,3) == 'pdb'){
            dropsample <- c(dropsample, ID)
            next
        }
        keepsample <- c(keepsample, ID)
        search_result <- entrez_search(db = "protein", term = ID, retmax = 1)
        accession <- search_result$ids[[1]]
        
        info <- entrez_fetch(db = "protein", id = accession, rettype = "gb", retmode = "text")
        info <- tryCatch(strsplit(info, '\n')[[1]], error = function(e) "")
        idx1 <- grep('/country=',  info); idx2 <- grep('/collection_date=',  info)
        info1 <- info[idx1]; info2 <- info[idx2]

        country <- tryCatch(strsplit(info1, '\\"')[[1]][2], error = function(e) "NA")
        if (grepl(':', country)) {
            country <- strsplit(country,':')[[1]][1]
        }

        date <- tryCatch(strsplit(info2, '\\"')[[1]][2], error = function(e) "NA")
        
        countries <- c(countries, country); dates <- c(dates, date)
    }
    tmp <- data.frame('ID' = keepsample, 'country' = countries, 'date' = dates)

    for(i in 1:nrow(tmp)){
        if(!is.na(as.Date(tmp$date[i],format='%Y-%m-%d'))){
            tt <- 1
        } else if(!is.na(as.Date(tmp$date[i],format='%d-%B-%Y'))){
            tmp$date[i] <- as.character(as.Date(tmp$date[i],format="%d-%B-%Y"))
        } else {
            dropsample <- c(dropsample, tmp$ID[i])
        }
    }
    warninfo <- paste(c('\nExcluded records:\n',dropsample, '\nCountry/date (%Y-%m-%d OR %d-%B-%Y format) are not provided.'), collapse = " ")
    warning(warninfo)
    tmp <- tmp[! tmp$ID %in% dropsample, ]
    return(tmp)
}

