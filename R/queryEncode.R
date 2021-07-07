#' query data from ENCODE by predefined criteria
#' @description Search ENCODE data by querying the ENCODE Portal using its
#' REST API.
#' @importFrom httr GET add_headers http_error http_status http_type content
#' @importFrom jsonlite fromJSON
#' @export
#' @param API_url ENCODE REST API url.
#' @param exactMatch character. Exact-match keywords refer to search results
#' that perfectly match all the keywords in the search query, exactly as
#' entered. It is a named character vector. The names will be the keys
#' and characters will be the values for search.
#' @param partialMatch character. Partial-match refer to search results that
#' contain the search query. It is a named character vector.
#' The names will be the keys and characters will be the values for search.
#' The value starting from '!' indicates logical negation(NOT).
#' The value starting from '>', '>=', '<', '==', '<=' indicates string
#' comparison.
#' @param API_url character. The ENCODE REST API url.
#' @param ... Not used.
#' @return A list of search results
#' @examples
#' res <- queryEncode(
#'  exactMatch=c(
#'   target.label="H3K4me1",
#'   replicates.library.biosample.donor.organism.scientific_name="Homo sapiens",
#'   assembly="GRCh38",
#'   assay_term_name="ChIP-seq"),
#'   partialMatch=c(biosample_summary="heart"))
queryEncode <- function(exactMatch,
                        partialMatch=character(0),
                        API_url="https://www.encodeproject.org/search/",
                        ...){
  stopifnot("exactMatch must be a named characters"=
              length(exactMatch)>0 &&
              length(names(exactMatch))==length(exactMatch) &&
              is.character(exactMatch))
  if(length(partialMatch)>0){
    stopifnot("partialMatch must be a named characters"=
                length(names(partialMatch))==length(partialMatch) &&
                is.character(partialMatch))
  }
  fixed_keys <- c(status="released",
                  limit="all",
                  frame="object",
                  format="json")
  fixed_keys <- fixed_keys[!names(fixed_keys) %in% names(exactMatch)]
  header <- c(Accept='application/json')
  response <- GET(url=API_url,
                  query = as.list(c(fixed_keys, exactMatch)),
                  config = add_headers(header))
  e <- http_error(response)
  if(e){
    message("searching string:", content(response, encodeing="UTF-8")$`@id`)
    stop(http_status(response)$message)
  }
  if (http_type(response) != "application/json") {
    stop("API did not return json", call. = FALSE)
  }
  cont <- fromJSON(httr::content(response, type="text", encoding="UTF-8"),
                   simplifyVector = FALSE)
  res <- cont$`@graph`

  if(length(partialMatch)){
    searchInList <- function(.ele, value){
      if(length(.ele)==0) return(TRUE) # the key is not in returned list.
      if(is.list(.ele) || length(.ele)>1){
        any(vapply(.ele, FUN= searchInList, FUN.VALUE = logical(1),
                   value=value))
      }else{
        if(substr(value, 1, 1) %in% c('<', '=', '>')){
          if(substr(value, 1, 2) %in% c('<=', '==', '>=')){
            switch (substr(value, 1, 2),
                    '<=' = .ele <= substring(value, 3),
                    '==' = .ele == substring(value, 3),
                    '>=' = .ele >= substring(value, 3)
            )
          }else{
            switch (substr(value, 1, 1),
                    '<' = .ele < substring(value, 2),
                    '=' = .ele == substring(value, 2),
                    '>' = .ele > substring(value, 2)
            )
          }
        }else{
          if(grepl("^!", value)){
            value <- sub("^.", "", value)
            !grepl(tolower(value), tolower(.ele))
          }else{
            grepl(tolower(value), tolower(.ele))
          }
        }
      }
    }
    keep <- mapply(FUN = function(key, value){
      vapply(lapply(res, `[[`, i=key),
             FUN = searchInList, FUN.VALUE = logical(1),
             value=value)
    }, names(partialMatch), partialMatch, SIMPLIFY = FALSE)
    keep <- do.call(cbind, keep)
    keep <- apply(keep, 1, sum) >0
    res <- res[keep]
  }

  return(res)
}
