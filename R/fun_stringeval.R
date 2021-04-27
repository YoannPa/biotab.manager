##FUNCTIONS

#' Checks if an alphanumerical string contains any upper case character.
#'
#' @param str A \code{character} vector to be evaluated looking for any
#'                               uppercase character.
#' @return A \code{logical}.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

is.upper <- function(str){ grepl(pattern = "[A-Z]+", x = str) }

#' Checks if an array, a dataframe or a matrix contains any empty cell.
#'
#' @param data An \code{array}, a \code{dataframe} or a \code{matrix}.
#' @return A \code{logical}, TRUE if the array contains at least 1 empty cell,
#'         FALSE if it does not.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

any.empty <- function(data){
  any(apply(X = data, MARGIN = 2,FUN = grepl, pattern = "^$"))
}

#' Develops a data.table row if it contains information for more than 1 sample.
#'
#' @param dt.row A \code{data.table} of 1 single row.
#' @return A developped \code{data.table} of n rows, with 1 row per sample.
#' @author Yoann Pageaud.
#' @export

dt.develop <- function(dt.row){
  #Develop only rows with num_samples > 1
  if(dt.row$num_samples > 1){
    #Develop all columns from 50 to 76
    #Split following the "||" separator
    res.split.pipe <- apply(
      X = dt.row[, -c(1:48), ], MARGIN = 2,
      FUN = function(i){
        if(is.character(i)){ strsplit(x = i, split = "\\|\\|") } else { i }
      })
    #Unlist each result within the row
    res.split.pipe <- lapply(X = res.split.pipe, FUN = unlist)
    #Split following the "," separator
    res.split.comma <- lapply(X = res.split.pipe, FUN = strsplit, split = ",")
    res.split.comma <- lapply(X = res.split.comma, FUN = unlist)
    #Merge back values together in columns "wgs_rna_icgc_specimen_id" and
    # "wgs_rna_icgc_specimen_type" if there length is above 1
    if(length(res.split.comma$wgs_rna_icgc_specimen_id) > 1){
      res.split.comma$wgs_rna_icgc_specimen_id <- paste(
        res.split.comma$wgs_rna_icgc_specimen_id, collapse = " | ")
    }
    if(length(res.split.comma$wgs_rna_icgc_specimen_type) > 1){
      res.split.comma$wgs_rna_icgc_specimen_type <- paste(
        res.split.comma$wgs_rna_icgc_specimen_type, collapse = " | ")
    }
    #Compute length of each element
    res.length <- unlist(lapply(X = res.split.comma, FUN = length))

    #Check that every element is either of length 1 or length = num_samples
    if(all(res.length == 1 | res.length == dt.row$num_samples)){
      #Create new data.table
      new.half.dt <- as.data.table(res.split.comma)
      #Replicate num_samples times the first half of the row
      first.half.dt <- dt.row[rep(seq(nrow(dt.row)), num_samples)][, c(1:48), ]
      #Cbind first and second half
      new.dtrow <- cbind(first.half.dt, new.half.dt)
    } else { stop("Some elements length do not match requirements.") }
  } else { new.dtrow <- dt.row }
  #Return developped data.table
  return(new.dtrow)
}
