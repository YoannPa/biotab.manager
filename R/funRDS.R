##FUNCTIONS

#' @description Save a list of named object as independent RDS files.
#' 
#' @param list.object A \code{list} of named object to serialize. Object names
#'                    are used as file names.
#' @param dir.save    A \code{character} to specify the path to a folder where
#'                    the RDS files should be saved.
#' @param ascii       A \code{logical}. If TRUE or NA, an ASCII representation
#'                    is written; otherwise (default), a binary one is used.
#' @param version     The workspace format version to use. NULL specifies the
#'                    current default version (3). The only other supported
#'                    value is 2, the default from R 1.4.0 to R 3.5.0.
#' @param compress    A logical specifying whether saving to a named file is to
#'                    use "gzip" compression, or one of "gzip", "bzip2" or "xz"
#'                    to indicate the type of compression to be used. Ignored if
#'                    file is a connection.
#' @param refhook     A hook function for handling reference objects.
#' @value A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' @references

saveRDSls<-function(list.object, dir.save, ascii = FALSE, version = NULL,
                    compress = TRUE, refhook = NULL){
  if(!dir.exists(dir.save)){ dir.create(dir.save) }
  invisible(lapply(seq_along(list.object),function(i){
    saveRDS(
      object = list.object[[i]],
      file = file.path(dir.save,paste0(names(list.object[i]),".RDS")),
      ascii = ascii, version = version, compress = compress, refhook = refhook)
  }))
}
