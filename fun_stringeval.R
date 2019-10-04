##FUNCTIONS

#' @description Check if an alphanumerical string contains any upper case
#'              character.
#' 
#' @param str A \code{character} vector to be evaluated looking for any
#'                               uppercase character.
#' @value A \code{logical}.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' @references

is.upper<-function(str){ grepl(pattern = "[A-Z]+", x = str) }