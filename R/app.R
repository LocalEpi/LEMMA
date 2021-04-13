#' @title Shiny interface
#' @description write me!
#' @export
run_app <- function() {
  appDir <- system.file("shiny", package = "LEMMA")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `LEMMA`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}