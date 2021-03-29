server <- function(input, output, session) {
  output$summary <- shiny::renderPrint({
    dataset <- get(input$dataset, "package:datasets")
    summary(dataset)
  })
  
  output$table <- shiny::renderTable({
    dataset <- get(input$dataset, "package:datasets")
    dataset
  })
}