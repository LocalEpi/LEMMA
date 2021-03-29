ui <- shiny::fluidPage(
  shiny::selectInput("dataset", label = "Dataset", choices = ls("package:datasets")),
  shiny::verbatimTextOutput("summary"),
  shiny::tableOutput("table")
)