ui <- shiny::fluidPage(
  shiny::fluidRow(
    shiny::column(3, 
           shiny::numericInput("lambda1", label = "lambda1", value = 3),
           shiny::numericInput("lambda2", label = "lambda2", value = 5),
           shiny::numericInput("n", label = "n", value = 1e4, min = 0),
           shiny::actionButton("simulate", "Simulate!")
    ),
    shiny::column(9, plotOutput("hist"))
  )
)

# ui <- shiny::fluidPage(
#   # shiny::fileInput("upload", NULL) # this is how we'll input files
#   shiny::fluidRow(
#     shiny::selectInput("dataset", label = "Dataset", choices = ls("package:datasets")),
#     shiny::verbatimTextOutput("summary"),
#     shiny::tableOutput("table")
#   )
# )