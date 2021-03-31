# server <- function(input, output, session) {
#   
#   # Create a reactive expression
#   dataset <- shiny::reactive({
#     get(input$dataset, "package:datasets")
#   })
# 
#   output$summary <- shiny::renderPrint({
#     summary(dataset())
#   })
# 
#   output$table <- shiny::renderTable({
#     dataset()
#   })
# }


  freqpoly <- function(x1, x2, binwidth = 0.1, xlim = c(-3, 3)) {
    df <- data.frame(
      x = c(x1, x2),
      g = c(rep("x1", length(x1)), rep("x2", length(x2)))
    )
    
    ggplot(df, aes(x, colour = g)) +
      geom_freqpoly(binwidth = binwidth, size = 1) +
      coord_cartesian(xlim = xlim)
  }
  
  t_test <- function(x1, x2) {
    test <- t.test(x1, x2)
    
    # use sprintf() to format t.test() results compactly
    sprintf(
      "p value: %0.3f\n[%0.2f, %0.2f]",
      test$p.value, test$conf.int[1], test$conf.int[2]
    )
  }
  
  server <- function(input, output, session) {
    x1 <- reactive({
      input$simulate
      rpois(input$n, input$lambda1)
    })
    x2 <- reactive({
      input$simulate
      rpois(input$n, input$lambda2)
    })
    output$hist <- renderPlot({
      freqpoly(x1(), x2(), binwidth = 1, xlim = c(0, 40))
    }, res = 96)
  }