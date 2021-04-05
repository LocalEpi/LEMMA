expected_sheets <- c(
  "Parameters with Distributions","Interventions","Model Inputs","Data",
  "Vaccine Distribution","Vaccine Doses - Observed","Vaccine Doses - Future","Variants",
  "PUI Details","Internal"
)

server <- function(input, output, session) {
  
  # reactive: excel upload
  xlsx_input <- reactive({
    req(input$upload)
    ext <- tools::file_ext(input$upload$name)
    if(ext != "xlsx"){
      validate("Invalid file; Please upload a .xlsx file")
    }
    sheets <- readxl::excel_sheets(path = normalizePath(input$upload$datapath))
    if(!identical(expected_sheets,sheets)){
      validate(cat("Invalid file; File needs 10 named sheets: ",paste0(expected_sheets,collapse = ", ")))
    }
    id <- showNotification("Reading data", duration = NULL, closeButton = FALSE,type = "message")
    on.exit(removeNotification(id), add = TRUE)
    LEMMA:::ReadInputs(path = input$upload$datapath)
  })
  
  # reactive: inputs
  LEMMA_inputs <- reactive({
    req(xlsx_input())
    id <- showNotification("Generating LEMMA parameters", duration = NULL, closeButton = FALSE,type = "message")
    on.exit(removeNotification(id), add = TRUE)
    LEMMA:::ProcessSheets(xlsx_input())
  })
  
  # reactive: LEMMA run from excel upload
  LEMMA_excel_run <- eventReactive(input$LEMMA_xlsx, {
    req(LEMMA_inputs())
    id <- showNotification("Running LEMMA", duration = NULL, closeButton = FALSE,type = "message")
    on.exit(removeNotification(id), add = TRUE)
    LEMMA:::CredibilityIntervalData(inputs = LEMMA_inputs(),fit.to.data = NULL)
  })
  
  # reactive: LEMMA excel output
  LEMMA_excel_out <- reactive({
    req(LEMMA_excel_run())
    id <- showNotification("Creating .xlsx output file", duration = NULL, closeButton = FALSE,type = "message")
    on.exit(removeNotification(id), add = TRUE)
    LEMMA:::GetExcelOutputData(LEMMA_excel_run()$projection, LEMMA_excel_run()$fit.to.data, LEMMA_excel_run()$inputs)
  })
  
  # output: checker to let users know excel uploaded
  output$xlsx_check_txt <- renderText({
    req(xlsx_input())
    paste0("File successfully uploaded, size ",signif(input$upload$size/1e3,digits = 6),"Kb")
  })
  
  # output: download the sample template
  output$download_template <- downloadHandler(
    filename = function() {
      return("example.xlsx")
    },
    content = function(file) {
      file.copy(from = normalizePath(path = paste0(path.package("LEMMA"),"/extdata/template.xlsx")),to = file)
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )
  
  # output: downloald excel output
  output$download_pdf_out <- downloadHandler(
    filename = function() {
      # req(LEMMA_inputs())
      return("output.xlsx")
    },
    content = function(file) {
      req(LEMMA_excel_out())
      openxlsx::write.xlsx(LEMMA_excel_out(), file = file)
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )
  
  
  
  # DEBUGGING
  output$table <- renderTable({
    LEMMA_excel_run()$projection
  })
}