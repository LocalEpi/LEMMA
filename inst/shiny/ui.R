ui <- shiny::navbarPage(
  "LEMMA (Local Epidemic Modeling for Management and Action)",   
  tabPanel("Model Structure", 
           includeMarkdown(path = normalizePath(path = paste0(path.package("LEMMA"),"/shiny/src/SEIRModel.md")))
  ),
  tabPanel("Excel Interface",
           fluidRow(
             column(4,
                    fileInput("upload", "Upload a spreadsheet"),
                    textOutput("xlsx_check_txt"),
                    br(),
                    HTML(r"(<label class="control-label" id="upload-label" for="upload">Download template spreadsheet</label>)"),
                    br(),
                    downloadButton("download_template", "Download")
              ),
             column(8,
                    includeMarkdown(path = normalizePath(path = paste0(path.package("LEMMA"),"/shiny/src/excel_input.md")))
              )
           )
  ),
  tabPanel("debugging", 
           fluidRow(
             tableOutput("files")
           )
  ),
  navbarMenu("subpanels", 
             tabPanel("panel 4a", "four-a"),
             tabPanel("panel 4b", "four-b"),
             tabPanel("panel 4c", "four-c")
  )
)