#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinydashboard)
library(readr)
#library(DT)
library(dplyr)
library(tidyr)
library(reactable)
library(pracma)
library(jsonlite)

options(shiny.maxRequestSize = 300*1024^2)

source('hepQuan.R')
source('qualSearch.R')
source('quanSummary.R')
source('quanDetails.R')
source('shiftGen.R')
source('formatDatabase.R')

# Define UI for application that draws a histogram
ui <- dashboardPage(
  dashboardHeader(
    title = "HS Search Tool"
  ),
  dashboardSidebar(
    #Customize the sidebar in this section of code
    sidebarMenu(
      menuItem("Home", tabName='home',icon = icon("home")),
      tags$hr(style="border-color: grey;"),
      
      menuItem("HepQual",tabName='MS',icon=icon("search")),
      
      menuItem("HepQuant",tabName='HepQuan',icon=icon("search")),
      
      tags$hr(style="border-color: grey;")
      
      #menuItem('Contact',tabName='help')
    )
    #END OF SIDEBAR CODE
  ),
  dashboardBody(
    #Customize the body in this section of code
    tabItems(
      #home-----
      tabItem(tabName = "home", 
        fluidRow(
          column(
            width = 10,
            h4("Hello. This is home.")
          )
        )
      ),
      tabItem(tabName = "MS",
        #hep qual searching
        fluidRow(
          box(title=strong("HepQual"),status = "info",width = 3,
              selectizeInput(
                 "db", "Select DB:",
                 choices = c(
                   "4-30 mer db" = "db/hs_4_30mer_backbone.tsv"
                 ),
                 multiple = FALSE
              ),
              checkboxGroupInput("qual_reducing_end", "Reducing End:",
                c(
                  "Normal" = "normal",
                  "pNP" = "pnp",
                  "1,6-Anhydro" = "1-6-anhydro"
                )
              ),
              checkboxGroupInput("qual_non_end", "Non-Reducing End:",
                c(
                  "Normal" = "normal",
                  "ΔUA" = "dua"
                )
              ),
              numericInput( 
                "mz", 
                "m/z", 
                min = 0,
                max = 5000,
                value = 1195.3236
              ),
              numericInput( 
                "charge", 
                "Charges", 
                min = 1,
                max = 15,
                value = 3
              ),
              numericInput( 
                "ppm", 
                "PPM", 
                min = 0,
                max = 30,
                value = 5
              ),
              selectizeInput('iso_peak', 'Is input m/z a monoisotopic peak?', 
                choices = c('Yes' = 'yes', 'No' = 'no'),
                multiple = FALSE
              ),
              h4('Adducts'),
              checkboxInput("show_qualna", "Show Na", value = FALSE),
              conditionalPanel(
                condition = "input.show_qualna == true",
                numericInput("qual_Na", "Na", min = 0, max = 10, value = 0)
              ),
              checkboxInput("show_qualnh", "Show NH3", value = FALSE),
              conditionalPanel(
                condition = "input.show_qualnh == true",
                numericInput("qual_NH3", "NH3", min = 0, max = 10, value = 0)
              ),
              checkboxInput("show_qualmn", "Show Mn", value = FALSE),
              conditionalPanel(
                condition = "input.show_qualmn == true",
                numericInput("qual_Mn", "Mn", min = 0, max = 10, value = 0)
              ),
              checkboxInput("show_qualfa", "Show FA", value = FALSE),
              conditionalPanel(
                condition = "input.show_qualfa == true",
                numericInput("qual_FA", "FA", min = 0, max = 10, value = 0)
              ),
              actionButton("qual_search", "Search")
          ),
          box(title="Result:",width = 9,
              solidHeader = FALSE,status = "success",
              reactableOutput("db_table")
          )
        )
      ),
      #Quantative Filter & Input
      tabItem(tabName ="HepQuan",
        fluidRow(
          box(title=strong("HepQuan"),status = "info", width = 3,
              fileInput("importConfig", "Import config", accept = ".json"),
              fileInput("scan", "Upload scan.csv",
                multiple = FALSE,
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
              ),
              fileInput("iso", "Upload iso.csv",
                multiple = FALSE,
                accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")
              ),
              selectizeInput("db", "Select DB:",
                choices = c(
                            "4-30 mer db" = "db/hs_4_30mer_backbone.tsv"
                          ),
                multiple = FALSE
              ),
              checkboxGroupInput("quan_reducing_end", "Reducing End:",
               c(
                 "Normal" = "normal",
                 "pNP" = "pnp",
                 "1,6-Anhydro" = "1-6-anhydro"
               )
              ),
              checkboxGroupInput("quan_non_end", "Non-Reducing End:",
               c(
                 "Normal" = "normal",
                 "ΔUA" = "dua"
               )
              ),
              h4('Adducts'),
              checkboxInput("show_nh", "Show NH3", value = FALSE),
              conditionalPanel(
                condition = "input.show_nh == true",
                numericInput("quan_NH3", "NH3", min = 0, max = 10, value = 0)
              ),
              checkboxInput("show_na", "Show Na", value = FALSE),
              conditionalPanel(
                condition = "input.show_na == true",
                numericInput("quan_Na", "Na", min = 0, max = 10, value = 0)
              ),
              checkboxInput("show_mn", "Show Mn", value = FALSE),
              conditionalPanel(
                condition = "input.show_mn == true",
                numericInput("quan_Mn", "Mn", min = 0, max = 10, value = 0)
              ),
              checkboxInput("show_fa", "Show FA", value = FALSE),
              conditionalPanel(
                condition = "input.show_fa == true",
                numericInput("quan_FA", "FA", min = 0, max = 10, value = 0)
              ),
              numericInput("quan_ppm", "PPM", min = 0, max = 100, value = 15),
              h4('DP Range'),
              splitLayout(
                numericInput("dp_lwr", "Start", min = 0, max = 1000, value = 4),
                numericInput("dp_upr", "End", min = 0, max = 1000, value = 30)
              ),
              h4('Elution Time (Min)'),
              splitLayout(
                numericInput("start", "Start", min = 0, max = 1000, value = 25),
                numericInput("end", "End", min = 0, max = 1000, value = 50)
              ),
              splitLayout(
                numericInput("minscan", "Min Scan Number", min = 0, max = 200, value = 15)
              ),
              checkboxInput("show_adv", "Advanced Settings", value = FALSE),
              conditionalPanel(
                condition = "input.show_adv == true",
                numericInput("window", "Gaussian Smoothing Window", min = 0, max = 10000, value = 100),
                numericInput("nups", "findpeaks nups", min = 0, max = 100, value = 3),
                numericInput("ndowns", "findpeaks ndowns", min = 0, max = 100, value = 3),
                numericInput("threshold", "findpeaks threshold", min = 0, max = 100, value = 2),
                numericInput("ticfilter", "TIC grouping filter", min = 0, max = 1, step = 0.01, value = 0.95)
              ),
              downloadButton("exportConfig", "Export config"),
              actionButton("quan_search", "Search")
          ),
          box(title=strong("Result:"),width = 9,
            downloadButton("downloadRawData", "Download Raw Data"),
            downloadButton("downloadMergeData", "Download Merge Data"),
            reactableOutput("quanResults"),
            reactableOutput("quanDetails")
          )
        )
      )
    )
  )
)

# Define server logic here
server <- function(input, output, session) {

    #Quantitative search reactive - only updates on quan_search event activation
    quan_search <- eventReactive(input$quan_search, {
      withProgress(                
        message = 'Calculating',
        detail = 'This may take a while...', 
        value = 0, {
          req(input$scan)
          scan <-read.csv(input$scan$datapath, header = TRUE)
          req(input$iso)
          iso <- read.csv(input$iso$datapath, header = TRUE)
          quan_outp <- hepQuan(
            scan, iso, input$quan_ppm, input$db, 
            input$minscan, input$start, input$end, 
            input$dp_lwr, input$dp_upr,
            input$quan_Na, input$quan_NH3, input$quan_Mn, input$quan_FA,
            input$window, input$nups, input$ndowns, input$threshold, input$ticfilter,
            input$quan_reducing_end, input$quan_non_end
          )
        }
      )
    })
    
    #### QUAL SEARCH EVENT FUNCTION ####
    qual_search <- eventReactive(input$qual_search, {
      withProgress(                
        message = 'Calculating',
        detail = 'This may take a while...', 
        value = 0, {
          qual_outp <- qualSearch(
            input$mz, input$charge, 
            input$ppm, input$iso_peak, input$db, 
            input$qual_Na, input$qual_NH3, input$qual_Mn, input$qual_FA,
            input$qual_reducing_end, input$qual_non_end
          )
        }
      )
    })
    
    #Quan Search Event Listener
    observeEvent(input$quan_search, {
      #get the quan search summary
      summary <- quanSummary(quan_search())
      
      output$quanResults <- renderReactable(({
        reactable(summary,
          columns = list(
            name = colDef(
              minWidth = 120, 
              show = T
            ),
            abundance = colDef(
              show = T
            ),
            score = colDef(
              show = T
            )
          ),
          defaultColDef = colDef(
            show = F, 
            minWidth = 60
          ),
          onClick = JS("function(rowInfo, colInfo) {
            Shiny.setInputValue('row_click', { 
              index: rowInfo.index, 
              name: rowInfo.values.name
            }, { 
              priority: 'event' 
            })
          }")
        )
      }))
    })
    
    
    #### QUAL REDUCING END CHECKBOX CONTROL ####
    prev_qual_end <- reactiveVal(NULL)
    
    observeEvent(input$qual_reducing_end, {
      current <- input$qual_reducing_end
      previous <- prev_qual_end()
      
      if ("pnp" %in% current && !("pnp" %in% previous)) {
        if("normal" %in% current) current <- c("normal", "pnp") else current <- "pnp"
      } else if("1-6-anhydro" %in% current && !("1-6-anhydro" %in% previous)) {
        if("normal" %in% current) current <- c("normal", "1-6-anhydro") else current <- "1-6-anhydro"
      }
      updateCheckboxGroupInput(session, "qual_reducing_end", selected = current)
      prev_qual_end(current)
    }, ignoreNULL = FALSE)
    #### QUAL REDUCING END CHECKBOX CONTROL ####
    
    
    #### QUAL NON-REDUCING END CHECKBOX CONTROL ####
    prev_qual_non <- reactiveVal(NULL)
    
    observeEvent(input$qual_non_end, {
      current <- input$qual_non_end
      previous <- prev_qual_non()
      
      if ("dua" %in% current && !("dua" %in% previous)) {
        if("normal" %in% current) current <- c("normal", "dua") else current <- "dua"
      }
      updateCheckboxGroupInput(session, "qual_non_end", selected = current)
      prev_qual_non(current)
    }, ignoreNULL = FALSE)
    #### QUAL NON-REDUCING END CHECKBOX CONTROL ####
    
    
    #### QUANT REDUCING END CHECKBOX CONTROL ####
    prev_quan_end <- reactiveVal(NULL)
    
    observeEvent(input$quan_reducing_end, {
      current <- input$quan_reducing_end
      previous <- prev_quan_end()

      if ("pnp" %in% current && !("pnp" %in% previous)) {
        if("normal" %in% current) current <- c("normal", "pnp") else current <- "pnp"
      } else if("1-6-anhydro" %in% current && !("1-6-anhydro" %in% previous)) {
        if("normal" %in% current) current <- c("normal", "1-6-anhydro") else current <- "1-6-anhydro"
      }
      updateCheckboxGroupInput(session, "quan_reducing_end", selected = current)
      prev_quan_end(current)
    }, ignoreNULL = FALSE)
    #### QUANT REDUCING END CHECKBOX CONTROL ####
    
    
    #### QUANT NON-REDUCING END CHECKBOX CONTROL ####
    prev_quan_non <- reactiveVal(NULL)
    
    observeEvent(input$quan_non_end, {
      current <- input$quan_non_end
      previous <- prev_quan_non()
      
      if ("dua" %in% current && !("dua" %in% previous)) {
        if("normal" %in% current) current <- c("normal", "dua") else current <- "dua"
      }
      
      updateCheckboxGroupInput(session, "quan_non_end", selected = current)
      prev_quan_non(current)
    }, ignoreNULL = FALSE)
    #### QUANT NON-REDUCING END CHECKBOX CONTROL ####
    
    
    observeEvent(input$row_click, {
      details <- getDetails(input$row_click$name, quan_search())
      output$quanDetails <- renderReactable({
        reactable(details, 
          columns = list(
            abundance = colDef(
              format = colFormat(digits = 0)
            ),
            time = colDef(
              format = colFormat(digits = 2)
            )
          )
        )
      })
    })
    
    #Qualitative Search Event Listener
    observeEvent(input$qual_search, {
      output$db_table <- renderReactable({
        qual_search()
      })
    })
    
    output$downloadRawData <- downloadHandler(
      filename = function() {
        paste0("raw_result_", Sys.Date(), ".csv")
      },
      content = function(fname) {
        outp_db <- quan_search()
        write.table(outp_db, fname, sep=",", row.names=FALSE)
      },
      contentType = 'text/csv'
    )
    
    output$downloadMergeData <- downloadHandler(
      filename = function() {
        paste0("merge_result_", Sys.Date(), ".csv")
      },
      content = function(fname) {
        outp_db <- quanSummary(quan_search())
        write.table(outp_db, fname, sep=",", row.names=FALSE)
      },
      contentType = 'text/csv'
    )
    
    cfg_num  <- c("quan_ppm", "dp_lwr", "dp_upr", "start", "end", "minscan",
                  "quan_Na", "quan_NH3", "quan_Mn", "quan_FA",
                  "window", "nups", "ndowns", "threshold", "ticfilter")
    cfg_text <- c("quan_reducing_end", "quan_non_end")
    
    output$exportConfig <- downloadHandler(
      filename = function() paste0("hepquan_config_", Sys.Date(), ".json"),
      content = function(file) {
        ids <- c(cfg_num, cfg_text)
        cfg <- setNames(lapply(ids, function(id) input[[id]]), ids)
        cfg$version <- 1
        write_json(cfg, file, auto_unbox = TRUE, pretty = TRUE, null = "null")
      }
    )
    
    observeEvent(input$importConfig, {
      cfg <- tryCatch(read_json(input$importConfig$datapath, simplifyVector = TRUE),
                      error = function(e) NULL)
      if (is.null(cfg)) {
        showNotification("Could not read config file", type = "error")
        return()
      }
      
      for (id in intersect(names(cfg), cfg_num))
        if (is.numeric(cfg[[id]])) updateNumericInput(session, id, value = cfg[[id]])
      for (id in intersect(names(cfg), cfg_text))
        updateCheckboxGroupInput(session, id, selected = cfg[[id]])
      
      # Reveal the conditional adduct panels that now hold non-zero values
      updateCheckboxInput(session, "show_na", value = isTRUE(cfg$quan_Na > 0))
      updateCheckboxInput(session, "show_nh", value = isTRUE(cfg$quan_NH3 > 0))
      updateCheckboxInput(session, "show_mn", value = isTRUE(cfg$quan_Mn > 0))
      updateCheckboxInput(session, "show_fa", value = isTRUE(cfg$quan_FA > 0))
      updateCheckboxInput(session, "show_adv", value = isTRUE(cfg$window != 100 || 
                                                              cfg$nups != 3 ||
                                                              cfg$ndowns != 3 ||
                                                              cfg$threshold != 2 ||
                                                              cfg$ticfilter != .95))
      
      showNotification("Filters loaded", type = "message")
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
