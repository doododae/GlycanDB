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
library(reactable)
library(pracma)

options(shiny.maxRequestSize = 300*1024^2)

source('hepQuan.R')
source('qualSearch.R')
source('quanSummary.R')
source('quanDetails.R')

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
          box(title=strong('User Guide'),status='warning',width = 10,
              solidHeader = FALSE,
              collapsible = TRUE,
              collapsed = FALSE,
              closable = FALSE
          )
        ),
        fluidRow(
          box(title=strong("HepQual"),status = "info",width = 3,
              h4('Qual Filter:'),
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
              selectizeInput(
                'iso_peak', 
                'Is input m/z a monoisotopic peak?', 
                choices = c('Yes' = 'yes', 'No' = 'no'),
                multiple = FALSE
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
              fileInput("scan", "Upload scan.csv",
                multiple = FALSE,
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
              ),
              fileInput("iso", "Upload iso.csv",
                multiple = FALSE,
                accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")
              ),
              h4('PPM'),
              numericInput("quan_ppm", "PPM",min = 0, max = 30, value = 15),
              h4('DP Range'),
              splitLayout(
                numericInput("dp_lwr", "Start",min = 0, max = 1000, value = 4),
                numericInput("dp_upr", "End",min = 0, max = 1000, value = 30)
              ),
              h4('Elution Time (Min)'),
              splitLayout(
                numericInput("start", "Start",min = 0, max = 1000, value = 25),
                numericInput("end", "End",min = 0, max = 1000, value = 50)
              ),
              splitLayout(
                numericInput("minscan", "Min Scan Number",min = 0, max = 200, value = 15)
              ),
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
    db_location = 'db/hs-library.tsv'
    
    db <- read.table(file = db_location, sep = '\t', header = TRUE)

    outp_db <- db
    
    #initial qualitative table render - i might be able to put it elsewhere?
    output$db_table <- renderReactable({ 
      reactable(outp_db, 
        # ALL COLUMNS (name, HexA, HexN, Ac, S, formula, neutral_mass, floating_Na, floating_NH3)
        columns = list(
          DP = colDef(show = F),
          formula = colDef(show = F),
          floating_Na = colDef(show = F),
          floating_NH3 = colDef(show = F)
        ),
        defaultColDef = colDef(show = T), 
        details = colDef(
          name = "More",
          details = JS("function(rowInfo) {
            return `Details for row: ${rowInfo.index}` +
              `<pre>${JSON.stringify(rowInfo.values, null, 2)}</pre>`
          }"),
          html = TRUE,
          width = 60
        )
      )
    })
    
    #Quantitative search reactive - only updates on quan_search event activation
    quan_search <- eventReactive(input$quan_search, {
      withProgress(                
        message = 'Calculating',
        detail = 'This may take a while...', 
        value = 0, {
          for (i in 1:10) {
            req(input$scan)
            scan <-read.csv(input$scan$datapath, header = TRUE)
            req(input$iso)
            iso <- read.csv(input$iso$datapath, header = TRUE)
            quan_outp <- hepQuan(
              scan, iso, input$quan_ppm, db, 
              input$minscan, input$start, input$end, 
              input$dp_lwr, input$dp_upr
            )
            return(quan_outp)
            incProgress(1/15)
            Sys.sleep(0.25)
          }
        }
      )
    })
    
    #Qualitative search reactive - only updates on qual_search event activation
    qual_search <- eventReactive(input$qual_search, {
      return(qualSearch(input$mz, input$charge, input$ppm, input$iso_peak, db))
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
        outp_db <- getDetails(quan_search())
        write.table(outp_db, fname, sep=",", row.names=FALSE)
      },
      contentType = 'text/csv'
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
