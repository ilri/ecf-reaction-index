#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# The input data is required to have the following columns:
# experimentName
# animalID
# date, dayPostChallenge
# temp
# schizontsLocal, schizontsContra
# piroplasms
# WBC
# exitDay, exitType
#
# optional column: group

library(shiny)
library(DT)
library(readr)
library(readxl)
library(dplyr)
library(forcats) # for as_factor
library(plotly) # for interactive graphics
library(scales) # for nicer scales
library(magrittr) # for %<>% only


# According to https://stackoverflow.com/questions/74562346/prevent-ggplotly-to-change-the-legends-style-created-in-ggplot,
# we have to fix plotly's handling of legends for geom_line(),
# because it ignores its show.legend option
solid_lines_legend <- function(plotly_obj) {
  # the input to this function is a plotly output.
  # this fix borrows heavily from the one by https://stackoverflow.com/users/5329073/kat
  # here: https://stackoverflow.com/questions/74562346/prevent-ggplotly-to-change-the-legends-style-created-in-ggplot
  # BEWARE: lines that are dash-only WILL NOT appear in the legend
  lapply(1:length(plotly_obj$x$data),
         function(j) {
           if(plotly_obj$x$data[[j]]$mode == "lines") {
             if(plotly_obj$x$data[[j]]$line$dash == "dash" |
                nchar(plotly_obj$x$data[[j]]$name) == 0) # anonymous line: do not legend
               plotly_obj$x$data[[j]]$showlegend <<- F
             else
               plotly_obj$x$data[[j]]$showlegend <<- T
           } # endif
         }) #endfunction j #end lapply
  plotly_obj
} #endfunction solid_lines_legend



# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("ECF Reaction Index calculation"),
  
  # Sidebar with inputs
  sidebarLayout(
    sidebarPanel(
      fileInput(
        "dataSourceFile",
        "Data Source (csv, xls or xlsx)",
        accept = c(".csv", ".xls", ".xlsx")
      ),
      selectInput(
        "experiment",
        "Which experiment? (challenge date and trial duration displayed)",
        choices = character(),
        multiple = F
      ),
      radioButtons(
        "displayType",
        "Table display",
        choices = c(
          `Condensed (clinical observations only)` = "Condensed",
          `Full (all columns, including calculated variables)` = "Full"),
        selected = "Condensed"
      ),
      radioButtons(
        "colourMapping",
        "Colour mapping in graphs",
        choices = c(`Individual animals` = "animalID", `Experimental groups` = "group"),
        selected = "animalID"
      ),
      #textOutput("warningOnGroups"),
      #br(),
      
      selectInput(
        "whichAnimal",
        "Which animal(s)?",
        choices = c("(all)"),
        multiple = T
      ),
      numericInput(
        "RI_threshold",
        "Threshold to display on the ECF reaction index graphs (e.g. for humane endpoint)",
        value = 6,
        min = 0,
        max = 10,
        step = 0.1
      ),
      selectInput(
        "additionalVar",
        "What to visualize apart from reaction indices?",
        choices = character(),
        multiple = F
      ),
      width = 3
    ),
    
    
    # Show a plot of the generated distribution
    mainPanel(dataTableOutput("mainDataTable"),
              plotlyOutput("mainPlot"),
              plotlyOutput("additionalPlot"),
              width = 9)
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  RI_col <- "newRI" # the main reaction index to display
  
  output$warningOnGroups <- renderText({
    "(the above will trigger an error when \"Groups\" is selected but the input\
    data doesn't have the information on experimental groups: \ 
    please check your input file)"
  })
  
  # the following filters trim the table display when
  # input$displayType == "Condensed"
  data_cols_to_display <- c("experimentName", "date", "dayPostChallenge", "animalID", "temp",
                            "schizontsLocal", "schizontsContra", "piroplasms", "RBC", "Hgb", "PCV%", "newRI", "oldCALPR1")
  
  # we first instantiate the all_dataset variable, as a reactive (tibble) object:
  all_datasets <- reactive({
    req(input$dataSourceFile) # says what we require: silent output if not present
    file <- input$dataSourceFile
    extension <- tools::file_ext(file$datapath)
    validate(
      need(extension %in% c("csv", "xls", "xlsx"), "Please provide input data as a csv, xls or xslx file.")
    )
    
    if(extension == "csv")
      df <- read_csv(file$datapath, show_col_types = F)
    else
      df <- read_excel(file$datapath)
      
    # BEWARE: hardcoded "experimentName" and "date" column names
    df$experimentName <- factor(df$experimentName)
    # and useless to keep everything of the date, including the time:
    df$date <- lubridate::date(df$date)
    return(df)
  })
  
  # we dynamically recompute the list of experiments to pick from:
  observe({
    all_datasets() %>% group_by(experimentName) %>% summarize(
      cha_date = as.character(lubridate::date(min(date, na.rm = T) - 1)),
      num_days = first(experimentDuration, na_rm = T)) %>%
      mutate(
        fullname = paste0(experimentName, " (", cha_date, " + ", num_days, " days)"),
        .keep = "all") -> temp_table
    temp_table %>% pull(experimentName) -> named_vec
    temp_table %>% pull(fullname) -> names(named_vec)
    # the display values will contain the challenge date (the lowest date minus one day)
    updateSelectInput(session, "experiment", choices = named_vec)
  })
  
  # reactively set the dataset
  dataset <- reactive({
    req(input$experiment)
    return(all_datasets() %>% filter(experimentName == input$experiment))
  })
  
  # and set the list of columns to pick from, to filter on a column
  observe({
    req(input$experiment)
    updateSelectInput(session, "additionalVar", choices = setdiff(colnames(dataset()), c(RI_col, "animalID")), selected = "temp")
  })
  
  # we dynamically recompute the list of animal IDs to pick from:
  observe({
    req(input$experiment)
    updateSelectInput(session, "whichAnimal",
                      choices = c("(all)", unique(as.character(
      `[[`(dataset(), "animalID")
    ))))
  })
  
  
  # main data table output
  output$mainDataTable <- renderDataTable({
    req(input$experiment)
    
    d <- dataset()
    if(input$displayType == "Condensed")
      d %<>% select(any_of(data_cols_to_display))

    if(is.null(input$whichAnimal) | "(all)" %in% input$whichAnimal)
        d
    else
        d %>%
        filter(as.character(.data[["animalID"]]) %in% input$whichAnimal)
  
  }, options = list(pageLength = 5))
  
  
  # careful in the following: we use hardcoded column names, including animalID
  # TRICK: there is a computed column experimentDuration, NA only for the virtual
  # observations. We use that to filter the stuff we are plotting.
  # We plot dashed lines underneath, and then solid lines on top.
  # In all the graphs, setting group = animalID makes sure we keep one line per animal
  # no matter what the other aesthetics are.
  output$mainPlot <- renderPlotly({
    req(input$experiment)
    validate(
      need(input$colourMapping == "animalID" | ("group" %in% colnames(dataset()) & any(!is.na(dataset() %>% pull(group)))),
           "You asked for a colour mapping on groups, but your input table doesn't contain \
           any information on groups for this trial.")
    )
    
    if (is.null(input$whichAnimal) | "(all)" %in% input$whichAnimal)
      { dataset() %>% mutate(animalID = as_factor(animalID)) %>%
        ggplot(aes(group = animalID, x = dayPostChallenge, y = newRI, color = !!sym(input$colourMapping))) +
        geom_line(linetype = 2, show.legend = FALSE) + # all lines dashed first
        geom_line(data = ~filter(.x, !is.na(experimentDuration))) + # solid lines
        labs(title = paste0("Reaction indices (", RI_col, ") for all animals"),
                            x = "Day post challenge", y = "ECF reaction index") +
        scale_x_continuous(breaks = breaks_width(1)) -> p
      if(isTruthy(input$RI_threshold)) p + geom_hline(yintercept = input$RI_threshold) -> p
      ggplotly(p) %>% solid_lines_legend()
    } else {
      dataset() %>% mutate(animalID = as_factor(animalID)) %>%
        filter(animalID %in% input$whichAnimal) %>%
        ggplot(aes(group = animalID, x = dayPostChallenge, y = newRI, color = !!sym(input$colourMapping))) +
        geom_line(linetype = 2, show.legend = FALSE) + # all lines dashed first
        geom_line(data = ~filter(.x, !is.na(experimentDuration))) + # solid lines
        labs(
        title = paste0("Reaction index (", RI_col, ") for animal(s) ", paste(input$whichAnimal, collapse = ", ")),
        x = "Day post challenge",
        y = "ECF reaction index") +
        scale_x_continuous(breaks = breaks_width(1)) -> p
      if(isTruthy(input$RI_threshold)) p + geom_hline(yintercept = input$RI_threshold) -> p
      ggplotly(p) %>% solid_lines_legend()
    }
  })
  
  # careful in the following: we use hardcoded column names, including animalID
  output$additionalPlot <- renderPlotly({
    req(input$experiment, input$additionalVar)
    validate(
      need(input$colourMapping == "animalID" | ("group" %in% colnames(dataset()) & any(!is.na(dataset() %>% pull(group)))),
           "You asked for a colour mapping on groups, but your input table doesn't contain \
           any information on groups for this trial.")
    )
    
    if (is.null(input$whichAnimal) | "(all)" %in% input$whichAnimal)
    { dataset() %>% mutate(animalID = as_factor(animalID)) %>%
        ggplot(aes(group = animalID, x = dayPostChallenge, y = !!sym(input$additionalVar), color = !!sym(input$colourMapping))) +
        geom_line(linetype = 2, show.legend = FALSE) + # all lines dashed first
        geom_line(data = ~filter(.x, !is.na(experimentDuration))) + # solid lines
        labs(title = paste(input$additionalVar, "for all animals"), x = "Day post challenge", y = input$additionalVar) +
        scale_x_continuous(breaks = breaks_width(1)) -> p
      ggplotly(p) %>% solid_lines_legend()
    } else {
      dataset() %>% mutate(animalID = as_factor(animalID)) %>%
        filter(animalID %in% input$whichAnimal) %>%
        ggplot(aes(group = animalID, x = dayPostChallenge, y = !!sym(input$additionalVar), color = !!sym(input$colourMapping))) +
        geom_line(linetype = 2, show.legend = FALSE) + # all lines dashed first
        geom_line(data = ~filter(.x, !is.na(experimentDuration))) + # solid lines
        labs(
        title = paste(input$additionalVar, "for animal(s)", paste(input$whichAnimal, collapse = ", ")),
        x = "Day post challenge",
        y = input$additionalVar) +
        scale_x_continuous(breaks = breaks_width(1)) -> p
      ggplotly(p) %>% solid_lines_legend()
    }
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
