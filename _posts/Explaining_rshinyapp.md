---
layout: post
title: "Building an rshiny SMR map"
author: "Rasim M Musal"
date: "2023-10-07"
output:
  html_document:
   theme: darkly
   highlight: espresso
   toc: true
   keep_md: yes
   toc_float: true
   toc_collapsed: true
   toc_depth: 1
   number_sections: true
tags: [rshiny, maps,SMR]
always_allow_html: true

---



# Discussing Rshiny

R allows you to use the rshiny package in order to build web applications. An excellent introduction to building rshiny apps are given by [Posit](https://shiny.posit.co/r/articles/start/build/) and [Bookdown](https://bookdown.org/hadrien/how_to_build_a_shiny_app_from_scratch/).

Here we assume the reader has working knowledge of rshiny and present the code that my graduate student Anjali Goel have presented in 2022. The application is a map of California  that shows how the standardized mortality ratios (SMR) of Covid-19 mortality change over 77 biweeks.  
The application is at [shiny server](https://mmusal.shinyapps.io/timeseriesofSMR/).  

# Application code 


```r
knitr::opts_chunk$set(eval = FALSE)
library(shiny)
library(rgdal)
library(DT)
library(dygraphs)
library(xts)
library(leaflet)
library(shinyjs)

#It is important to comment out the lines starting with wd and setwd after running otherwise the app will not deploy.
#wd='C:/Users/rm84/Documents/timeseriesofSMR'
#setwd(wd)

data <- as.data.frame(read.csv("./data/SMR.csv"))
names(data)[names(data) == "County"] <- "county"
##Reading shape file
map <- readOGR("./map/tl_2021_us_county.shp")

#Filtering to include only counties in CA.
map=map[(map$STATEFP %in% '06'),]

# ui object
ui <- fluidPage(
  useShinyjs(),
  titlePanel(p("Spatial app", style = "color:#3474A7")),
  actionButton(inputId = "button", label = "Show/Hide Panel"),
  sidebarLayout(
    div(id = "sidebar", sidebarPanel(
      selectInput(
        inputId = "variableselected",
        label = "Select variable",
        choices = c("SMR")
      ),

      sliderInput("periods", 
                  "Periods (1 - 77):",
                  min = 1, max = 77,
                  value = 1, step = 1,
                  animate =
                    animationOptions(interval = 1000, loop = FALSE)
      ),
      
      
    )),
    
    mainPanel(
      leafletOutput(outputId = "map", width = "100%", height = 800)
    )
  )
)

# server()
server <- function(input, output) {
  shinyjs::hide(id = "variableselected")
  ## observe the button being pressed
  observeEvent(input$button, {
   
    if(input$button %% 2 == 1){
      shinyjs::hide(id = "sidebar")
    }else{
      shinyjs::show(id = "sidebar")
    }
  })

  output$map <- renderLeaflet({
    
    # Add data to map
    datafiltered <- data[which(data$Time == input$periods), ]
    ordercounties <- match(map@data$NAME, datafiltered$county)
    map@data <- datafiltered[ordercounties, ]
    
    # Create variableplot
    # ADD this to create variableplot
    map$variableplot <- as.numeric(
      map@data[, input$variableselected])
    
    # Create leaflet
    # CHANGE map$cases by map$variableplot
    pal <- colorBin("Reds", domain = map$variableplot , bins = c(0,1,2,3,4,5,6,7,51)
                        )
    
    # CHANGE map$cases by map$variableplot
    labels <- sprintf("%s: %g", map$county, map$variableplot) %>%
      lapply(htmltools::HTML)
    
    # CHANGE cases by variableplot
    l <- leaflet(map) %>%
      addTiles() %>%
      addPolygons(
        fillColor = ~ pal(variableplot),
        color = "white",
        dashArray = "3",
        fillOpacity = 0.4,
        label = labels
      ) %>%
      # CHANGE cases by variableplot
      leaflet::addLegend(
        pal = pal, values = ~variableplot,
        opacity = 0.7, title = NULL
      )
    
  })
}

# shinyApp()
shinyApp(ui = ui, server = server)
```

# Discussing the calculation for SMR

The application itself resides in the [shiny server](https://mmusal.shinyapps.io/timeseriesofSMR/)
and shows 77 biweeks (t=1,...T=77) of Covid-19 Standardized Mortality Ratio (SMR) in the 58 counties of California. The SMR in this application is calculated for a given time period "t".
Below $Pop_{i}$ is population of county "i". We used the 2020 and 2021 population (there will be more explanation when I document the data work).
$\omega$ represents the population weight of county "i" among N (58) counties in California.
\[\omega_{i}=\frac{Pop_{i}}{\sum_{i}^{N}Pop_{i}}\]
In other applications expected mortality at time period $t$, county $i$, is adjusted by age/sex but we prefer to use these solely as model variables and therefore $E_{t,i}$ is, 
\[E_{t,i}=(\sum_{i}y_{t,i})*\omega_{i}\]
Therefore you can think about $E_{t,i}$ as the expected number of people who would die of Covid-19 if mortality risk was uniformly distributed across the population (death from Covid-19 was equally likely).

In conclusion SMR is calculated as

\[SMR_{t,i}=\frac{y_{t,i}}{E_{t,i}}\]

# Noting issues on SMR values

There are a few things to note about SMR among these 58 counties across 77 biweeks. 

0. One of the first things people note is that we are not using SIR models where population changes at every time period with deaths and births. The main reason we are not using SIR models is that we do not have such data directly available. On the other hand since we are working with ratios the main assumption we make is that the county population ratios $\omega$ do not change.     

1. Each SMR value is calculated based on the total observed deaths in a biweek. Especially in the first set of biweeks most of the counties have no recorded deaths. This means a single death occurring in a  county with an $E_{t,i}$ of 0.05 lead to an $SMR_{t,i}$ of 20. Drawing inferences based on a few biweeks of data would be an example of [the law of small numbers](https://en.wikipedia.org/wiki/Law_of_small_numbers). We can not draw healthy inferences on what affects SMR (our main goal) based on few weeks of data. 

2. Connected to 1, we have to remember smaller counties will have higher variances. 

3. 1 and 2 will require us to come up with some sort of smoothing mechanism if we are to draw inferences of how various variables effect SMR.   

4. When we create visualizations we need to be aware as to what we show. Every biweek SMR is calculated independent of the previous biweeks. A county that does not perform well will create higher E values therefore reducing the SMR for the other counties. We could develop a cumulative SMR value that can be in another post. 








