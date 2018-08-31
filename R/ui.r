

library(shiny)


# Define UI
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Model parameters"),

  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    selectInput("variable", "Model type:",
                c("Numerical Model"        = "Numer", 
                  "Individual-Based Model" = "IBM")),

    sliderInput("time",
                "Number of generations",
                min = 10,
                max = 200,
                step  = 10,
                value = 100),

    sliderInput("N_initial",
                "Initial prey density",
                min   = 0,
                max   = 1,
                step  = 0.01,
                value = 0.4),

    sliderInput("P_initial",
                "Initial predator density",
                min   = 0,
                max   = 1,
                step  = 0.01,
                value = 0.4),

    sliderInput("lambda_N",
                "Prey population growth rate",
                min   = 0,
                max   = 2,
                step  = 0.01,
                value = 0.53),

    sliderInput("lambda_P",
                "Predator population decline rate",
                min   = 0,
                max   = 2,
                step  = 0.01,
                value = 0.5),

    sliderInput("attack",
                "Attack rate of predators",
                min   = 0,
                max   = 2,
                step  = 0.01,
                value = 1.3),

    sliderInput("K",
                "Prey density carrying capacity",
                min   = 0,
                max   = 1,
                step  = 0.01,
                value = 0.8),

    sliderInput("b",
                "Predator births per attack",
                min = 0,
                max = 2,
                step  = 0.01,
                value = 1.6),

    sliderInput("dispY",
                "Prey dispersal (IBM only)",
                min = 0.01,
                max = 1,
                step  = 0.01,
                value = 1),

    sliderInput("dispP",
                "Predator dispersal (IBM only)",
                min = 0.01,
                max = 1,
                step  = 0.01,
                value = 1),

    actionButton("go", "Run"),

    actionButton("finish", "Finish!"),

    actionButton("reset", "Reset")

  ),

  # Show the caption and plot of the requested variable against mpg
   mainPanel(plotOutput(outputId="distPlot"))

))


