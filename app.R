#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(deSolve)
library(zeallot)

# Defining parameters
c(rA, rB, k1, k2, k3, k4) %<-% list(4,5,1,2,3,4)

ui <- fluidPage(

  # Application title
  titlePanel("A Damped Biological Oscillator"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("alpha1", "Degradation X (\u03B11)", min = 1, max = 20, value = 5),
      sliderInput("alpha2", "Degradation Y (\u03B12)", min = 1, max = 20, value = 5),
      sliderInput("beta1", "Formation X (\u03B21)", min = 1, max = 20, value = 5),
      sliderInput("beta2", "Formation Y (\u03B22)", min = 1, max = 20, value = 5),
      sliderInput("len", "Simulation Time", min = 1, max = 50, value = 2)
    ),
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
)



brusselator <- function (t,y,p) {
  with(as.list(c(y,p)), {
    dX <- rA*k1 - k2*rB*X + k3*(X^2)*Y - k4*X
    dY <- k2*rB*X - k3*(X^2)*Y
    list(c(X=dX, Y=dY))
  })
}

brussC <- cOde::funC(
  c(
    X = "k1*A - k2*B*X + k3*(X^2)*Y - k4*X",
    Y = "k2*B*X - k3*(X^2)*Y"
  ), modelname = "brusselC", compile=TRUE
)

k3 <- rnorm(0.1)

osc2 <- cOde::funC(
  c(
    X = "-b1*Y - a1*X",
    Y = "b2*X - a2*Y"
  ), modelname = "osc_Alon", compile=TRUE
)

forc = data.frame(name="k3", value = rnorm(10), times = seq(0.2,20,0.1))
forcings

osc3 <- cOde::funC(
  c(
    X = "-b1*Y - a1*X",
    Y = "b2*X - a2*Y"
  ), modelname = "osc_Alon", compile=TRUE
)


# Define server logic required to draw a histogram
server <- function(input, output) {

  output$distPlot <- renderPlot({
    par(mfrow=c(1,3))
    #k_list = c(A = input$rA, B = input$rB, k1=k1, k2=k2, k3=k3, k4=k4)
    y0 <- c(X=1,Y=1)
    #example_ODE <- ode(y=y0, times=seq(0,100, 0.1), brusselator, k_list)
    k_list <- c(b1=input$beta1, b2=input$beta2, a1=input$alpha1, a2=input$alpha2)
    example_ODE <- odeC(y=y0, times=seq(0,input$len,0.01), func=osc3, parms = k_list)
    plot(example_ODE[,1], example_ODE[,2], type='l')
    plot(example_ODE[,1], example_ODE[,3], type='l')
    plot(example_ODE[,2], example_ODE[,3], type='l')
  })
}

# Run the application
shinyApp(ui = ui, server = server)
