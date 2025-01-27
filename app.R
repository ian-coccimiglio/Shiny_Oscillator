library(shiny)
library(cOde)
library(deSolve)
library(bslib)

# Feature 1, Creating a damped biological oscillator. This compiles an object which
# simulates a biological system that tends towards damped oscillation.
# Users can change the INDIVIDUAL parameters and the simulation time.
# Outcomes are changes in the phase-space diagram and the time-concentration plot.

#compile <- !file.exists(paste0("BrusselC", .Platform$dynlib.ext))
osc2 <- cOde::funC(
    c(
        X = "-b1*Y - a1*X",
        Y = "b2*X - a2*Y"
    ), modelname = "Alon_damped", compile=TRUE
)

forcings <- c("noise")
noise <- cOde::funC(
    c(
        X = "-b1*Y - a1*X",
        Y = "b2*(X+noise) - a2*Y"
    ),
    forcings = forcings,
    modelname = "damped_noise",
    compile=TRUE,
    fcontrol = "nospline",
    nGridpoints = 10
)

## Feature 2, Creating a 3-part repressilator. This compiles an object which
# simulates a 3 component system that tends towards infinite oscillations.
# Users are able to change the global degradation, creation, and concentration constants.
# Additionally, they can change the simulation length.
# As above, this produces a phase-space diagram and the time-concentration plot.
osc3 <- cOde::funC(
    c(
        # 1 and 3 are rate constants for the purposes of this demo.
        X = "beta/(1+((Z/k)^3)) - g2*X",
        Y = "beta/(1+((X/k)^3)) - g2*Y",
        Z = "beta/(1+((Y/k)^3)) - g2*Z"
    ),
    modelname = "Osc_3", compile=TRUE
)

## Feature 3, Creating the brusselator This compiles an object which
# simulates a 2 component system that tends towards infinite oscillations.
# Users are able to manipulate the different inflows and outflow rates,
# as well as concentrations of products A, and B
# As above, this produces a phase-space diagram and the time-concentration plot.
brussC <- cOde::funC(
    c(
        X = "k1*A - k2*B*X + k3*(X^2)*Y - k4*X",
        Y = "k2*B*X - k3*(X^2)*Y"
    ), modelname = "brusselC", compile=TRUE
)

# Extra, this is a repression function that allows modification of N
repress_R <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
        dX    <- beta/(1+((Z/k)^n)) - g2*X
        dY   <- beta/(1+((X/k)^n)) - g2*Y
        dZ <- beta/(1+((Y/k)^n)) - g2*Z

        return(list(c(dX, dY, dZ)))
    })
}


# Define UI for application that draws a histogram
ui <- fluidPage(
    theme = bs_theme(version = 4, bootswatch = "minty"),
    titlePanel("Biological Oscillations"),
    tabsetPanel(
        tabPanel("Damped Biological Oscillator",
            sidebarLayout(
                sidebarPanel(
                    "Parameters for X",
                    hr(),
                    sliderInput("alpha1", "Degradation X (\u03B11)", min = 0, max = 20, value = 1),
                    sliderInput("beta1", "Formation X (\u03B21)", min = 0, max = 20, value = 8),
                    "Parameters for Y",
                    hr(),
                    sliderInput("alpha2", "Degradation Y (\u03B12)", min = 0, max = 20, value = 1),
                    sliderInput("beta2", "Formation Y (\u03B22)", min = 0, max = 20, value = 8),
                    "Simulation Parameters",
                    hr(),
                    sliderInput("len", "Simulation Time", min = 1, max = 50, value = 15),
                    "Click to add a small amount of random, normal noise (mu=0, sd=1).",
                    checkboxInput("addnoise", "Add Random Noise", value = FALSE)
                ),

                mainPanel(
                         h5("Damped oscillations occur in specific circumstances"),
                         "This model shows a 2-component system, where a protein/factor X activates another protein/factor Y, but Y also represses X. These are controlled by the ordinary differential equations below - change the beta and alpha parameters to see which conditions lead to oscillations.",
                         withMathJax("$$\\text{Damped Oscillation System }$$ \
                                     $$\\frac{dX}{dt}=-\\beta_1Y - \\alpha_1X$$ \
                                     $$\\frac{dY}{dt}=\\beta_2X - \\alpha_2Y$$"),
                         plotOutput("dampOsc"),
                         htmlOutput("eigenvalue"),
                         htmlOutput("jaxeign"),
                         htmlOutput("feedback"),
                         htmlOutput("mismatch"),
                         "Depending on the above inequality, if the \u03B2 feedback strength exceeds the squared mismatch in \u03B1 degradation, damped oscillations are become likely to occur (relative to monotonic behaviour).",
                         imageOutput("damp", height="200px"),
                         br(),
                         "If a moderate amount of noise is added, a damped oscillation system will start to oscillate indefinitely",
                         uiOutput("alonlink")
                )
            )
        ),

        tabPanel("Brusselator",
             sidebarLayout(
                 sidebarPanel(
                     "Generative parameters",
                     hr(),
                     sliderInput("A", "Autocatalytic A", min = 1, max = 20, value = 1),
                     sliderInput("B", "Autocatalytic B", min = 1, max = 20, value = 3),
                     "Rate constants",
                     hr(),
                     sliderInput("k1", "Catalytic Rate A (k1)", min = 1, max = 20, value = 1),
                     sliderInput("k2", "Catalytic Rate B (k2)", min = 1, max = 20, value = 1),
                     sliderInput("k3", "XY binding rate (k3)", min = 1, max = 20, value = 1),
                     sliderInput("k4", "X degradation rate (k4)", min = 1, max = 20, value = 1),
                     "Simulation parameters",
                     hr(),
                     sliderInput("len_3", "Simulation Time", min = 1, max = 200, value = 50)
                 ),

                 # Show a plot of the generated distribution
                 mainPanel(
                     h5("The brusselator system shows a noise-free method of achieving indefinite oscillation of X and Y, by incorporating additional dynamics of new components A and B"),
                     div("$$\\text{Brusselator System }$$ \
                                     $$\\frac{dX}{dt}= k_1A - k_2B + k_3X^2Y - k_4X - \\gamma{X}$$ \
                                     $$\\frac{dY}{dt}= k_2BX - k_3X^2Y - \\gamma{Y}$$"),
                     plotOutput("Brussel"),
                     imageOutput("brussel"),
                     uiOutput("brussellink")
                 )
             )
        ),

        tabPanel("Repressilator",
                 sidebarLayout(
                     sidebarPanel(
                         "Rate constants",
                         hr(),
                         sliderInput("Beta", "Protein Generation Rate (\u03B2)", min = 0, max = 20, value = 1),
                         sliderInput("Gamma", "Protein Degradation Rate (\u03B3)", min = 0, max = 5, value = 0.1, step=0.1),
                         "Simulation parameters",
                         hr(),
                         sliderInput("k", "Protein Concentration Factor", min = 1, max = 20, value = 1),
                         sliderInput("len_2", "Simulation Length", min = 1, max = 1000, value = 300)
                     ),

                     # Show a plot of the generated distribution
                     mainPanel(
                         h5("This first repressilation system shows how oscillation can occur using a cycle of repression."),
                         br(),
                         "In this system, 3 components each repress the following one in the circuit, such that X represses Y, Y represses Z, and Z represses X.",
                         br(),
                         "If the other parameters are set correctly (betas are large, gammas are small), then constant, undamped oscillation will occur",
                         br(),
                         "Here, the exponential 'Hill cooperativity' parameter is fixed at 3, so we can observe the effect of altering the formation and degradation rates",
                         div("$$\\text{Repressilation System }$$ \
                                     $$\\frac{dX}{dt}=\\frac{\\beta}{1+(Z/k)^3} - \\gamma{X}$$ \
                                     $$\\frac{dY}{dt}=\\frac{\\beta}{1+(X/k)^3} - \\gamma{Y}$$ \
                                     $$\\frac{dZ}{dt}=\\frac{\\beta}{1+(Y/k)^3} - \\gamma{Z}$$"),
                         plotOutput("Osc3"),
                         imageOutput("repress")
                     )
                 )
        ),

        tabPanel("Repressilation Cooperativity",
                 sidebarLayout(
                     sidebarPanel(
                         "Rate constants",
                         hr(),
                         sliderInput("Beta_2", "Maximal Rate (\u03B2)", min = 0, max = 20, value = 1),
                         sliderInput("Gamma_2", "Degradation Rate (\u03B3)", min = 0, max = 5, value = 0.1, step=0.1),
                         sliderInput("N", "Cooperativity (n)", min = 0, max = 5, value = 1, step=0.1),
                         "Simulation parameters",
                         hr(),
                         sliderInput("k_2", "Protein Concentration Factor", min = 1, max = 20, value = 1),
                         sliderInput("len_4", "Simulation Length", min = 1, max = 1000, value = 300)
                     ),

                     # Show a plot of the generated distribution
                     mainPanel(
                         h5("This repressilation system is the same as before, but allows you to see the impact of cooperativity"),
                         "If you increase cooperativity (n) above 2, a system of stable oscillations occurs. It is positively impacted by the beta parameters, and negatively by the gamma parameters",
                         div("$$\\text{Repressilation System }$$ \
                                     $$\\frac{dX}{dt}=\\frac{\\beta}{1+(Z/k)^n} - \\gamma{X}$$ \
                                     $$\\frac{dY}{dt}=\\frac{\\beta}{1+(X/k)^n} - \\gamma{Y}$$ \
                                     $$\\frac{dZ}{dt}=\\frac{\\beta}{1+(Y/k)^n} - \\gamma{Z}$$"),
                         plotOutput("Rep2"),
                         uiOutput("represslink")
                     )
                 )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$damp <- renderImage({
        list(src = "Damped_Oscillator.png",
        width = 200,
        height = 200)
    }, deleteFile = FALSE)

    output$repress <- renderImage({
        list(src = "Repressilator.png",
             width = 200,
             height = 250)
    }, deleteFile = FALSE)

    output$brussel <- renderImage({
        list(src = "Brusselator.png",
             width = 250,
             height = 250)
    }, deleteFile = FALSE)
    url_alon <- a("An Introduction to Systems Biology", href="https://www.amazon.ca/Introduction-Systems-Biology-Principles-Biological-dp-1439837171/dp/1439837171")
    output$alonlink <- renderUI({
        tagList("This model was inspired from Chapters 6.1 and 6.2 of the book: ", url_alon)
    })

    url_brussel <- a("Interactive web-based simulation models with deSolve", href="https://tpetzoldt.github.io/deSolve-shiny/deSolve-shiny.html")
    output$brussellink <- renderUI({
        tagList("This model was inspired from the writing of Thomas Petzoldt on using deSolve with Shiny, here: ", url_brussel)
    })

    url_repress <- a("Blinking bacteria", href="https://biocircuits.github.io/chapters/09_repressilator.html")
    output$represslink <- renderUI({
        tagList("This model was inspired from the Python code of Michael Elowitz, here: ", url_repress)
    })

    output$dampOsc <- renderPlot({

        noisefreq = 10
        step = 0.01
        times_noise = seq(0,input$len,1/noisefreq)
        times = seq(0,input$len,step)

        y0 <- c(X=1,Y=1)
        k_list <- c(b1=input$beta1, b2=input$beta2, a1=input$alpha1, a2=input$alpha2)

        if(input$addnoise == TRUE) {
            noisy <- data.frame(name=rep("noise",
                                input$len*noisefreq+1),
                                time=times_noise,
                                value=rnorm(input$len*noisefreq+1))
        }
        else {
            noisy <- data.frame(name=rep("noise",
                                input$len*noisefreq+1),
                                time=times_noise,
                                value=rep(0, input$len*noisefreq+1))
        }


        forc <- setForcings(noise, noisy)

        damp_osc <- odeC(y=y0, times=times, func=noise, parms = k_list, forcings=forc)

        par(mfrow=c(1,2))
        plot(damp_osc[,1], damp_osc[,2], type='l',
             xlab = "Time",
             ylab = "Concentration")
        lines(damp_osc[,1], damp_osc[,3], col='red', type='l')
        legend("topright", legend=c("Concentration X", "Concentration Y"), col=c("black", "red"), lty=1, cex=0.8)
        plot(damp_osc[,2], damp_osc[,3], type='l',
             xlab = "Concentration X",
             ylab = "Concentration Y"
             )
    }
    )

    output$Osc3 <- renderPlot({
        par(mfrow=c(1,2))
        Z_Conc = 3
        Y_Conc = 2
        X_Conc = 1
        y0_2 <- c(X=X_Conc,Y=Y_Conc,Z=Z_Conc)
        k_list_2 <- c(beta=input$Beta, k=input$k, g2=input$Gamma)
        step=0.2
        osc_3 <- odeC(y=y0_2, times=seq(0,input$len_2,step), func=osc3, parms = k_list_2)
        plot(osc_3[,1], osc_3[,2], col='black', type='l',
             xlab = "Time",
             ylab = "Concentration")
        lines(osc_3[,1], osc_3[,3], col='red')
        lines(osc_3[,1], osc_3[,4], col='blue')
        legend("topright", legend=c("Concentration X", "Concentration Y", "Concentration Z"), col=c("black", "red", "blue"), lty=1, cex=0.8)
        plot(osc_3[,2], osc_3[,3], col='black', type='l',
             xlab = "Concentration X",
             ylab = "Concentration Y"
             )
    }
    )

    output$Rep2 <- renderPlot({
        par(mfrow=c(1,2))
        Z_Conc = 3
        Y_Conc = 2
        X_Conc = 1
        yini <- c(X=X_Conc,Y=Y_Conc,Z=Z_Conc)
        pars <- c(beta=input$Beta_2, k=input$k_2, n=input$N, g2=input$Gamma_2)
        step=0.2
        times <- seq(0,input$len_4,step)
        out_rep <- ode(yini, times, repress_R, pars)
        plot(out_rep[,1], out_rep[,2], col='black', type='l',
             xlab = "Time",
             ylab = "Concentration",
             ylim = c(0, max(out_rep[,2:4])))
        lines(out_rep[,1], out_rep[,3], col='red')
        lines(out_rep[,1], out_rep[,4], col='blue')
        legend("topright", legend=c("Concentration X", "Concentration Y", "Concentration Z"), col=c("black", "red", "blue"), lty=1, cex=0.8)
        plot(out_rep[,2], out_rep[,3], col='black', type='l',
             xlab = "Concentration X",
             ylab = "Concentration Y"
        )
    }
    )

    output$Brussel <- renderPlot({
        par(mfrow=c(1,2))
        X_Conc_Bruss <- 2
        Y_Conc_Bruss <- 2
        y0_3 <- c(X=X_Conc_Bruss,Y=Y_Conc_Bruss)
        k_list3 <- c(A=input$A, B=input$B, k1=input$k1, k2=input$k2, k3=input$k3, k4=input$k4)
        step=0.1
        ode_bruss <- odeC(y=y0_3, times=seq(0,input$len_3,step), func=brussC, parms = k_list3)
        height = max(ode_bruss[,2:3])
        plot(ode_bruss[,1], ode_bruss[,2],col='black', type='l',ylim=c(0,height),
             xlab = "Time",
             ylab = "Concentration")
        lines(ode_bruss[,1], ode_bruss[,3],col='blue')
        legend("topright", legend=c("Concentration X", "Concentration Y"), col=c("black", "blue"), lty=1, cex=0.8)
        plot(ode_bruss[,2], ode_bruss[,3],type='l',
             xlab = "Concentration X",
             ylab = "Concentration Y")
    }
    )

    RHS=reactive(4*(input$beta1*input$beta2))
    LHS=reactive((input$alpha1-input$alpha2)^2)
    eign <- reactive({RHS()-LHS()})

    output$jaxeign <- renderText({
        if( eign() > 0){
            return(paste(withMathJax('$$(\\alpha_1-\\alpha_2)^2 < 4\\beta_1\\beta_2$$')))
        } else{
            return(paste(withMathJax('$$(\\alpha_1-\\alpha_2)^2 > 4\\beta_1\\beta_2$$')))
        }
    }
    )

    output$feedback <- renderText({paste("Feedback Strength = (\u03B11-\u03B12)^2 = ", RHS())})
    output$mismatch <- renderText({paste("Squared mismatch = 4*(\u03B21*\u03B22) = ", LHS())})

    output$eigenvalue<- renderText({
        if(eign() > 0){
            return(paste('Damped oscillations <span style=\"color:green\">are</span> occurring'))
        }else{
            return(paste("Damped oscillations <span style=\"color:red\">are not</span> occurring"))
    }})
}

# Run the application
shinyApp(ui = ui, server = server)
