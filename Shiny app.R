library(shiny)
library(deSolve)
library(RColorBrewer)

# Define the MAPK pathway function
MAPK_pathway <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Sinusoidal PTEN with time-dependent oscillations
    PTEN <- PTEN_base + A * sin(omega * time)
    
    # ODEs for MAPK pathway
    dRaf <- 1 / (1 + (ERK / 1.3)^30) * beta_ras * Ras + beta_nras * NRAS - alpha_raf * Raf
    dMEK <- beta_braf * BRAF + beta_raf * Raf - alpha_mek * MEK
    dERK <- beta_mek * MEK - alpha_erk * ERK
    
    # ODEs for Akt pathway
    dPI3K <- beta_nras1 * NRAS + beta_ras1 * Ras - alpha_pi3k * PI3K
    dPIP3 <- beta_pi3k * PI3K - beta_pten * PTEN - alpha_pip3 * PIP3
    
    return(list(c(dRaf, dMEK, dPI3K, dPIP3, dERK)))
  })
}

# Define the initial conditions and time sequence
initial_conditions <- c(Raf = 1, MEK = 1.2, PI3K = 1, PIP3 = 1, ERK = 1)
times <- seq(0, 1000, by = 0.1)

# Default parameters
default_params <- c(
  BRAF = 0.005, NRAS = 0.01, Ras = 0.5, PTEN_base = 0.8, A = 0.2, omega = 0.02, 
  beta_ras = 0.4, beta_nras = 0.25, beta_braf = 6, beta_raf = 0.1,
  beta_nras1 = 0.1, beta_ras1 = 0.15, beta_pi3k = 0.25, beta_pten = 0.1, 
  beta_mek = 0.25, alpha_raf = 0.05, alpha_mek = 0.1, alpha_erk = 0.1, 
  alpha_pi3k = 0.08, alpha_pip3 = 0.11
)

# Shiny UI
ui <- fluidPage(
  titlePanel("MAPK/Akt Pathway Simulation"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("BRAF", "BRAF:", min = 0.003, max = 0.01, value = 0.005, step = 0.001),
      sliderInput("NRAS", "NRAS:", min = 0.005, max = 0.1, value = 0.01, step = 0.005),
      sliderInput("PTEN_base", "PTEN_base:", min = 0.1, max = 1, value = 0.8, step = 0.05),
      sliderInput("A", "A (oscillation amplitude):", min = 0, max = 0.3, value = 0.2, step = 0.01)
    ),
    mainPanel(
      plotOutput("plot")
    )
  )
)

# Shiny Server
server <- function(input, output) {
  output$plot <- renderPlot({
    # Update parameters based on slider inputs
    params <- default_params
    params["BRAF"] <- input$BRAF
    params["NRAS"] <- input$NRAS
    params["PTEN_base"] <- input$PTEN_base
    params["A"] <- input$A
    
    # Solve the ODEs
    out <- ode(y = initial_conditions, times = times, func = MAPK_pathway, parms = params)
    out_df <- as.data.frame(out)
    
    # Plot the results
    colors <- c("dodgerblue3", "deeppink", "black")  # Colors for the plot
    plot(out_df$time, out_df$ERK, type = "l", col = colors[1], lwd = 2, 
         main = "Concentration of ERK and PIP3 Over Time",
         xlab = "Time (hours)", ylab = "Concentration (microM)", 
         ylim = c(0, 3.2), 
         cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1, font.main = 2)
    grid(col = "gray85", lty = "dotted")
    lines(out_df$time, out_df$PIP3, col = colors[2], lwd = 2)
    abline(h = 1.5, col = colors[3], lty = 2, lwd = 1.5)
    
    legend("bottomright", 
           legend = c("ERK (proliferation)", "PIP3 (survival)", "Threshold"), 
           col = colors, 
           lty = c(1, 1, 2), 
           lwd = c(2, 2, 1.5), 
           bg = "transparent", 
           box.lwd = 0, 
           box.col = "transparent",
           cex = 0.9, 
           inset = 0.001, 
           y.intersp = 0.9, 
           x.intersp = 0.5)
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
