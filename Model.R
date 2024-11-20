rm(list=ls())

library(deSolve)
library(ggplot2)

########### MAPK pathway and Akt pathway #############

# Values from articles: BRAF, 

params <- c(BRAF = 0.005, NRAS =0.01, Ras = 0.5, 
            PTEN_base = 0.8, A = 0.2, omega = 0.02, 
            beta_ras = 0.4,
            beta_nras = 0.25,
            beta_braf = 6, beta_raf = 0.1,
            beta_nras1 = 0.1, beta_ras1 = 0.15, beta_pi3k = 0.25,
            beta_pten = 0.1,
            beta_mek = 0.25,
            
            alpha_raf = 0.05, alpha_mek = 0.1, alpha_erk = 0.1,
            alpha_pi3k = 0.08, alpha_pip3 = 0.11
)

# Define the initial conditions
initial_conditions <- c(Raf = 1, 
                        MEK = 1.2, 
                        PI3K = 1, 
                        PIP3 = 1,
                        ERK = 1)

# Define the time sequence for the simulation
times <- seq(0, 1000, by = 0.1)

# Define the ODE function
MAPK_pathway <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Sinusoidal PTEN with time-dependent oscillations
    PTEN <- PTEN_base + A * sin(omega * time)
    
    # ODEs for MAPK pathway
  
    dRaf <- 1/(1+(ERK/1.3)^30)* beta_ras * Ras + beta_nras * NRAS - alpha_raf * Raf
    
    dMEK <- beta_braf * BRAF + beta_raf * Raf - alpha_mek * MEK
    
    dERK <- beta_mek * MEK - alpha_erk * ERK
    
    # ODEs for Akt pathway
    dPI3K <- beta_nras1 * NRAS + beta_ras1 * Ras - alpha_pi3k * PI3K
    
    dPIP3 <- beta_pi3k * PI3K - beta_pten * PTEN - alpha_pip3 * PIP3
    
    return(list(c(dRaf, dMEK, dPI3K, dPIP3, dERK)))
  })
}

out <- ode(y = initial_conditions, times = times, func = MAPK_pathway, parms = params)

out_df <- as.data.frame(out)


############ Nice plot #############
# Load necessary package for color palettes if not already installed
library(RColorBrewer)

# Set up colors for the lines
colors <- c("dodgerblue3", "deeppink", "black")  # Black for the threshold line

# Start plotting
par(mfrow = c(1, 1))
plot(out_df$time, out_df$ERK, type = "l", col = colors[1], lwd = 2, 
     main = "Concentration of ERK and PIP3 Over Time",
     xlab = "Time (hours)", ylab = "Concentration (µM)", 
     ylim = c(0, 3.2), 
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1, font.main = 2)

# Add a grid for better readability
grid(col = "gray85", lty = "dotted")

# Add second line for PIP3
lines(out_df$time, out_df$PIP3, col = colors[2], lwd = 2)

# Add horizontal dotted line at y = 1.3
abline(h = 1.5, col = colors[3], lty = 2, lwd = 1.5)


legend(x=750,y=1, 
       legend = c("ERK (proliferation)", "PIP3 (survival)", "Threshold"), 
       col = colors, 
       lty = c(1, 1, 2), 
       lwd = c(2, 2, 1.5), 
       bg = "transparent", 
       box.lwd = 0, 
       box.col = "transparent",
       cex = 0.9,               # Adjust text size
       inset = 0.001,            # Position closer to the edge
       y.intersp = 0.5,         # Reduce vertical spacing between lines
       x.intersp = 0.5          # Reduce horizontal padding
)


##############################################################################
########################### Mutation of BRAF #################################
##############################################################################
# BRAF = 0.01

params <- c(BRAF = 0.01, NRAS =0.01, Ras = 0.5, 
            PTEN_base = 0.8, A = 0.2, omega = 0.02, 
            beta_ras = 0.4,
            beta_nras = 0.25,
            beta_braf = 6, beta_raf = 0.1,
            beta_nras1 = 0.1, beta_ras1 = 0.15, beta_pi3k = 0.25,
            beta_pten = 0.1,
            beta_mek = 0.25,
            
            alpha_raf = 0.05, alpha_mek = 0.1, alpha_erk = 0.1,
            alpha_pi3k = 0.08, alpha_pip3 = 0.11
)

out <- ode(y = initial_conditions, times = times, func = MAPK_pathway, parms = params)

out_df <- as.data.frame(out)

par(mfrow = c(1, 1))
plot(out_df$time, out_df$ERK, type = "l", col = colors[1], lwd = 2, 
     main = "Mutation of BRAF",
     xlab = "Time (hours)", ylab = "Concentration (µM)", 
     ylim = c(0, 3.2), 
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1, font.main = 2)

# Add a grid for better readability
grid(col = "gray85", lty = "dotted")

# Add second line for PIP3
lines(out_df$time, out_df$PIP3, col = colors[2], lwd = 2)

# Add horizontal dotted line at y = 1.3
abline(h = 1.5, col = colors[3], lty = 2, lwd = 1.5)

# Add legend including the Threshold line
legend(x=750,y=1, 
       legend = c("ERK (proliferation)", "PIP3 (survival)", "Threshold"), 
       col = colors, 
       lty = c(1, 1, 2), 
       lwd = c(2, 2, 1.5), 
       bg = "transparent", 
       box.lwd = 0, 
       box.col = "transparent",
       cex = 0.9,               # Adjust text size
       inset = 0.001,            # Position closer to the edge
       y.intersp = 0.5,         # Reduce vertical spacing between lines
       x.intersp = 0.5          # Reduce horizontal padding
)

##############################################################################
########################### Mutation of NRAS #################################
##############################################################################
# NRAS = 0.08

params <- c(BRAF = 0.005, NRAS =0.08, Ras = 0.5, 
            PTEN_base = 0.8, A = 0.2, omega = 0.02, 
            beta_ras = 0.4,
            beta_nras = 0.25,
            beta_braf = 6, beta_raf = 0.1,
            beta_nras1 = 0.1, beta_ras1 = 0.15, beta_pi3k = 0.25,
            beta_pten = 0.1,
            beta_mek = 0.25,
            
            alpha_raf = 0.05, alpha_mek = 0.1, alpha_erk = 0.1,
            alpha_pi3k = 0.08, alpha_pip3 = 0.11
)

out <- ode(y = initial_conditions, times = times, func = MAPK_pathway, parms = params)

out_df <- as.data.frame(out)

par(mfrow = c(1, 1))
plot(out_df$time, out_df$ERK, type = "l", col = colors[1], lwd = 2, 
     main = "Mutation of NRAS",
     xlab = "Time (hours)", ylab = "Concentration (µM)", 
     ylim = c(0, 3.2), 
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1, font.main = 2)

# Add a grid for better readability
grid(col = "gray85", lty = "dotted")

# Add second line for PIP3
lines(out_df$time, out_df$PIP3, col = colors[2], lwd = 2)

# Add horizontal dotted line at y = 1.3
abline(h = 1.5, col = colors[3], lty = 2, lwd = 1.5)

# Add legend including the Threshold line
legend(x=750,y=1, 
       legend = c("ERK (proliferation)", "PIP3 (survival)", "Threshold"), 
       col = colors, 
       lty = c(1, 1, 2), 
       lwd = c(2, 2, 1.5), 
       bg = "transparent", 
       box.lwd = 0, 
       box.col = "transparent",
       cex = 0.9,               # Adjust text size
       inset = 0.001,            # Position closer to the edge
       y.intersp = 0.5,         # Reduce vertical spacing between lines
       x.intersp = 0.5          # Reduce horizontal padding
)


##############################################################################
################## "Mutation"/downregulation of PTEN #########################
##############################################################################
# PTEN_base = 0.5, A = 0

params <- c(BRAF = 0.005, NRAS =0.01, Ras = 0.5, 
            PTEN_base = 0.5, A = 0, omega = 0.02, 
            beta_ras = 0.4,
            beta_nras = 0.25,
            beta_braf = 6, beta_raf = 0.1,
            beta_nras1 = 0.1, beta_ras1 = 0.15, beta_pi3k = 0.25,
            beta_pten = 0.1,
            beta_mek = 0.25,
            
            alpha_raf = 0.05, alpha_mek = 0.1, alpha_erk = 0.1,
            alpha_pi3k = 0.08, alpha_pip3 = 0.11
)

out <- ode(y = initial_conditions, times = times, func = MAPK_pathway, parms = params)

out_df <- as.data.frame(out)

par(mfrow = c(1, 1))
plot(out_df$time, out_df$ERK, type = "l", col = colors[1], lwd = 2, 
     main = "Mutation/downregulation of PTEN",
     xlab = "Time (hours)", ylab = "Concentration (µM)", 
     ylim = c(0, 3.2), 
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1, font.main = 2)

# Add a grid for better readability
grid(col = "gray85", lty = "dotted")

# Add second line for PIP3
lines(out_df$time, out_df$PIP3, col = colors[2], lwd = 2)

# Add horizontal dotted line at y = 1.3
abline(h = 1.5, col = colors[3], lty = 2, lwd = 1.5)

# Add legend including the Threshold line
legend(x=750,y=1, 
       legend = c("ERK (proliferation)", "PIP3 (survival)", "Threshold"), 
       col = colors, 
       lty = c(1, 1, 2), 
       lwd = c(2, 2, 1.5), 
       bg = "transparent", 
       box.lwd = 0, 
       box.col = "transparent",
       cex = 0.9,        
       inset = 0.001,        
       y.intersp = 0.5,         
       x.intersp = 0.5          
)
