####################### Bifurcations ##########################

rm(list=ls())

library(deSolve)
library(ggplot2)

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
#times <- seq(0, 86400, by = 10)

# Define the ODE function
MAPK_pathway <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Sinusoidal PTEN with time-dependent oscillations
    PTEN <- PTEN_base + A * sin(omega * time)
    
    # ODEs for MAPK pathway
    # dRaf <- beta_nras * NRAS + beta_ras * Ras - alpha_raf * Raf
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


###########################################################################
############################ Amplitude vs BRAF ############################ 
###########################################################################
# Initialize a vector to store the amplitudes for each BRAF value
braf_values <- seq(0.003, 0.01, by = 0.001)
amplitudes <- numeric(length(braf_values))

# Loop through the different BRAF values
for (i in 1:length(braf_values)) {
  # Update the BRAF parameter
  params["BRAF"] <- braf_values[i]
  
  # Run the ODE simulation
  out <- ode(y = initial_conditions, times = times, func = MAPK_pathway, parms = params)
  
  out_df <- as.data.frame(out)
  
  # Step 1: Remove the initial transient data (e.g., the first 500 time points)
  steady_state_erk <- out_df$ERK[out_df$time > 500]  # Adjust time as needed to skip transients
  
  # Step 2: Calculate the amplitude of the steady-state oscillations
  amplitude_erk <- max(steady_state_erk) - min(steady_state_erk)
  
  # Store the amplitude for the current BRAF value
  amplitudes[i] <- amplitude_erk
  
  # Optional: Plot ERK for each BRAF value
  plot(out_df$time, out_df$ERK, type = "l", col = "blue", lwd = 2, main = paste("BRAF =", round(braf_values[i], 4)),
       xlab = "Time", ylab = "ERK", ylim = c(0, 3))
  # lines(out_df$time, out_df$PIP3, col = "magenta", lwd = 2)
  legend("bottomright", legend = c("ERK"), col = c("blue"), lty = 1, lwd = 2)
  
  # Pause between plots (optional)
  Sys.sleep(0.5)  # Adds a small delay so each plot is visible
}

# Print the amplitudes for each BRAF value
print("Amplitudes for each BRAF value:")
print(amplitudes)

# Plot the amplitude vs BRAF value
plot(braf_values, amplitudes, type = "b", col = "blue", pch = 19, main = "Amplitude of ERK vs BRAF",
     xlab = "BRAF", ylab = "Amplitude", ylim = c(0, max(amplitudes)))

### Nice plot
# Enhanced plot of Amplitude vs BRAF value
par(mfrow = c(1, 1))
plot(braf_values, amplitudes, type = "b", col = "dodgerblue3", pch = 19, lwd = 2,
     main = "Amplitude of ERK vs BRAF",
     xlab = "BRAF", ylab = "Amplitude", ylim = c(0, max(amplitudes) * 1.1),
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1, font.main = 2)

# Add a subtle grid for better readability
grid(col = "gray85", lty = "dotted")

# Optional: Add a smoother line to emphasize trends if it suits your data
lines(smooth.spline(braf_values, amplitudes), col = "darkblue", lwd = 1.5, lty = "solid")

# Highlight data points
points(braf_values, amplitudes, col = "dodgerblue3", pch = 19, cex = 1.3)


###########################################################################
############################ Amplitude vs NRAS ############################ 
###########################################################################
# Initialize a vector to store the amplitudes for each BRAF value
nras_values <- seq(0.01, 0.08, by = 0.005)
amplitudes <- numeric(length(nras_values))

# Loop through the different BRAF values
for (i in 1:length(nras_values)) {
  # Update the NRAS parameter
  params["NRAS"] <- nras_values[i]
  
  # Run the ODE simulation
  out <- ode(y = initial_conditions, times = times, func = MAPK_pathway, parms = params)
  
  out_df <- as.data.frame(out)
  
  # Step 1: Remove the initial transient data (e.g., the first 500 time points)
  steady_state_erk <- out_df$ERK[out_df$time > 500]  # Adjust time as needed to skip transients
  
  # Step 2: Calculate the amplitude of the steady-state oscillations
  amplitude_erk <- max(steady_state_erk) - min(steady_state_erk)
  
  # Store the amplitude for the current BRAF value
  amplitudes[i] <- amplitude_erk
  
  # Optional: Plot ERK for each BRAF value
  plot(out_df$time, out_df$ERK, type = "l", col = "blue", lwd = 2, main = paste("NRAS =", round(nras_values[i], 4)),
       xlab = "Time", ylab = "ERK", ylim = c(0, 3))
  lines(out_df$time, out_df$PIP3, col = "magenta", lwd = 2)
  legend("bottomright", legend = c("ERK"), col = c("blue"), lty = 1, lwd = 2)
  
  # Pause between plots (optional)
  Sys.sleep(0.5)  # Adds a small delay so each plot is visible
}

# Print the amplitudes for each BRAF value
print("Amplitudes for each NRAS value:")
print(amplitudes)

### Nice plot
# Enhanced plot of Amplitude vs BRAF value
par(mfrow = c(1, 1))
plot(nras_values, amplitudes, type = "b", col = "forestgreen", pch = 19, lwd = 2,
     main = "Amplitude of ERK vs NRAS",
     xlab = "NRAS", ylab = "Amplitude", ylim = c(0, max(amplitudes) * 1.1),
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1, font.main = 2)

# Add a subtle grid for better readability
grid(col = "gray85", lty = "dotted")

# Optional: Add a smoother line to emphasize trends if it suits your data
lines(smooth.spline(nras_values, amplitudes), col = "darkgreen", lwd = 1.5, lty = "solid")

# Highlight data points
points(nras_values, amplitudes, col = "forestgreen", pch = 19, cex = 1.3)



###########################################################################
######################### Phase plot of ERK vs MEK ######################## 
###########################################################################

rm(list=ls())

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
#times <- seq(0, 86400, by = 10)

# Define the ODE function
MAPK_pathway <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Sinusoidal PTEN with time-dependent oscillations
    PTEN <- PTEN_base + A * sin(omega * time)
    
    # ODEs for MAPK pathway
    # dRaf <- beta_nras * NRAS + beta_ras * Ras - alpha_raf * Raf
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

# Remove the initial transient data (e.g., the first 500 time points)
steady_state_erk <- out_df$ERK[out_df$time > 500]  # Adjust time as needed to skip transients
steady_state_mek <- out_df$MEK[out_df$time > 500]  # MEK values after transient

# Create a phase plot (ERK vs MEK)
par(mfrow = c(1, 1))
plot(steady_state_erk, steady_state_mek, type = "l", col = "dodgerblue", lwd = 2,
     xlab = "ERK", ylab = "MEK", main = "Phase Plot: ERK vs MEK", 
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1, font.main = 2,
     xlim = c(min(steady_state_erk), max(steady_state_erk)),  # Adjust X-axis limits
     ylim = c(min(steady_state_mek), max(steady_state_mek))) # Adjust Y-axis limits

# Add a subtle grid for better readability
grid(col = "gray85", lty = "dotted")




###########################################################################
####################### Phase plot of PI3K vs PIP3 ######################## 
###########################################################################


