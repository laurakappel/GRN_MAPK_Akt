rm(list=ls())

library(deSolve)
library(ggplot2)


########### MAPK pathway and Akt pathway #############


# Define the parameters for the simulation
params <- c(NRAS = 0.01, Ras = 0.5, PTEN = 0.8,
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
times <- seq(0, 1000, by = 1)

# Define the ODE function
MAPK_pathway <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
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
  
  # Plot ERK for each BRAF value
  plot(out_df$time, out_df$ERK, type = "l", col = "blue", lwd = 2, main = paste("BRAF =", round(braf_values[i], 4)),
       xlab = "Time", ylab = "ERK", ylim = c(0, 3))
  lines(out_df$time, out_df$PIP3, col = "magenta", lwd = 2)
  legend("bottomright", legend = c("ERK", "PIP3"), col = c("blue", "magenta"), lty = 1, lwd = 2)
  
  # Pause between plots
  Sys.sleep(0.5)  # Adds a small delay so each plot is visible
}

# Print the amplitudes for each BRAF value
print("Amplitudes for each BRAF value:")
print(amplitudes)

# Plot the amplitude vs BRAF value
plot(braf_values, amplitudes, type = "b", col = "blue", pch = 19, main = "Amplitude of ERK vs BRAF",
     xlab = "BRAF", ylab = "Amplitude", ylim = c(0, max(amplitudes)))







########### Phase plot #############

rm(list=ls())

library(deSolve)
library(ggplot2)

########### MAPK pathway and Akt pathway #############
#Y --> X

# Define the parameters for the simulation
params <- c(NRAS = 0.01, Ras = 0.5, PTEN = 0.8,
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
times <- seq(0, 1000, by = 1)

# Define the ODE function
MAPK_pathway <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
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

# Initialize a vector to store the amplitudes for each BRAF value
braf_values <- seq(0.003, 0.01, by = 0.001)

# Loop through the different BRAF values
for (i in 1:length(braf_values)) {
  # Update the BRAF parameter
  params["BRAF"] <- braf_values[i]
  
  # Run the ODE simulation
  out <- ode(y = initial_conditions, times = times, func = MAPK_pathway, parms = params)
  
  out_df <- as.data.frame(out)
  
  #Remove the initial transient data (e.g., the first 500 time points)
  steady_state_erk <- out_df$ERK[out_df$time > 500]  # Adjust time as needed to skip transients
  steady_state_meek <- out_df$MEK[out_df$time > 500]  # MEK values after transient
  
  #Create a phase plot (ERK vs MEK)
  plot(steady_state_erk, steady_state_meek, type = "l", col = "blue", 
       xlab = "ERK", ylab = "MEK", main = paste("Phase Plot (ERK vs MEK) for BRAF =", round(braf_values[i], 4)))
  
  points(steady_state_erk, steady_state_meek, col = "red", pch = 16)
  
  # Pause between plots
  Sys.sleep(0.5)  # Adds a small delay so each plot is visible
}


# Plot MEK vs. ERK to explore phase behavior
ggplot(out_df, aes(x = MEK, y = ERK)) +
  geom_path(color = "blue") +
  labs(title = "Phase Plot of MEK vs ERK",
       x = "MEK",
       y = "ERK")




### Saddle node bifurcation in BRAF vs ERK
# Set up a range of BRAF values to test
braf_values <- seq(0.003, 0.01, by = 0.0005)

# Initialize vectors to store the bifurcation data
braf_list <- numeric()
max_erk <- numeric()
min_erk <- numeric()

# Loop over each BRAF value to perform the bifurcation analysis
for (BRAF_value in braf_values) {
  # Update BRAF in parameters
  params["BRAF"] <- BRAF_value
  
  # Run the ODE simulation
  out <- ode(y = initial_conditions, times = times, func = MAPK_pathway, parms = params)
  out_df <- as.data.frame(out)
  
  # Remove transient data (e.g., the first 500 time points) to focus on steady-state behavior
  steady_state_erk <- out_df$ERK[out_df$time > 500]
  
  # Calculate the maximum and minimum ERK values after transient removal
  max_erk_value <- max(steady_state_erk)
  min_erk_value <- min(steady_state_erk)
  
  # Store the results for the bifurcation plot
  braf_list <- c(braf_list, BRAF_value)
  max_erk <- c(max_erk, max_erk_value)
  min_erk <- c(min_erk, min_erk_value)
}

# Create a data frame to store the bifurcation results
bifurcation_df <- data.frame(BRAF = braf_list, Max_ERK = max_erk, Min_ERK = min_erk)

# Plotting the bifurcation diagram for ERK vs BRAF with threshold line
ggplot(bifurcation_df, aes(x = BRAF)) +
  geom_line(aes(y = Max_ERK, color = "Max ERK"), size = 1.2) +  # Plot max ERK as stable branch
  geom_line(aes(y = Min_ERK, color = "Min ERK"), size = 1.2) + # Plot min ERK as unstable branch
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "black", size = 1) + # Threshold line at ERK = 1.5
  labs(title = "Bifurcation Diagram for ERK vs BRAF with Threshold",
       x = "BRAF (nM)",
       y = "ERK (nM)",
       color = "") +  # Legend title
  scale_color_manual(values = c("Max ERK" = "lightgreen", "Min ERK" = "darkgreen", "Threshold" = "black" )) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
