
# Load required libraries
library(deSolve)
library(ggplot2)
library(dplyr)
library(gridExtra)  # For arranging multiple plots in one window

# Define parameters and initial conditions
base_params <- c(BRAF = 0.005, NRAS = 0.01, Ras = 0.5, 
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

initial_conditions <- c(Raf = 1, MEK = 1.2, PI3K = 1, PIP3 = 1, ERK = 1)
times <- seq(0, 1000, by = 1)

MAPK_pathway <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Sinusoidal PTEN_base with time-dependent oscillations
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

# Function to run model and extract results for sensitivity
run_simulation <- function(param_name, variation) {
  params <- base_params
  params[param_name] <- base_params[param_name] * variation
  out <- ode(y = initial_conditions, times = times, func = MAPK_pathway, parms = params)
  out_df <- as.data.frame(out)
  out_df$variation <- paste(param_name, variation, sep = "_")
  out_df$param <- param_name
  out_df$variation_value <- variation
  return(out_df)
}

# Define list for which plots to include for each parameter
plot_inclusion <- list(
  BRAF = "ERK", NRAS = c("ERK", "PIP3"), Ras = c("ERK", "PIP3"),
  PTEN_base = "PIP3", A = "PIP3", omega = "PIP3", 
  beta_ras = "ERK", beta_nras = "ERK", beta_braf = "ERK",
  beta_raf = "ERK", beta_nras1 = "PIP3", beta_ras1 = "PIP3",
  beta_pi3k = "PIP3", beta_pten = "PIP3", beta_mek = "ERK",
  alpha_raf = "ERK", alpha_mek = "ERK", alpha_erk = "ERK",
  alpha_pi3k = "PIP3", alpha_pip3 = "PIP3"
)

# List of all parameters to vary and sensitivity percentages
sensitivity_params <- names(base_params)
variations <- c(0.7, 1.0, 1.3)  # 70%, 100%, and 130% of the base value

# Run simulations for each parameter and variation
results <- bind_rows(
  lapply(sensitivity_params, function(param) {
    lapply(variations, function(var) {
      run_simulation(param, var)
    })
  })
)

# Create and store plots for each parameter in a list
plot_list <- list()
for (param in sensitivity_params) {
  
  # Filter results for the current parameter
  param_results <- results %>% filter(param == !!param)
  
  # Check inclusion list and create the specified plots
  include_plots <- plot_inclusion[[param]]
  
  if ("ERK" %in% include_plots) {
    # Plot for ERK
    p_erk <- ggplot(param_results, aes(x = time, y = ERK, color = as.factor(variation_value))) +
      geom_line(lwd = 1) +
      labs(title = paste("ERK with varying", param),
           y = "ERK (Proliferation Pathway)", x = "Time",
           color = "Variation (%)") +
      theme_minimal() +
      scale_color_manual(values = c("blue", "black", "red"))
    plot_list[[length(plot_list) + 1]] <- p_erk
  }
  
  if ("PIP3" %in% include_plots) {
    # Plot for PIP3
    p_pip3 <- ggplot(param_results, aes(x = time, y = PIP3, color = as.factor(variation_value))) +
      geom_line(lwd = 1) +
      labs(title = paste("PIP3 with varying", param),
           y = "PIP3 (Survival Pathway)", x = "Time",
           color = "Variation (%)") +
      theme_minimal() +
      scale_color_manual(values = c("blue", "black", "red"))
    plot_list[[length(plot_list) + 1]] <- p_pip3
  }
}

# Set up grid parameters and paginate plots if needed
plots_per_page <- 4  # Adjust number of plots per page
num_pages <- ceiling(length(plot_list) / plots_per_page)

# Display plots in paginated grid layout
for (i in seq_len(num_pages)) {
  # Select a subset of plots for the current page
  start_idx <- (i - 1) * plots_per_page + 1
  end_idx <- min(i * plots_per_page, length(plot_list))
  
  # Display the current set of plots
  do.call(grid.arrange, c(plot_list[start_idx:end_idx], ncol = 2, nrow = 2))
}





dev.off()



# Function to calculate maximum value in the last time steps (steady-state values)
get_max_steady_state_value <- function(param_results, pathway) {
  # Get the last few timesteps (for steady state, consider the last 100 time points)
  last_timesteps <- tail(param_results, 100)
  
  # Return the maximum value for the specified pathway (PIP3 or ERK)
  max_value <- max(last_timesteps[[pathway]], na.rm = TRUE)
  
  return(max_value)
}





# Create a list to store the max values for each parameter and variation
max_values <- list()

# Loop over each parameter, variation, and pathway
for (param in sensitivity_params) {
  for (var in variations) {
    # Filter results for the current parameter and variation
    param_results <- results %>% filter(param == !!param, variation_value == var)
    
    # Check the inclusion list and calculate max steady-state values for PIP3 or ERK
    include_plots <- plot_inclusion[[param]]
    
    if ("ERK" %in% include_plots) {
      # Calculate max value of ERK
      max_erk <- get_max_steady_state_value(param_results, "ERK")
      max_values[[length(max_values) + 1]] <- data.frame(
        Parameter = param, 
        Variation = paste0(var * 100, "%"),  # Convert to percentage format
        Pathway = "ERK", 
        Max_Value = max_erk
      )
    }
    
    if ("PIP3" %in% include_plots) {
      # Calculate max value of PIP3
      max_pip3 <- get_max_steady_state_value(param_results, "PIP3")
      max_values[[length(max_values) + 1]] <- data.frame(
        Parameter = param, 
        Variation = paste0(var * 100, "%"),  # Convert to percentage format
        Pathway = "PIP3", 
        Max_Value = max_pip3
      )
    }
  }
}

# Combine all max values into a single data frame
max_values_df <- bind_rows(max_values)

# Display the resulting table
print(max_values_df)




















# Function to calculate the amplitude (difference between max and min) in the last time steps
get_amplitude_last_period <- function(param_results, pathway) {
  # Get the last few timesteps (for the last period, consider the last 100 time points)
  last_timesteps <- tail(param_results, 500)
  
  # Calculate the amplitude (max - min) for the specified pathway (PIP3 or ERK)
  amplitude <- max(last_timesteps[[pathway]], na.rm = TRUE) - min(last_timesteps[[pathway]], na.rm = TRUE)
  
  return(amplitude)
}

# Create a list to store the max values and amplitudes for each parameter and variation
max_values_with_amplitude <- list()

# Loop over each parameter, variation, and pathway
for (param in sensitivity_params) {
  for (var in variations) {
    # Filter results for the current parameter and variation
    param_results <- results %>% filter(param == !!param, variation_value == var)
    
    # Check the inclusion list and calculate max steady-state values and amplitude for PIP3 or ERK
    include_plots <- plot_inclusion[[param]]
    
    if ("ERK" %in% include_plots) {
      # Calculate max value of ERK
      max_erk <- get_max_steady_state_value(param_results, "ERK")
      
      # Calculate amplitude of ERK in the last 100 time steps
      amplitude_erk <- get_amplitude_last_period(param_results, "ERK")
      
      max_values_with_amplitude[[length(max_values_with_amplitude) + 1]] <- data.frame(
        Parameter = param, 
        Variation = paste0(var * 100, "%"),  # Convert to percentage format
        Pathway = "ERK", 
        Max_Value = round(max_erk,2),
        Amplitude = round(amplitude_erk,2)  # Add the amplitude
      )
    }
    
    if ("PIP3" %in% include_plots) {
      # Calculate max value of PIP3
      max_pip3 <- get_max_steady_state_value(param_results, "PIP3")
      
      # Calculate amplitude of PIP3 in the last 100 time steps
      amplitude_pip3 <- get_amplitude_last_period(param_results, "PIP3")
      
      max_values_with_amplitude[[length(max_values_with_amplitude) + 1]] <- data.frame(
        Parameter = param, 
        Variation = paste0(var * 100, "%"),  # Convert to percentage format
        Pathway = "PIP3", 
        Max_Value = round(max_pip3,2),
        Amplitude = round(amplitude_pip3,2)  # Add the amplitude
      )
    }
  }
}

# Combine all max values and amplitudes into a single data frame
max_values_with_amplitude_df <- bind_rows(max_values_with_amplitude)

# Display the resulting table
print(max_values_with_amplitude_df)


write.csv(max_values_with_amplitude_df, "max_values_with_amplitude2.csv", row.names = FALSE)


