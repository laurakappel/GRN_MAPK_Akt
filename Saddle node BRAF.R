

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
       x = "BRAF (µM)",
       y = "ERK (µM)",
       color = "") +  # Legend title
  scale_color_manual(values = c("Max ERK" = "lightgreen", "Min ERK" = "darkgreen", "Threshold" = "black" )) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

