rm(list=ls())

library(deSolve)
library(ggplot2)


########### MAPK pathway and Akt pathway #############
#Y --> X

params <- c(BRAF = 0.005, NRAS =0.01, Ras = 0.5, PTEN = 0.8,
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



par(mfrow = c(1,1))
#plot(out_df$time, out_df$MEK, type = "l", col = "blue",lwd=2, main = "Plot", xlab = "Time", ylab = "gene", ylim=c(0,4))   
plot(out_df$time, out_df$ERK, type = "l", col = "blue",lwd=2, main = "Plot", 
     xlab = "Time", ylab = "gene", ylim=c(0,3))
#lines(out_df$time, out_df$PI3K, col = "green",lwd=1)
lines(out_df$time, out_df$PIP3, col = "magenta",lwd=2)
#lines(out_df$time, out_df$Raf, col = "orange",lwd=1)
#lines(out_df$time, out_df$ERK, col = "red",lwd=2)
legend("bottomright", legend=c("ERK (proliferation)","PiP3 (survival)"), col = c("blue","magenta"), lty = 1,lwd=2)

# Plot MEK vs. ERK to explore phase behavior
ggplot(out_df, aes(x = MEK, y = ERK)) +
  geom_path(color = "blue") +
  labs(title = "Phase Plot of MEK vs ERK",
       x = "MEK",
       y = "ERK")

# Plot MEK vs. ERK to explore phase behavior
ggplot(out_df, aes(x = PIP3, y = PI3K)) +
  geom_path(color = "blue") +
  labs(title = "Phase Plot of MEK vs ERK",
       x = "PIP3",
       y = "PI3K")
