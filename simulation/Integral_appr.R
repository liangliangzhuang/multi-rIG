# ========================================================================
# ===== Multivariate inverse Gaussian process with common effects ========
# ======= Comparison of Different Integral Approximation Methods =========
# ===== Monte Carlo, empirical and Gauss-Legendre quadrature methods =====
# ========================================================================
# load packages
library(gaussquad)
library(ggplot2)
library(SuppDists)
library(tidyverse)
library(patchwork)
# load functions
source("utility/fct.R")
source("utility/appr.R")

# Simulation study =========
# Parameter settings
K <- 3 # Number of PCs
m <- 10 # Number of measurements
n <- 5  # Number of units
t <- 0:m; time = 0:m
beta <- c(0.4, 0.5, 0.6)*2
gamma <- 4

# Scen.I -- pp ============
types <- "pp"
para = c(1.0, 1.0, 0.8, 1, 1.2, beta, gamma)
# Generate simulated data
sim_dat1 <- sim.dat.path(type = types, m, K, n, t, par = para)
D = c(3.6, 4.8, 7.2) # Threshold
# Comparison of different approximation methods
pp_app = scen_appr_CDF(types=types, D=D, sim_dat=sim_dat1, para = para, leg.pos = c(0.15,0.65))
# Plot 
ppp1 = pp_app$plot + coord_cartesian(xlim = c(9.5, 9.9), ylim = c(0.6, 0.86)) + 
  theme(legend.position = 'none') +
  xlab("") + ylab("")
# Add a zoomed in figure
pp1 = pp_app$plot + 
  geom_rect(aes(xmin = 9.5, xmax = 9.9, ymin = 0.6, ymax = 0.86),
                fill = "transparent", color = "black", alpha = 0, 
                linetype = "dashed", linewidth =0.2) + # Add a box to select the zoom position (Manually adjust based on results)
  geom_segment(aes(x = 9.4, xend = 8.6, y = 0.75, yend = 0.62), 
               col = "gray60", linewidth =0.2,linetype = "dashed",
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + # Add a pointing arrow
  inset_element(ppp1, 0.25, 0.4, 0.7, 0.9, on_top = TRUE) 
  
# Scen.II -- pl ============
types <- "pl"
para = c(1.0, 1.0, 0.25, 0.33, 0.37, beta, gamma) 
# Generate simulated data
sim_dat2 <- sim.dat.path(type = types, m, K, n, t, par = para)
D = c(4, 8, 13) # Threshold
pl_app = scen_appr_CDF(types=types, D=D, sim_dat=sim_dat2, para = para, limits=c(8, 10))
ppp2 = pl_app$plot + coord_cartesian(xlim = c(9.5, 9.9), ylim = c(0.6, 1.5)) + 
  theme(legend.position = 'none') +
  xlab("") + ylab("")
# Add a zoomed in figure
pl1 = pl_app$plot + 
  geom_rect(aes(xmin = 9.5, xmax = 9.9, ymin = 0.6, ymax = 1.5),
            fill = "transparent", color = "black", alpha = 0, 
            linetype = "dashed", linewidth =0.2) + # Add a box to select the zoom position (Manually adjust based on results)
  geom_segment(aes(x = 9.45, xend = 9, y = 1, yend = 1.04), 
               col = "gray60", linewidth =0.2,linetype = "dashed",
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + # Add a pointing arrow
  inset_element(ppp2, 0.02, 0.4, 0.5, 0.9, on_top = TRUE) 

# Scen.III -- lp ============
types <- "lp"
para = c(0.05, 1.0, 0.3, 0.4, 0.4, beta, gamma) 
# Generate simulated data
sim_dat3 <- sim.dat.path(type = types, m, K, n, t, par = para)
D = c(0.44, 0.66, 0.88) # Threshold
lp_app = scen_appr_CDF(types=types, D=D, sim_dat=sim_dat3, para = para, limits=c(0, 10))
# Add a zoomed in figure
ppp3 = lp_app$plot + coord_cartesian(xlim = c(9, 9.9), ylim = c(0.83, 0.95)) + 
  theme(legend.position = 'none') +
  xlab("") + ylab("")
lp1 = lp_app$plot + 
  geom_rect(aes(xmin = 9.5, xmax = 9.9, ymin = 0.8, ymax = 0.95),
            fill = "transparent", color = "black", alpha = 0, 
            linetype = "dashed", linewidth =0.2) + # Add a box to select the zoom position (Manually adjust based on results)
  geom_segment(aes(x = 9.5, xend = 5.2, y = 0.86, yend = 0.78), 
               col = "gray60", linewidth =0.2,linetype = "dashed",
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + # Add a pointing arrow
  inset_element(ppp3, 0.02, 0.5, 0.5, 0.98, on_top = TRUE) 

# Scen.IV -- ll ============
types <- "ll"
para = c(0.05, 1.0, 0.1, 0.1, 0.1, beta, gamma)
# Generate simulated data
sim_dat4 <- sim.dat.path(type = types, m, K, n, t, par = para)
D = c(0.42, 0.56, 0.42) # Threshold
ll_app = scen_appr_CDF(types=types, D=D, sim_dat=sim_dat4, para = para, limits=c(2, 10))
# Add a zoomed in figure
ppp4 = ll_app$plot + coord_cartesian(xlim = c(8, 9), ylim = c(0.65, 0.9)) + 
  theme(legend.position = 'none') +
  xlab("") + ylab("")
ll1 = ll_app$plot + 
  geom_rect(aes(xmin = 8, xmax = 9, ymin = 0.65, ymax = 0.9),
            fill = "transparent", color = "black", alpha = 0, 
            linetype = "dashed", linewidth =0.2) + # Add a box to select the zoom position (Manually adjust based on results)
  geom_segment(aes(x = 7.9, xend = 6, y = 0.78, yend = 0.75), 
               col = "gray60", linewidth =0.2,linetype = "dashed",
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + # Add a pointing arrow
  inset_element(ppp4, 0.02, 0.5, 0.5, 0.98, on_top = TRUE) 

## Time comparison =====
cost_time_sum = cbind(pp_app$cost_time, pl_app$cost_time,lp_app$cost_time,ll_app$cost_time) 
round(apply(cost_time_sum,1,mean),3)

# Figure 2
cowplot::plot_grid(pp1, pl1, lp1, ll1,nrow=2,labels = c("I","II","III","IV"))
# ggsave("simulation/CDF_appr.pdf", width = 12, height = 8)

# save.image(file = paste("simulation/Integral_appr.RData", sep=''))



