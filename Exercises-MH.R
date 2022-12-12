# PACKAGES ----
library("tidyverse")

# PLOTTING COLOUR PALETTE ----
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")

# FUNCTIONS ----

## interaction coefficient
eta_coef <- function(f1, f2, eta) {
  return(exp(log(f1) * log(f2) * eta)) # strictly interactive effects
}

## interactions
eta_out <- function(f1, f2, eta) {
  return(f1 * f2 * eta_coef(f1, f2, eta)) # interactive eff. mult. with noninteractive
}

## caclulating equilibria

### n1
n1_eq <- function(a = 1, b1 = 1, b2 = 0.1, # species parameters
                  f11 = 1, f12 = 1, # stressor effects on resource
                  f21 = 1, f22 = 1, # stressor effects on consumer
                  eta1 = 0, eta2 = 0) { # interactive effects
  
  return(b2 / (a * eta_out(f21, f22, eta2)))
  
}

### n2
n2_eq <- function(a = 1, b1 = 1, b2 = 0.1, # species parameters
                  f11 = 1, f12 = 1, # stressor effects on resource
                  f21 = 1, f22 = 1, # stressor effects on consumer
                  eta1 = 0, eta2 = 0) { # interactive effects
  
  # numerator
  num <- a * eta_out(f21, f22, eta2) * b1 * eta_out(f11, f12, eta1) - b2
  
  # denominator
  den <- (a * eta_out(f21, f22, eta2))^2
  
  return(num / den)
}

## return both in a tibble
compute_equilibria <- function(a = 1, b1 = 1, b2 = 0.1, # species parameters
                               f11 = 1, f12 = 1, # stressor effects on resource
                               f21 = 1, f22 = 1, # stressor effects on consumer
                               eta1 = 0, eta2 = 0) { # interactive effects
  
  return(tibble(n1 = n1_eq(a, b1, b2, f11, f12, f21, f22, eta1, eta2), 
                n2 = n2_eq(a, b1, b2, f11, f12, f21, f22, eta1, eta2)))
}

## to get rho when there are interactions
rhofunc <- function(a = 1, b1 = 1, b2 = 0.1, 
                    f11 = 0.5, f12 = 0.5, 
                    f21 = 0.5, f22 = 0.5, 
                    eta1 = 10, eta2 = 10) {
  
  # numerator
  num <- (b2 - a * eta_out(f21, f22, eta2) * b1 * eta_out(f11, f12, eta1)) * (b2 - a * b1)
  
  # denominator
  den <- eta_coef(f21, f22, eta2) * (b2 - a * b1 * f11 * f21) * (b2 - a * b1 * f12 * f22)
  
  return(num / den)
}

# ANALYSES ----

## BASIC PLOT: How does rho behave at various scenarios of who's affected by what ----
runs <- 100

### creating data showing the zones ----
scenarios <- tibble(start_x = c(-1, 0, 0, -1, -0.02, -1), # bottom left of rectangle
                    end_x = c(0, 1, 1, 0, 0.02, 1), # bottom right of rectangle
                    start_y = c(-1, 0, -1, 0, -1, -0.02), # top left of rectangle
                    end_y = c(0, 1, 0, 1, 1, 0.02), # top right of rectangle
                    interaction = rep(c("antagonism", "synergism", "additive"), each = 2))

### creating data showing the line ----
scen_2 <- tibble(f_1 = 10^seq(-1, 1, len = runs)) %>% # eff. on resource species
  mutate(f_2 = (2 - f_1)) %>% # eff. on consumer species
  filter(f_2 >= 0.1) %>%
  mutate(f_1 = log10(f_1),
         f_2 = log10(f_2))

### plotting ----
ggplot(scenarios) + 
 geom_rect(aes(xmin = start_x, xmax = end_x, ymin = start_y, ymax = end_y, # zones
                fill = interaction)) + 
  geom_abline(intercept = 0, slope = 1, lty = "dashed") + # identity line
  geom_line(data = scen_2, aes(x = f_1, y = f_2), lty = "solid") + # scenario
  geom_point(data = NULL, aes(x = 0, y = 0), col = "black", cex = 4) + # centre
  geom_label(
    data = data.frame(x = c(-1, -1, -0.5, 0.5, -0.5, 0.5, -0.2), # labels
                      y = c(-1, 0.25, 0.85, -0.85, -0.85, 0.85, 0),
                      label = c(1, 2, 3, 3, 3, 3, 4)),
    aes(x = x, y = y, label = label)) +
  scale_fill_manual(values = cbPalette[c(6,7,5)]) + # color scheme
  coord_cartesian(xlim = c(-1,1), ylim = c(-1,1)) + # x and y limits
  labs(x = expression(paste(log[10],"(",f[1],")")), # labels
       y = expression(paste(log[10],"(",f[2],")")),
       fill = "") +
  theme_classic() + # plot appearance
  theme(legend.position = "bottom")

### save plot ----
ggsave(filename = "figures/plot.pdf", width = 3.5, height = 4, device = "pdf")

## SUPP INFO: What about interactions? ----
int_data <- expand_grid(
  f11 = 10^runif(10,-0.1,0.1), # independent stressor effects
  f12 = 10^runif(10,-0.1,0.1), 
  f21 = 10^runif(10,-0.1,0.1), 
  f22 = 10^runif(10,-0.1,0.1),
  b1 = c(0.1, 1, 10), # resource growth rate
  b2 = 0.1, # consumer mortality rate
  a = c(0.1, 1, 10), # consumption rate
  eta = c(-10, 0, 10)) %>% # stressor interactions, same for both species
  # computing ratios of N2:N1 when stressors are present/absent
  ## eq. when both stressors present
  mutate(n_hat_mix = compute_equilibria(a, b1, b2, f11, f12, f21, f22, eta, eta)) %>% 
  unnest(n_hat_mix, names_sep = "") %>%
  ## eq. when first stressor present
  mutate(n_hat_1 = compute_equilibria(a, b1, b2, f11, f12 = 1, f21, f22 = 1)) %>% 
  unnest(n_hat_1, names_sep = "") %>%
  # eq. when second stressor present
  mutate(n_hat_2 = compute_equilibria(a, b1, b2, f11 = 1, f12, f21 = 1, f22)) %>%
  unnest(n_hat_2, names_sep = "") %>%
  # eq. when no stressors present i.e. control
  mutate(n_hat_0 = compute_equilibria(a, b1, b2, f11 = 1, f12 = 1, f21 = 1, f22 = 1)) %>% 
  unnest(n_hat_0, names_sep = "") %>%
  filter( # only keep parameter settings for which all equilibria are >0
    n_hat_mixn1 > 0, n_hat_mixn2 > 0, n_hat_1n1 > 0, n_hat_1n2 > 0, 
         n_hat_2n1 > 0, n_hat_2n2 > 0, n_hat_0n1 > 0, n_hat_0n2 > 0) %>%
  #apply rho with interactions
  mutate(rho_int = rhofunc(a, b1, b2, f11, f12, f21, f22, eta1 = eta, eta2 = eta)) %>%
  #apply rho without interactions
  mutate(rho_ref = rhofunc(a, b1, b2, f11, f12, f21, f22, eta1 = 0, eta2 = 0)) 

### plot ----
ggplot(int_data, aes(x = rho_ref, y = rho_int, col = as.factor(eta))) + 
  geom_point() +
  geom_abline(slope = 1) + 
  scale_colour_manual(values = rev(cbPalette)) + 
  scale_x_log10() +
  scale_y_log10() +
  labs(x = expression(paste("log(", rho[0], ")")),
       y = expression(paste("log(", rho, ")")),
       col = expression(paste(eta))) + 
  theme_bw()

### save ----
ggsave("figures/interactions.pdf", width = 4, height = 2, device = "pdf")
