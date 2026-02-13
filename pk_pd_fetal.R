# ==========================================
# PKPD FETAL PBPK — INTEGRATED ANALYSIS
# ==========================================

library(deSolve)
library(ggplot2)
library(dplyr)
library(purrr)
library(scales)
library(tidyr)

# ------------------------------------------
# 1. CORE MODEL FUNCTION
# ------------------------------------------
# Incorporates filtration, secretion, swallowing, and GI absorption
mf_pbpk_full <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # GA-dependent maturation functions
    fGFR <- GFR_term * exp(k_GFR * (GA - 40))
    if (fGFR < 0) fGFR <- 0
    Q_sw_ga <- Q_swallow_term * ((GA/40)^2.5)
    
    # Concentrations
    C_mat <- A_mat / V_mat
    C_fet <- A_fet / V_fet
    C_amn <- A_amn / V_amn
    
    # Total Fetal Renal Clearance (Filtration + Active Secretion)
    # Adding CL_sec_f addresses under-prediction in standard literature models
    Total_fCL <- (fu_f * fGFR) + CL_sec_f
    
    # ODE Mass Balance
    dA_mat <- -CL_mat * C_mat - Q_plac * (C_mat - C_fet)
    dA_fet <- Q_plac * (C_mat - C_fet) - Total_fCL * C_fet + Ka * A_GI
    dA_amn <- Total_fCL * C_fet - Q_sw_ga * C_amn
    dA_GI  <- Q_sw_ga * C_amn - Ka * A_GI
    
    list(c(dA_mat, dA_fet, dA_amn, dA_GI))
  })
}

# ------------------------------------------
# 2. PARAMETERS (Cefazolin Term Pregnancy)
# ------------------------------------------
params_base <- list(
  V_mat = 3.5, CL_mat = 8.5, Q_plac = 0.5, 
  V_fet = 0.36, V_amn = 0.8,
  fu_f = 0.8, GFR_term = 0.00047, k_GFR = 0.15, 
  GA = 40, Q_swallow_term = 0.23/24, Ka = 2,
  CL_sec_f = 0.0025 # Estimated secretion factor
)

state0 <- c(A_mat = 1000, A_fet = 0, A_amn = 0, A_GI = 0)
times <- seq(0, 48, by = 0.1)

# ------------------------------------------
# 3. BASE SIMULATION: PASSIVE VS. ACTIVE
# ------------------------------------------
# Scenario 1: Passive Only (Original Gap)
out_passive <- ode(y=state0, times=times, func=mf_pbpk_full, parms=modifyList(params_base, list(CL_sec_f=0))) %>% as.data.frame()

# Scenario 2: Active Secretion Included (Refitted)
out_active <- ode(y=state0, times=times, func=mf_pbpk_full, parms=params_base) %>% as.data.frame()

# Visualization of the under-prediction gap
plot_gap <- data.frame(
  time = out_passive$time,
  Passive = out_passive$A_amn / params_base$V_amn,
  Active = out_active$A_amn / params_base$V_amn
) %>% pivot_longer(-time)

ggplot(plot_gap, aes(x=time, y=value, color=name)) +
  geom_line(size=1.2) + scale_y_log10(labels = scientific) +
  theme_bw() + labs(title="Amniotic Fluid: Resolving the Under-prediction Gap", 
                    subtitle="Active secretion accounts for high observed levels in CZ",
                    y="Concentration (mg/L)", color="Model Type")


# ------------------------------------------
# 4. GESTATIONAL AGE (GA) SENSITIVITY
# ------------------------------------------
ga_range <- seq(25, 40, by = 5)
ga_sims <- map_dfr(ga_range, function(g){
  p <- modifyList(params_base, list(GA = g))
  ode(y=state0, times=times, func=mf_pbpk_full, parms=p) %>% as.data.frame() %>% mutate(GA = g)
})

ggplot(ga_sims, aes(x=time, y=A_amn/params_base$V_amn, color=factor(GA))) +
  geom_line() + scale_y_log10() + theme_bw() +
  labs(title="Impact of Gestational Age on Amniotic Accumulation", 
       y="Amniotic Conc (mg/L)", color="GA (Weeks)")

# ------------------------------------------
# 5. POPULATION VARIABILITY (PPI)
# ------------------------------------------
set.seed(42)
pop_sims <- map_dfr(1:100, function(i){
  p_rand <- params_base
  p_rand$CL_sec_f <- rnorm(1, 0.0025, 0.0005) # Variations in transporter ontogeny
  ode(y=state0, times=times, func=mf_pbpk_full, parms=p_rand) %>% as.data.frame()
})

ci_data <- pop_sims %>% group_by(time) %>%
  summarise(p05 = quantile(A_amn/0.8, 0.05), p50 = median(A_amn/0.8), p95 = quantile(A_amn/0.8, 0.95))

ggplot(ci_data, aes(x=time)) +
  geom_ribbon(aes(ymin=p05, ymax=p95), fill="blue", alpha=0.2) +
  geom_line(aes(y=p50), size=1.2) + scale_y_log10() + theme_bw() +
  labs(title="90% Population Prediction Interval", y="Amniotic Conc (mg/L)")


# ------------------------------------------
# 6. FINAL AUC & RATIO ANALYSIS
# ------------------------------------------
calc_auc <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

final_out <- out_active %>% mutate(C_mat = A_mat/3.5, C_fet = A_fet/0.36, C_amn = A_amn/0.8)

auc_mat <- calc_auc(final_out$time, final_out$C_mat)
auc_fet <- calc_auc(final_out$time, final_out$C_fet)
auc_amn <- calc_auc(final_out$time, final_out$C_amn)

auc_summary <- data.frame(
  Compartment = c("Maternal", "Fetal", "Amniotic"),
  AUC_mg_h_L = c(auc_mat, auc_fet, auc_amn),
  Ratio_to_Maternal = c(1, auc_fet/auc_mat, auc_amn/auc_mat)
)

print(auc_summary)