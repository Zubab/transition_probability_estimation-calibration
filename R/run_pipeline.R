############################################################
# Title: Bayesian Estimation + Calibration of Transition Probs
# Purpose:
#   1) Estimate a 5x5 transition matrix (rows=from, cols=to)
#      via a Bayesian Dirichlet–Multinomial model in JAGS.
#   2) Calibrate that matrix via simulated annealing (SA) so the
#      projected distribution matches a target (“real”) distribution.
#
# Outputs (in-memory):
#   - final_transition_probs (5x5): calibrated, row-stochastic matrix
#   - sim_mat (2x5): simulated distribution across cycles
#   - test (scalar): RMSE between simulated and target (“real”)
############################################################

# ---------------------------
# 0) Packages
# ---------------------------
library(rjags)         # JAGS engine
library(runjags)       # (optional helpers)
library(coda)          # MCMC handling
library(dplyr)         # data wrangling
library(stats)         # base stats
library(optimization)  # optim_sa() simulated annealing
library(R2jags)        # jags() front-end

# ---------------------------
# 1) Transition counts (5x5)
#    Rows = from-state, Cols = to-state
# ---------------------------
tt <- matrix(c(
  3905, 517, 196, 146, 137,
  255, 229, 107, 122, 102,
  83,  62,  80,  73,  73,
  60,  44,  55,  91, 118,
  71,  16,  37,  77, 209
), ncol = 5, nrow = 5, byrow = TRUE)

# ---------------------------
# 2) JAGS model (Dirichlet–Multinomial)
#    - Each row p[i,] ~ Dirichlet(alpha_i)
#    - Row counts tt[i,] ~ Multinomial(p[i,], n[i])
# ---------------------------
modelstring = "
model{
  p[1, 1:S] ~ ddirch(alpha1[1:S])
  p[2, 1:S] ~ ddirch(alpha2[1:S])
  p[3, 1:S] ~ ddirch(alpha3[1:S])
  p[4, 1:S] ~ ddirch(alpha4[1:S])
  p[5, 1:S] ~ ddirch(alpha5[1:S])

  tt[1, 1:S] ~ dmulti(p[1, 1:S], n[1])
  tt[2, 1:S] ~ dmulti(p[2, 1:S], n[2])
  tt[3, 1:S] ~ dmulti(p[3, 1:S], n[3])
  tt[4, 1:S] ~ dmulti(p[4, 1:S], n[4])
  tt[5, 1:S] ~ dmulti(p[5, 1:S], n[5])
}
"
writeLines(modelstring, con = "TEMPmodel.txt")

# ---------------------------
# 3) Data & priors for JAGS
# ---------------------------
S <- 5
n <- apply(tt, 1, sum)    # row totals
alpha1 <- c(1,1,1,1,1)
alpha2 <- c(1,1,1,1,1)
alpha3 <- c(1,1,1,1,1)
alpha4 <- c(1,1,1,1,1)
alpha5 <- c(1,1,1,1,1)

params <- c("p")
# R2jags supports this character-list style:
jags.data <- list("S", "n", "alpha1","alpha2","alpha3","alpha4","alpha5","tt")

# ---------------------------
# 4) Run JAGS (Bayesian estimation)
# ---------------------------
set.seed(100)  # reproducibility for JAGS
jagsfit <- jags(
  data = jags.data,
  parameters.to.save = params,
  model.file = "TEMPmodel.txt",
  n.chains = 3,
  n.iter   = 100000,
  n.thin   = 5,
  progress.bar = "none"
)
print(jagsfit)  # quick console summary

# Posterior mean transition matrix (5x5; row-stochastic)
final_transition_probs_jags <- jagsfit$BUGSoutput$mean$p

# ---------------------------
# 5) Target (“real”) distributions across cycles (2 x 5)
#    Row 1 = initial (%), Row 2 = target after one cycle (%)
#    (Rows already sum to 100.)
# ---------------------------
real <- matrix(c(
  86.678, 5.998, 2.750, 2.463, 2.111,
  78.175, 8.294, 4.419, 4.987, 4.125
), ncol = 5, nrow = 2, byrow = TRUE)

# ---------------------------
# 6) Calibration objective (RMSE)
#    - Vectorized 5x5 matrix → row-normalize → simulate forward
#    - Minimize RMSE vs. 'real'
# ---------------------------
objectivefuntion <- function(Sim){
  Sim_mat <- matrix(Sim, nrow = 5, ncol = 5)
  Sim_mat <- Sim_mat / (rowSums(Sim_mat) %x% t(rep(1,5)))  # enforce row-stochastic
  
  # Simulated prevalence has same shape as 'real'
  sim_data <- 0 * real
  
  # Initial cohort % distribution (must sum to 100)
  sim_data[1,] <- c(86.678, 5.998, 2.750, 2.463, 2.111)
  
  # Project forward across rows of 'real'
  j <- 1
  for (i in 2:dim(real)[[1]]) {
    sim_data[i,] <- sim_data[j,] %*% Sim_mat
    j <- i
  }
  
  sqrt(mean((as.matrix(sim_data) - as.matrix(real))^2))
}

# ---------------------------
# 7) Start matrix for calibration + bounds
#    - Start at JAGS posterior mean
#    - Bounds are 0.75x to 1.25x elementwise
# ---------------------------
parameter <- final_transition_probs_jags
upper <- 1.25 * parameter
lower <- 0.75 * parameter

# ---------------------------
# 8) Simulated annealing calibration
#    - rf length must equal number of parameters (5*5 = 25)
# ---------------------------
set.seed(123)  # reproducible SA
optimized_probs <- optim_sa(
  fun   = objectivefuntion,
  start = parameter,
  lower = lower,
  upper = upper,
  control = list(
    t0 = 1000,
    rf = rep(0.01, length(parameter)),  # length = 25
    nlimit = 1000,
    stopac = 99999,
    r = 0.8
  )
)

# ---------------------------
# 9) Final calibrated TP matrix (row-stochastic)
# ---------------------------
final_transition_probs <- matrix(optimized_probs$par, nrow = 5)
final_transition_probs <- final_transition_probs / (rowSums(final_transition_probs) %x% t(rep(1,5)))

# ---------------------------
# 10) One-step simulation check + RMSE
# ---------------------------
n_t <- 2
n_s <- 5
v_state_names <- c("Abstinent","Low","Medium","High","Very_High")
m_p <- final_transition_probs

sim <- array(NA_real_, dim = c(n_t, n_s),
             dimnames = list(cycle = 1:n_t, state = v_state_names))
# initial distribution (must match the objective)
sim[1,] <- c(86.678, 5.998, 2.750, 2.463, 2.111)

for (i in 2:n_t) {
  sim[i,] <- sim[i-1,] %*% m_p
}

sim_mat <- as.matrix(sim)
test <- sqrt(mean((as.matrix(sim) - as.matrix(real))^2))

# ---------------------------
# 11) Console summaries
# ---------------------------
cat("\nRow sums of calibrated TP (should be 1):\n")
print(round(rowSums(final_transition_probs), 5))

cat("\nRMSE after calibration: ", round(test, 5), "\n", sep = "")
