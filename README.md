# Transition Probabilities (TPs) in R — Alcohol Consumption Example

This repository shows how to **calculate transition probabilities (TPs)** with a simple example based on **alcohol consumption states**. It:
1) Estimates a 5×5 TP matrix using a **Bayesian Dirichlet–Multinomial** model in **JAGS**.
2) **Calibrates** the matrix via **simulated annealing** so the projected distribution matches a target (“real”) distribution.

## Requirements
- **JAGS** installed: https://mcmc-jags.sourceforge.io/
- R packages: `rjags`, `R2jags`, `runjags`, `coda`, `dplyr`, `optimization`

Install packages:
```r
install.packages(c("rjags","R2jags","runjags","coda","dplyr","optimization"))
