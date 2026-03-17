import numpy as np

alpha = 0.05
mu_S1 = 1.18
mu_SX = 0.35

Gamma_S1 = 1 + alpha * np.log10(mu_S1)
Gamma_SX = 1 + alpha * np.log10(mu_SX)

print(f"Gamma_S1: {Gamma_S1:.4f}")
print(f"Gamma_SX: {Gamma_SX:.4f}")

# Delay is dt = t_SX - t_S1 (positive, ~376 days)
# t_obs = t_geom * Gamma
# dt_obs = t_SX_geom * Gamma_SX - t_S1_geom * Gamma_S1
# Let's say t_S1_geom = 100 days (absolute delay from source), t_SX_geom = 476 days.
t_S1_geom = 100
t_SX_geom = 476
dt_geom = t_SX_geom - t_S1_geom
print(f"dt_geom: {dt_geom} days")

dt_obs = t_SX_geom * Gamma_SX - t_S1_geom * Gamma_S1
print(f"dt_obs: {dt_obs:.2f} days")

print(f"TEP effect on dt: {dt_obs - dt_geom:.2f} days")
