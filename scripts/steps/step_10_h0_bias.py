import numpy as np
import json

# TEP alpha
alpha = -0.05

def gamma(mu, mu_mean):
    return 1.0 + alpha * np.log10(mu / mu_mean)

# 1. SN Refsdal
# H0 measured: 66.6 +4.1 -3.3 (Kelly+2023)
# Delays are roughly: SX - S1 = 376 days
# mu: S1=1.158, S2=0.887, S3=0.716, S4=1.793, SX=0.347
# mean_mu = 0.980
mu_refsdal = {"S1": 1.158, "S2": 0.887, "S3": 0.716, "S4": 1.793, "SX": 0.347}
mean_mu_refsdal = np.mean(list(mu_refsdal.values()))
# GR expected delay for SX-S1 (from models) is around 363 days
# observed is 376 days.
# So observed delay is ~3.5% LARGER than GR delay.
# Thus H0 inferred is ~3.5% SMALLER than true H0.
# If H0_inferred = 66.6, then H0_true = 66.6 * (376 / 362.8) = 69.0 ?

# 2. SN H0pe
# H0 measured: 75.4 +8.1 -5.5
# Delays: A-B = -116.6, C-B = -48.6. So B-A = 116.6 days.
# mu: A=5.4, B=2.5, C=2.0
mu_h0pe = {"A": 5.4, "B": 2.5, "C": 2.0}
mean_mu_h0pe = np.mean(list(mu_h0pe.values()))
gamma_A = gamma(mu_h0pe["A"], mean_mu_h0pe)
gamma_B = gamma(mu_h0pe["B"], mean_mu_h0pe)
gamma_C = gamma(mu_h0pe["C"], mean_mu_h0pe)

# If B-A is the main delay determining H0:
# observed B-A = 116.6 days.
# Under TEP: dt_obs = Gamma_B * t_B - Gamma_A * t_A
# Since t_B = t_A + dt_GR
# dt_obs = Gamma_B * (t_A + dt_GR) - Gamma_A * t_A
# dt_obs = Gamma_B * dt_GR + t_A * (Gamma_B - Gamma_A)
print(f"H0pe: Gamma_A={gamma_A:.4f}, Gamma_B={gamma_B:.4f}, Gamma_C={gamma_C:.4f}")

# 3. SN Encore
# H0 measured: 66.9 +11.2 -8.1
# Delays: Pierel+2024 measures ~39.8 days.
# mu: What are the magnifications for Encore?
