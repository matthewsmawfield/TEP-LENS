"""Central TEP-LENS configuration constants.

All pipeline scripts should import values from here rather than hardcoding
them. This module is the single source of truth for empirical constants.
"""

# Phenomenological lensing-sector proxy parameter (formerly ALPHA_LENS).
# This is the sensitivity of the log-magnification proxy ansatz to the
# observed blind-prediction residuals. It is NOT a fundamental coupling
# constant — no solved TEP lensing transfer function exists yet.
ALPHA_PROXY = -0.055
SIGMA_ALPHA_PROXY = 0.044

# H0 literature values
# Refsdal: Kelly+2023 TD-only measurement
H0_REFSDAL_GR = 66.6
H0_REFSDAL_ERR_PLUS = 4.1
H0_REFSDAL_ERR_MINUS = 3.3

# Planck 2018 TTTEEE+lowE+lensing (final)
H0_PLANCK = 67.4
H0_PLANCK_ERR_PLUS = 0.5
H0_PLANCK_ERR_MINUS = 0.5

# H0pe & Encore measured GR H0 (Pierel+2024 TD-only; Pierel+2026 TD-only)
# These are the low-H0 cluster values, NOT Planck.
H0_H0PE_ENC_GR = 60.9
H0_H0PE_ENC_ERR_PLUS = 5.1
H0_H0PE_ENC_ERR_MINUS = 4.6

# Step-03 S1-S4-SX loop uncertainty (days)
SIGMA_R_TEP_STEP03 = 0.21
