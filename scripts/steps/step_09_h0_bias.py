import numpy as np

# We have pairs of images. For a pair (i, j), dt_obs = t_j - t_i
# H0 is inversely proportional to the time delay: H0 \propto 1 / dt_geom
# If we observe dt_obs and assume it is dt_geom, we infer H0_inferred \propto 1 / dt_obs
# But under TEP, dt_obs = t_j_geom * Gamma_j - t_i_geom * Gamma_i
# Let's define the fractional bias in H0:
# H0_inferred / H0_true = dt_geom / dt_obs

def calculate_h0_bias(t_ref_geom, dt_geom, mu_ref, mu_i, alpha=-0.05):
    # t_ref is the arrival time of the reference image (must be large, e.g. 100 years = 36500 days)
    # This is the absolute geometric arrival time from the Big Bang / emission?
    # NO! t_geom is the transit time from source to observer!
    # Transit time is distance / c. D ~ billions of light years!
    # t_geom ~ 10^10 years ~ 3.65e12 days.
    pass

