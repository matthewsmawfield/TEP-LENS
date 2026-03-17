import pytest
import numpy as np
import sys
from pathlib import Path

# Add project root to path so we can import step scripts
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

# Import specific functions
from scripts.steps.step_11_model_dependence import weighted_mean_and_sigma, exact_sign_flip_pvalue, beta_binom_tail_p
from scripts.steps.step_16_independence_tier_significance import fisher_from_ps, stouffer_from_ps

def test_weighted_mean_and_sigma():
    # Test identical weights
    values = np.array([10.0, 10.0, 10.0])
    sigmas = np.array([1.0, 1.0, 1.0])
    mean, sigma = weighted_mean_and_sigma(values, sigmas)
    assert np.isclose(mean, 10.0)
    assert np.isclose(sigma, 1.0 / np.sqrt(3))

    # Test heavily skewed weights
    values = np.array([10.0, 20.0])
    sigmas = np.array([1.0, 100.0]) # second value has huge error
    mean, sigma = weighted_mean_and_sigma(values, sigmas)
    assert np.isclose(mean, 10.0, rtol=1e-3)
    assert np.isclose(sigma, 1.0, rtol=1e-3)

def test_exact_sign_flip_pvalue():
    # If all values are positive and far from 0
    values = np.array([5.0, 6.0, 7.0])
    sigmas = np.array([1.0, 1.0, 1.0])
    # The actual R_obs is positive
    # With 3 items, there are 2^3 = 8 permutations. 
    # Only the all-positive permutation will yield a sum >= the actual sum
    # so p-value should be 1/8 = 0.125
    pval = exact_sign_flip_pvalue(values, sigmas)
    assert np.isclose(pval, 0.125)

def test_fisher_from_ps():
    # Fisher method: X^2 = -2 * sum(ln(p))
    # For p = [0.1, 0.1]
    pvals = [0.1, 0.1]
    res = fisher_from_ps(pvals)
    # chi2 = -2 * (ln(0.1) + ln(0.1)) = -2 * (-2.302 - 2.302) = 9.21
    # df = 4. 
    # 1 - chi2.cdf(9.21, 4) = 0.056
    assert np.isclose(res["p"], 0.056, atol=0.001)

def test_stouffer_from_ps():
    # Stouffer unweighted
    # Z_i = norm.ppf(1 - p_i)
    # Z = sum(Z_i) / sqrt(k)
    # For p = 0.5, Z_i = 0
    pvals = [0.5, 0.5, 0.5]
    res = stouffer_from_ps(pvals)
    assert np.isclose(res["z"], 0.0)
    assert np.isclose(res["p"], 0.5)

    # For highly significant values
    pvals = [0.02275, 0.02275] # corresponds to z ~ 2
    res = stouffer_from_ps(pvals)
    # Z_combined = (2 + 2) / sqrt(2) = 4 / 1.414 = 2.828
    assert np.isclose(res["z"], 2.828, atol=0.01)
