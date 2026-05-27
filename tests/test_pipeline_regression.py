"""
TEP-LENS Pipeline Regression Tests

Validates key numerical outputs from the analysis pipeline to catch
regressions when code or data are updated. These tests run against
pre-computed results in results/outputs/.

Run with:  pytest tests/ -v
"""

import json
from pathlib import Path

import numpy as np
import pytest

RESULTS_DIR = Path(__file__).resolve().parents[1] / "results" / "outputs"
DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "raw" / "sn_lensing"


# ---------- Fixtures ----------

@pytest.fixture(scope="module")
def step02():
    with open(RESULTS_DIR / "step_02_gr_closure.json") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def step03():
    with open(RESULTS_DIR / "step_03_tep_closure.json") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def step07():
    with open(RESULTS_DIR / "step_07_observed_vs_predicted.json") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def step13():
    path = RESULTS_DIR / "step_13_bayes_model_comparison.json"
    if not path.exists():
        pytest.skip("step_13 output not available")
    with open(path) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def step35():
    path = RESULTS_DIR / "step_35_single_contrast_dominance.json"
    if not path.exists():
        pytest.skip("step_35 output not available (run pipeline first)")
    with open(path) as f:
        return json.load(f)


# ---------- Data integrity ----------

class TestDataIntegrity:
    """Verify authoritative data files are self-consistent."""

    def test_catalog_exists(self):
        assert (DATA_DIR / "lensed_sn_catalog.json").exists()

    def test_blind_predictions_exist(self):
        assert (DATA_DIR / "blind_model_predictions.json").exists()

    def test_blind_predictions_count(self):
        with open(DATA_DIR / "blind_model_predictions.json") as f:
            data = json.load(f)
        models = data["models"]
        assert len(models) == 8
        assert sum(1 for m in models if m["blind"]) == 7

    def test_observed_delay_value(self):
        with open(DATA_DIR / "blind_model_predictions.json") as f:
            data = json.load(f)
        obs = data["metadata"]["observed_delay"]
        assert obs["value"] == 376.0
        assert obs["err"] == 5.6


# ---------- Step 02: GR algebraic loop sums ----------

class TestStep02GRClosure:
    """All GR loop sums must be identically zero."""

    def test_all_loops_zero(self, step02):
        loops = step02["gr_algebraic_loop_loops"]
        for name, loop in loops.items():
            assert loop["loop_sum_gr"] == 0.0, f"Loop {name} is not zero"

    def test_five_loops_present(self, step02):
        assert len(step02["gr_algebraic_loop_loops"]) == 5


# ---------- Step 03: TEP predicted discrepancies ----------

class TestStep03TEPClosure:
    """Verify TEP predicted residuals at alpha=-0.055."""

    def test_alpha_value(self, step03):
        assert step03["alpha_tep"] == -0.055

    def test_s1_s4_sx_residual(self, step03):
        r = step03["tep_predicted_discrepancies"]["S1_S4_SX"]
        assert pytest.approx(r["tep_gr_discrepancy_days"], abs=0.1) == -14.538

    def test_s1_s2_sx_residual(self, step03):
        r = step03["tep_predicted_discrepancies"]["S1_S2_SX"]
        assert pytest.approx(r["tep_gr_discrepancy_days"], abs=0.1) == -8.492

    def test_inner_cross_residuals_small(self, step03):
        for name in ["S1_S2_S3", "S1_S2_S4", "S1_S3_S4"]:
            r = step03["tep_predicted_discrepancies"][name]
            assert abs(r["tep_gr_discrepancy_days"]) < 0.5

    def test_sx_loops_high_snr(self, step03):
        for name in ["S1_S2_SX", "S1_S4_SX"]:
            r = step03["tep_predicted_discrepancies"][name]
            assert r["snr"] > 60.0


# ---------- Step 07: Observed vs. blind-predicted ----------

class TestStep07Evidence:
    """Core evidence test regression checks."""

    def test_weighted_mean_residual_positive(self, step07):
        assert step07["weighted_mean_residual"]["R_obs_days"] > 0

    def test_weighted_mean_residual_value(self, step07):
        r = step07["weighted_mean_residual"]["R_obs_days"]
        assert pytest.approx(r, abs=1.0) == 14.6

    def test_z_from_gr_above_one(self, step07):
        z = step07["weighted_mean_residual"]["z_from_gr_null"]
        assert z > 1.0

    def test_alpha_inferred_negative(self, step07):
        alpha = step07["weighted_mean_residual"]["alpha_inferred"]
        assert alpha < 0

    def test_wilcoxon_blind_p_below_0_02(self, step07):
        p = step07["binomial_sign_test"]["p_wilcoxon_signed_rank_blind"]
        assert p < 0.02, f"Wilcoxon blind p={p} exceeds 0.02"

    def test_wilcoxon_all_p_below_0_01(self, step07):
        p = step07["binomial_sign_test"]["p_wilcoxon_signed_rank_all"]
        assert p < 0.01, f"Wilcoxon all p={p} exceeds 0.01"

    def test_sign_test_7_of_8_positive(self, step07):
        assert step07["binomial_sign_test"]["n_positive"] == 7
        assert step07["binomial_sign_test"]["n_total"] == 8

    def test_chi2_tep_better_than_gr(self, step07):
        assert step07["chi2_model_comparison"]["delta_chi2"] > 0

    def test_eight_models_present(self, step07):
        assert len(step07["per_model_results"]) == 8


# ---------- Step 13: Bayesian model comparison ----------

class TestStep13BayesComparison:
    """Hierarchical Bayesian comparison must be non-decisive."""

    def test_bayes_factor_fixed_near_unity(self, step13):
        scenarios = step13.get("scenarios", {})
        for name, scenario in scenarios.items():
            if "evidence" in scenario:
                bf = scenario["evidence"]["bf_tep_fixed_over_gr"]
                assert 0.1 < bf < 10, f"{name}: Fixed BF={bf} is decisive"

    def test_posterior_alpha_spans_zero(self, step13):
        scenarios = step13.get("scenarios", {})
        for name, scenario in scenarios.items():
            if "posterior" in scenario:
                pa = scenario["posterior"]
                assert pa["alpha_p16"] < 0 < pa["alpha_p84"]


# ---------- Step 35: Signal-energy concentration ----------

class TestStep35SingleContrast:
    """Verify single-contrast dominance metrics."""

    def test_sx_energy_fraction_above_99(self, step35):
        frac = step35.get("sx_energy_fraction", 0)
        assert frac > 0.99, f"SX energy fraction {frac} < 0.99"

    def test_effective_dof_near_two(self, step35):
        d_eff = step35.get("effective_dof_participation_ratio", 0)
        assert pytest.approx(d_eff, abs=0.5) == 2.0


# ---------- Cross-step consistency ----------

class TestCrossStepConsistency:
    """Ensure outputs from different steps are mutually consistent."""

    def test_tep_residual_matches_prediction(self, step03, step07):
        r_pred = -step03["tep_predicted_discrepancies"]["S1_S4_SX"]["tep_gr_discrepancy_days"]
        r_obs = step07["weighted_mean_residual"]["R_obs_days"]
        # They should agree within ~1 sigma
        assert abs(r_pred - r_obs) < 12.0


# ---------- Utility tests ----------

class TestUtilities:
    """Verify shared utility modules work correctly."""

    def test_safe_json_default_numpy(self):
        from scripts.utils.logger import safe_json_default
        assert safe_json_default(np.float64(1.5)) == 1.5
        assert safe_json_default(np.int64(42)) == 42

    def test_safe_json_default_raises_for_unknown(self):
        from scripts.utils.logger import safe_json_default
        with pytest.raises(TypeError):
            safe_json_default(object())
