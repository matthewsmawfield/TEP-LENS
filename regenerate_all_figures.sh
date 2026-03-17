#!/bin/bash
# Regenerate all 9 manuscript figures with unified styling

cd /Users/matthewsmawfield/www/TEP-LENS

echo "Regenerating manuscript figures with unified styling..."
echo "=================================================="

# Figure 1: Route-closure residuals
echo ""
echo "Figure 1: Route-closure residuals"
python3 scripts/steps/step_04_plot_closure.py

# Figures 2-3: Observed vs predicted, residuals
echo ""
echo "Figures 2-3: Observed vs predicted, residuals"
python3 scripts/steps/step_07_observed_vs_predicted.py

# Figures 4-5: Delay vs mu, alpha inference
echo ""
echo "Figures 4-5: Delay vs mu, alpha inference"
python3 scripts/steps/step_08_new_evidence.py

# Figure 6: Model dependence
echo ""
echo "Figure 6: Model dependence"
python3 scripts/steps/step_11_model_dependence.py

# Figure 7: Significance synthesis
echo ""
echo "Figure 7: Significance synthesis"
python3 scripts/steps/step_16_independence_tier_significance.py

# Figure 8: H0 tension
echo ""
echo "Figure 8: H0 tension resolution"
python3 scripts/steps/step_10_h0_tension.py

# Figure 9: Evidence ladder (part of step_08)
echo ""
echo "Figure 9: Evidence ladder (included in step_08)"

echo ""
echo "=================================================="
echo "Figure regeneration complete!"
echo "Check results/figures/ for updated PNG files"
