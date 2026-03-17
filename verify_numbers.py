import json

outputs_dir = '/Users/matthewsmawfield/www/TEP-LENS/results/outputs/'

def load(f):
    with open(outputs_dir + f) as fh:
        return json.load(fh)

step_07 = load('step_07_observed_vs_predicted.json')
step_08 = load('step_08_new_evidence.json')
step_10 = load('step_10_h0_tension.json')

print("--- Step 07: Observed vs Predicted ---")
print(f"wmean_obs_minus_model: {step_07['weighted_mean_residual']['R_obs_days']:.2f} +- {step_07['weighted_mean_residual']['sigma_days']:.2f}")
print(f"wrms_reduction: {step_07['tep_corrected_consistency']['improvement_pct']:.1f}%")

print("\n--- Step 08: New Evidence ---")
print(f"Pearson r: {step_08['test_A_delay_mu_correlation']['pearson_r']:.3f}, p: {step_08['test_A_delay_mu_correlation']['pearson_p_onesided']:.3f}")
print(f"Alpha inferred mean: {step_08['test_C_alpha_inference']['weighted_mean_alpha']:.3f} +- {step_08['test_C_alpha_inference']['sigma_alpha']:.3f}")

print("\n--- Step 10: H0 Tension ---")
for sys in ['sn_refsdal', 'sn_encore', 'sn_h0pe']:
    print(f"{sys}: {step_10[sys]['gr']:.1f} -> {step_10[sys]['tep']:.1f}")

