import json

outputs_dir = '/Users/matthewsmawfield/www/TEP-LENS/results/outputs/'

def load(f):
    with open(outputs_dir + f) as fh:
        return json.load(fh)

step_07 = load('step_07_observed_vs_predicted.json')
step_08 = load('step_08_new_evidence.json')

print("Step 07 weighted mean residual alpha:")
print(f"Alpha inferred: {step_07['weighted_mean_residual']['alpha_inferred']:.3f} +- {step_07['weighted_mean_residual']['alpha_inferred_err']:.3f}")

print("\nStep 08 test C alpha inference:")
print(f"Alpha inferred: {step_08['test_C_alpha_inference']['weighted_mean_alpha']:.3f} +- {step_08['test_C_alpha_inference']['sigma_alpha']:.3f}")

