import os
import re

files_to_fix = [
    'step_04_plot_closure.py',
    'step_05_tdcosmo_shear.py',
    'step_06_alpha_sensitivity.py',
    'step_07_observed_vs_predicted.py',
    'step_08_new_evidence.py',
    'step_09_precision_roadmap.py',
    'step_10_h0_tension.py',
    'step_11_model_dependence.py',
    'step_12_microlensing_robustness.py',
    'step_13_bayes_model_comparison.py',
    'step_14_external_chain_ingestion.py',
    'step_15_external_informed_inflation.py',
    'step_16_independence_tier_significance.py',
    'step_17_directional_odds.py'
]

base_dir = '/Users/matthewsmawfield/www/TEP-LENS/scripts/steps'

for f_name in files_to_fix:
    f_path = os.path.join(base_dir, f_name)
    if not os.path.exists(f_path):
        continue
        
    with open(f_path, 'r') as f:
        content = f.read()
    
    # Ensure FIG_SIZE is imported
    if 'FIG_SIZE' not in content and 'import matplotlib' in content:
        content = content.replace('set_pub_style, COLORS', 'set_pub_style, COLORS, FIG_SIZE')
        if 'FIG_SIZE' not in content:
            content = content.replace('set_pub_style\n', 'set_pub_style, COLORS, FIG_SIZE\n')
             
    # Replace figsize=(...) with figsize=FIG_SIZE
    # Except if it's a 1,2 subplot, which might need wider like (11, 4.5) - we'll make them (12, 6)
    def repl_figsize(match):
        w, h = float(match.group(1)), float(match.group(2))
        if w > 10.5: # Wide subplots keep their width
            return f'figsize=({w}, 6)'
        return 'figsize=FIG_SIZE'
        
    content = re.sub(r'figsize=\(\s*([0-9\.]+)\s*,\s*([0-9\.]+)\s*\)', repl_figsize, content)
    
    # Fix old hardcoded colors
    content = content.replace('\"#d62728\"', 'COLORS[\"red\"]')
    content = content.replace('\"#d73027\"', 'COLORS[\"red\"]')
    content = content.replace('\"#4c78a8\"', 'COLORS[\"primary\"]')
    content = content.replace('\"#2166ac\"', 'COLORS[\"accent\"]')
    content = content.replace('\"#fc8d59\"', 'COLORS[\"success\"]')
    content = content.replace('\"darkorange\"', 'COLORS[\"red\"]')
    content = content.replace('\"#d9d9d9\"', 'COLORS[\"light_gray\"]')
    content = content.replace('\"#91bfdb\"', 'COLORS[\"success\"]')
    content = content.replace('\"purple\"', 'COLORS[\"primary\"]')
    content = content.replace('\"#54a24b\"', 'COLORS[\"info\"]')
    
    with open(f_path, 'w') as f:
        f.write(content)
    print(f'Fixed {f_name}')

print("\nDone applying figsize and color fixes!")
