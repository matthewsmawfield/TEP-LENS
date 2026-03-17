import os
import shutil
import re
import subprocess

print("1. Copying CSS styles from TEP-H0 to TEP-LENS...")
h0_styles = '/Users/matthewsmawfield/www/TEP-H0/site/styles'
lens_styles = '/Users/matthewsmawfield/www/TEP-LENS/site/styles'
os.makedirs(lens_styles, exist_ok=True)
for f in os.listdir(h0_styles):
    if f.endswith('.css'):
        shutil.copy2(os.path.join(h0_styles, f), os.path.join(lens_styles, f))

print("2. Syncing inline HTML styles...")
with open('/Users/matthewsmawfield/www/TEP-H0/site/index.html', 'r') as f:
    h0 = f.read()
with open('/Users/matthewsmawfield/www/TEP-LENS/site/index.html', 'r') as f:
    lens = f.read()

h0_style = re.search(r'<style>(.*?)</style>', h0, re.DOTALL)
if h0_style:
    lens = re.sub(r'<style>.*?</style>', f'<style>{h0_style.group(1)}</style>', lens, flags=re.DOTALL)
    with open('/Users/matthewsmawfield/www/TEP-LENS/site/index.html', 'w') as f:
        f.write(lens)

print("3. Updating figure plotting scripts with standardized colors and sizes...")
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
    
    if 'FIG_SIZE' not in content and 'import matplotlib' in content:
        content = content.replace('set_pub_style, COLORS', 'set_pub_style, COLORS, FIG_SIZE')
        if 'FIG_SIZE' not in content:
            content = content.replace('set_pub_style\n', 'set_pub_style, COLORS, FIG_SIZE\n')
             
    def repl_figsize(match):
        w, h = float(match.group(1)), float(match.group(2))
        if w > 10.5:
            return f'figsize=({w}, 6)'
        return 'figsize=FIG_SIZE'
        
    content = re.sub(r'figsize=\(\s*([0-9\.]+)\s*,\s*([0-9\.]+)\s*\)', repl_figsize, content)
    
    # Replace legacy hardcoded colors with standardized palette keys
    content = content.replace('\"#d62728\"', 'COLORS["red"]')
    content = content.replace('\"#d73027\"', 'COLORS["red"]')
    content = content.replace('\"#4c78a8\"', 'COLORS["primary"]')
    content = content.replace('\"#2166ac\"', 'COLORS["accent"]')
    content = content.replace('\"#fc8d59\"', 'COLORS["success"]')
    content = content.replace('\"darkorange\"', 'COLORS["red"]')
    content = content.replace('\"#d9d9d9\"', 'COLORS["light_gray"]')
    content = content.replace('\"#91bfdb\"', 'COLORS["success"]')
    content = content.replace('\"purple\"', 'COLORS["primary"]')
    content = content.replace('\"#54a24b\"', 'COLORS["info"]')
    
    with open(f_path, 'w') as f:
        f.write(content)

print("4. Regenerating figures...")
os.chdir('/Users/matthewsmawfield/www/TEP-LENS')
subprocess.run(['bash', 'regenerate_all_figures.sh'])

print("5. Copying updated figures to public folder...")
src_dir = '/Users/matthewsmawfield/www/TEP-LENS/results/figures'
dst_dir = '/Users/matthewsmawfield/www/TEP-LENS/site/public/figures'
os.makedirs(dst_dir, exist_ok=True)
for f in os.listdir(src_dir):
    if f.endswith('.png'):
        shutil.copy2(os.path.join(src_dir, f), dst_dir)

print("\nAll tasks completed successfully!")
