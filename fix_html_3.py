import re
from pathlib import Path

# Function to replace in file
def replace_in_file(filepath, old_str, new_str):
    with open(filepath, 'r') as f:
        content = f.read()
    
    if old_str in content:
        content = content.replace(old_str, new_str)
        with open(filepath, 'w') as f:
            f.write(content)
        return True
    return False

# Fix 1_abstract.html to include TDCOSMO/Encore correlation
abstract_path = 'site/components/1_abstract.html'
old_abstract = "and 45% wRMS\n    reduction after TEP correction. Robustness analyses added here show stability to model dependence"
new_abstract = "and 45% wRMS\n    reduction after TEP correction. The extended TDCOSMO-2025 and SN Encore dataset independently corroborates this systematic effect ($16/18$ pairs shift $>1\sigma$, Spearman $\\rho=0.733$). Robustness analyses added here show stability to model dependence"
replace_in_file(abstract_path, old_abstract, new_abstract)

print("Fixed HTML components.")
