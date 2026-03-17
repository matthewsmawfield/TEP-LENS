import re
import glob

# Alpha inferred in abstract is "0.057 \pm 0.060" according to zenodo.txt, let's fix the abstract and results.
# Wait, zenodo text actually showed:
# (4) per-model coupling inference - weighted mean $\bar{\alpha}_{\rm inferred} = 0.057 \pm 0.060$

# Let's see what's in the abstract.html
with open('/Users/matthewsmawfield/www/TEP-LENS/site/components/1_abstract.html', 'r') as f:
    abstract_content = f.read()
    
# Let's check the alpha value in results
with open('/Users/matthewsmawfield/www/TEP-LENS/site/components/4_results.html', 'r') as f:
    results_content = f.read()

# The abstract does not explicitly state the alpha inferred value, but results does.
print("Abstract alpha check:", re.search(r'alpha_{\\rm inferred}', abstract_content))

# Look in results for the 0.055 vs 0.057
for m in re.finditer(r'alpha_{\\rm inferred}(.*?)(0\.05[57])', results_content):
    print(f"Results matches: {m.group(0)}")

