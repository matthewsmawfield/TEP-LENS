with open('/Users/matthewsmawfield/www/TEP-LENS/site/components/4_results.html', 'r') as f:
    results_content = f.read()

import re

# Update caption of figure 5
results_content = re.sub(r'The weighted mean \$\\bar\{\\alpha\} = -0\.055 \\pm 0\.044\$ \(purple\)', 
                         r'The weighted mean $\\bar{\\alpha} = -0.057 \\pm 0.060$ (purple)', results_content)

with open('/Users/matthewsmawfield/www/TEP-LENS/site/components/4_results.html', 'w') as f:
    f.write(results_content)

print("Updated 4_results.html")
