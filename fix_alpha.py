with open('/Users/matthewsmawfield/www/TEP-LENS/site/components/4_results.html', 'r') as f:
    results_content = f.read()

import re

# Update the section 3.6.2
results_content = re.sub(r'\\bar\{\\alpha\}_\{\\rm inferred\} = -0\.055 \\pm 0\.044\\quad \(z = 1\.26 \\text\{ from zero\}\)', 
                         r'\\bar{\\alpha}_{\\rm inferred} = -0.057 \\pm 0.060\\quad (z = 0.95 \\text{ from zero})', results_content)

with open('/Users/matthewsmawfield/www/TEP-LENS/site/components/4_results.html', 'w') as f:
    f.write(results_content)

print("Updated 4_results.html")
