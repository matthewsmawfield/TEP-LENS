with open('/Users/matthewsmawfield/www/TEP-LENS/zenodo.txt', 'r') as f:
    zenodo_content = f.read()

import re

zenodo_content = zenodo_content.replace('0.057 \pm 0.060', '-0.057 \pm 0.060')

with open('/Users/matthewsmawfield/www/TEP-LENS/zenodo.txt', 'w') as f:
    f.write(zenodo_content)

with open('/Users/matthewsmawfield/www/TEP-LENS/README.md', 'r') as f:
    readme_content = f.read()

readme_content = readme_content.replace('0.057 \pm 0.060', '-0.057 \pm 0.060')

with open('/Users/matthewsmawfield/www/TEP-LENS/README.md', 'w') as f:
    f.write(readme_content)

with open('/Users/matthewsmawfield/www/TEP-LENS/site/components/1_abstract.html', 'r') as f:
    abstract_content = f.read()

abstract_content = abstract_content.replace('0.057 \pm 0.060', '-0.057 \pm 0.060')

with open('/Users/matthewsmawfield/www/TEP-LENS/site/components/1_abstract.html', 'w') as f:
    f.write(abstract_content)

print("Fixed negative sign")
