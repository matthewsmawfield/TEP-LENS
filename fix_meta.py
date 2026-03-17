with open('/Users/matthewsmawfield/www/TEP-LENS/site/components/1_abstract.html', 'r') as f:
    abstract_content = f.read()

import re
# Extract abstract text
abstract_text_match = re.search(r'<p><strong>Abstract\.</strong>\s*(.*?)\s*</p>\s*<p>\s*(.*?)\s*</p>', abstract_content, re.DOTALL)

abstract_part1 = abstract_text_match.group(1).replace('\n    ', ' ').replace('\n', ' ')
abstract_part2 = abstract_text_match.group(2).replace('\n    ', ' ').replace('\n', ' ')

full_abstract = abstract_part1 + " " + abstract_part2
full_abstract = re.sub(r'\s+', ' ', full_abstract)

# Also get keywords
keywords_match = re.search(r'<p><strong>Keywords:</strong>(.*?)</p>', abstract_content)
keywords = keywords_match.group(1).strip()

zenodo_text = f"""{full_abstract}

Keywords:{keywords}

Update note (v14 robustness extension): Steps 11-14 add model-dependence stress tests, microlensing nuisance Monte Carlo, hierarchical Bayesian GR-vs-TEP comparison with H0pe-informed prior sensitivity, and public H0LiCOW/TDCOSMO chain ingestion for immediate external-data integration. The directional evidence pattern remains robust; formal model-selection remains non-decisive at current lens-model uncertainty.

Open Science Statement:
This work is a preprint and is open to community review, ideas, and collaboration. All analysis code, data, and manuscripts are open source and available at https://github.com/matthewsmawfield/TEP-LENS. Feedback and contributions to further test these results are welcome."""

with open('/Users/matthewsmawfield/www/TEP-LENS/zenodo.txt', 'w') as f:
    f.write(zenodo_text)

print("Updated zenodo.txt")

with open('/Users/matthewsmawfield/www/TEP-LENS/README.md', 'r') as f:
    readme_content = f.read()

# Instead of using re.sub with replacement string that might have backslashes, use string replace or simple slicing
match = re.search(r'## Overview\n(.*?)\n\n## Repository Structure', readme_content, flags=re.DOTALL)
if match:
    readme_new = readme_content[:match.start(1)] + full_abstract + readme_content[match.end(1):]
    with open('/Users/matthewsmawfield/www/TEP-LENS/README.md', 'w') as f:
        f.write(readme_new)
    print("Updated README.md")

