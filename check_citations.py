import re
import glob

# Collect all references
with open('/Users/matthewsmawfield/www/TEP-LENS/site/components/7_references.html', 'r') as f:
    ref_content = f.read()

ref_items = re.findall(r'<p class="ref-item">(.*?)</p>', ref_content)
ref_authors_years = []
for ref in ref_items:
    match = re.search(r'^([^,]+).*?(\d{4})', ref)
    if match:
        ref_authors_years.append((match.group(1).strip(), match.group(2).strip()))
    else:
        print("Could not parse ref:", ref)

# Find all citations in text
citations = set()
for filepath in glob.glob('/Users/matthewsmawfield/www/TEP-LENS/site/components/*.html'):
    if '7_references.html' in filepath: continue
    with open(filepath, 'r') as f:
        content = f.read()
        
    # Pattern for (Author et al. 2023) or Author et al. (2023)
    cites1 = re.findall(r'\(([A-Za-z]+)\s+(?:et al\.\s*)?(\d{4})[a-z]?(?:,.*?)?\)', content)
    cites2 = re.findall(r'([A-Za-z]+)\s+(?:et al\.\s*)?\((\d{4})[a-z]?(?:,.*?)?\)', content)
    
    for author, year in cites1 + cites2:
        if author not in ['e', 'g', 'i', 'SNR', 'Fig', 'Table', 'Paper']: # Exclude some false positives
            citations.add((author, year))

print("Citations found in text:")
for c in sorted(citations):
    found = False
    for r_author, r_year in ref_authors_years:
        if c[0] in r_author and c[1] == r_year:
            found = True
            break
        # Sometimes there's a small mismatch like 'Planck' vs 'Planck Collaboration'
        if 'Planck' in c[0] and 'Planck' in r_author:
            found = True
            break
        if 'Gaia' in c[0] and 'Gaia' in r_author:
            found = True
            break
    
    status = "OK" if found else "MISSING IN REFERENCES"
    print(f"{c[0]} {c[1]} - {status}")

