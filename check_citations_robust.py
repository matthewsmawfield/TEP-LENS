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

print("References:", ref_authors_years)

# Check all papers cited
cites_raw = []
for filepath in glob.glob('/Users/matthewsmawfield/www/TEP-LENS/site/components/*.html'):
    if '7_references.html' in filepath: continue
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Catch formats like "Author et al. 2023" or "Author 2023" or "Author & Author 2023"
    matches = re.findall(r'([A-Z][a-z]+(?: et al\.)?(?:,)? (?:20\d\d|19\d\d)[a-z]?)', content)
    cites_raw.extend(matches)

# print("Raw cites:", set(cites_raw))
