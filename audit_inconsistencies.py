import os
import glob
import re

print("=== 1. Checking for missing files in TEP-LENS/site compared to TEP-H0/site ===")
h0_site = '/Users/matthewsmawfield/www/TEP-H0/site'
lens_site = '/Users/matthewsmawfield/www/TEP-LENS/site'

def get_all_files(base_dir):
    files = set()
    for root, _, filenames in os.walk(base_dir):
        if 'node_modules' in root or '.git' in root or 'dist' in root or '__pycache__' in root:
            continue
        for f in filenames:
            rel_path = os.path.relpath(os.path.join(root, f), base_dir)
            if not f.startswith('.'):
                files.add(rel_path)
    return files

h0_files = get_all_files(h0_site)
lens_files = get_all_files(lens_site)

missing_in_lens = h0_files - lens_files
print("Files missing in TEP-LENS:")
for f in sorted(list(missing_in_lens)):
    print(f"  - {f}")

print("\n=== 2. Checking for mid-sentence bolding in TEP-LENS manuscript components ===")
components_dir = os.path.join(lens_site, 'components')

# Match `**text**` or `<strong>text</strong>` or `<b>text</b>`
# Only flag if it's preceded by lowercase letters or commas/spaces, and followed by lowercase letters or punctuation
# This is a heuristic for "mid-sentence"
bold_pattern = re.compile(r'([a-z,]\s+)(?:\*\*|<strong>|<b>)(.*?)(?:\*\*|</strong>|</b>)(\s+[a-z,])')

found_bold = False
if os.path.exists(components_dir):
    for f in os.listdir(components_dir):
        if f.endswith('.html'):
            with open(os.path.join(components_dir, f), 'r') as file:
                lines = file.readlines()
                for i, line in enumerate(lines):
                    matches = bold_pattern.findall(line)
                    if matches:
                        found_bold = True
                        print(f"File {f}, Line {i+1}:")
                        for match in matches:
                            print(f"  Found bold: '{match[0]}**{match[1]}**{match[2]}'")

if not found_bold:
    print("No obvious mid-sentence bolding found.")

print("\n=== 3. Checking for missing metadata/links in index.html ===")
with open(os.path.join(lens_site, 'index.html'), 'r') as f:
    lens_index = f.read()

checks = {
    "Journal alternate link": r'<link rel="alternate" href="journal\.html">',
    "Manifest journal": r'<link rel="manifest" href="manifest-journal\.json">',
    "Citation JSON": r'citation\.json'
}

for name, pattern in checks.items():
    if re.search(pattern, lens_index):
        print(f"{name}: Found")
    else:
        print(f"{name}: MISSING")

