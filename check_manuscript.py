import re
import glob

def check_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    
    issues = []
    
    # Check for 'we'
    we_matches = re.finditer(r'\b(we|our|us)\b', content, re.IGNORECASE)
    for m in we_matches:
        if '<a ' not in content[max(0, m.start()-20):m.start()] and 'href=' not in content[max(0, m.start()-20):m.start()]: # try to ignore links
            issues.append(f"Found '{m.group(0)}' at {m.start()}")
            
    # Check for weasel words
    weasel_matches = re.finditer(r'\b(clearly|obviously|very|undeniably|indisputably)\b', content, re.IGNORECASE)
    for m in weasel_matches:
        issues.append(f"Found weasel word '{m.group(0)}' at {m.start()}")
        
    # Check for future timelines
    future_matches = re.finditer(r'\b(will|plan to|in the future|next year|upcoming)\b', content, re.IGNORECASE)
    for m in future_matches:
        issues.append(f"Found future reference '{m.group(0)}' at {m.start()}")
        
    # Check for bold/strong midway through sentence. 
    # A simple heuristic: <strong> not immediately after a > or newline or punctuation space
    strong_matches = re.finditer(r'([^>\n\.]\s+)(<strong>|<b>)(.*?)(</strong>|</b>)', content)
    for m in strong_matches:
        if m.group(1).strip() != '' and not m.group(1).endswith(':'):
            issues.append(f"Found bold mid-sentence: '{m.group(1)}{m.group(2)}{m.group(3)}{m.group(4)}'")
            
    if issues:
        print(f"\n--- Issues in {filepath} ---")
        for issue in issues[:10]: # print first 10
            print(issue)

for filepath in glob.glob('/Users/matthewsmawfield/www/TEP-LENS/site/components/*.html'):
    check_file(filepath)

