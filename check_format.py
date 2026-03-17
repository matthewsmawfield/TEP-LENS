import os
import re

def extract_elements(directory):
    elements = {'table': [], 'figure': [], 'callout': [], 'blockquote': []}
    if not os.path.exists(directory):
        return elements
    for fname in os.listdir(directory):
        if not fname.endswith('.html'): continue
        path = os.path.join(directory, fname)
        with open(path, 'r') as f:
            content = f.read()
        
        # Tables
        for m in re.finditer(r'<table[^>]*>', content):
            elements['table'].append((fname, m.group(0)))
            
        # Figures
        for m in re.finditer(r'<figure[^>]*>', content):
            elements['figure'].append((fname, m.group(0)))
            
        # Callouts
        for m in re.finditer(r'<div[^>]*class="[^"]*callout[^"]*"[^>]*>', content):
            elements['callout'].append((fname, m.group(0)))
        for m in re.finditer(r'<p[^>]*class="[^"]*callout[^"]*"[^>]*>', content):
            elements['callout'].append((fname, m.group(0)))
            
        # Blockquotes
        for m in re.finditer(r'<blockquote[^>]*>', content):
            elements['blockquote'].append((fname, m.group(0)))
            
    return elements

h0 = extract_elements('/Users/matthewsmawfield/www/TEP-H0/site/components')
lens = extract_elements('/Users/matthewsmawfield/www/TEP-LENS/site/components')

print("=== TEP-H0 Tables ===")
for e in h0['table']: print(f"  {e[0]}: {e[1]}")
print("\n=== TEP-LENS Tables ===")
for e in lens['table']: print(f"  {e[0]}: {e[1]}")

print("\n=== TEP-H0 Callouts ===")
for e in h0['callout']: print(f"  {e[0]}: {e[1]}")
print("\n=== TEP-LENS Callouts ===")
for e in lens['callout']: print(f"  {e[0]}: {e[1]}")

print("\n=== TEP-H0 Figures ===")
for e in h0['figure']: print(f"  {e[0]}: {e[1]}")
print("\n=== TEP-LENS Figures ===")
for e in lens['figure']: print(f"  {e[0]}: {e[1]}")
