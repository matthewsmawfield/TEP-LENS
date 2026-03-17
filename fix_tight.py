import glob
import re

for f in glob.glob('/Users/matthewsmawfield/www/TEP-LENS/scripts/steps/step_*.py'):
    with open(f, 'r') as file:
        content = file.read()
    
    # Restore tight_layout
    content = re.sub(r'#\s*(fig.*?\.tight_layout\(\))', r'\1', content)
    
    with open(f, 'w') as file:
        file.write(content)
