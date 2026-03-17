import json
import os

outputs_dir = '/Users/matthewsmawfield/www/TEP-LENS/results/outputs/'

if os.path.exists(outputs_dir):
    files = os.listdir(outputs_dir)
    print(f"Found {len(files)} files in {outputs_dir}")
    for file in sorted(files):
        if file.endswith('.json'):
            print(file)
else:
    print(f"Could not find {outputs_dir}")
