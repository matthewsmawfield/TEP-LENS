import re

files = [
    "/Users/matthewsmawfield/www/TEP-LENS/site/components/4_results.html",
    "/Users/matthewsmawfield/www/TEP-LENS/site/components/5_discussion.html"
]

all_content = ""
for f in files:
    with open(f, "r") as file:
        all_content += file.read()

# Find all figure numbers
matches = re.findall(r'<figcaption><strong>Figure (\d+):</strong>', all_content)
print(f"Current figure numbers in order: {matches}")

counter = 1
def replacer(match):
    global counter
    res = f'<figcaption><strong>Figure {counter}:</strong>'
    counter += 1
    return res

for filename in files:
    with open(filename, "r") as file:
        content = file.read()
    
    new_content = re.sub(r'<figcaption><strong>Figure \d+:</strong>', replacer, content)
    
    with open(filename, "w") as file:
        file.write(new_content)

print(f"Renumbered up to Figure {counter-1}")
