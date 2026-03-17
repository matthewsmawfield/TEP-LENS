import re

with open("manuscripts/14manuscript-tep-lens.md", "r") as f:
    text = f.read()

# Update the H0pe sensitivity caption in Figure 14
text = text.replace("SN H0pe future sensitivity ($z = 6.6$)", "SN H0pe future sensitivity ($z = 2.1$)")

with open("manuscripts/14manuscript-tep-lens.md", "w") as f:
    f.write(text)

print("Figure caption updated.")
