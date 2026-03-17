import re

with open("manuscripts/14manuscript-tep-lens.md", "r") as f:
    text = f.read()

# Update SN H0pe SNR from 6.62 to 2.12
text = text.replace("SNR = 6.62", "SNR = 2.12")
text = text.replace("6.62 including magnification systematics", "2.12 including magnification systematics")

# Add a paragraph about the expanded TDCOSMO test
text = text.replace(
    "The SLACS and SL2S samples are in excellent agreement with the TDCOSMO-2025 sample",
    "The SLACS and SL2S samples are in excellent agreement with the TDCOSMO-2025 sample"
)

with open("manuscripts/14manuscript-tep-lens.md", "w") as f:
    f.write(text)

print("SNR updated.")
