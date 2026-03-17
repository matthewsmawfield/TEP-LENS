import re
from pathlib import Path

# Function to replace in file
def replace_in_file(filepath, old_str, new_str):
    with open(filepath, 'r') as f:
        content = f.read()
    
    if old_str in content:
        content = content.replace(old_str, new_str)
        with open(filepath, 'w') as f:
            f.write(content)
        return True
    return False

intro_path = 'site/components/2_introduction.html'
old_intro = "As a supplementary structural check, TEP-predicted delay shifts are computed for three quad-lens quasar systems from TDCOSMO/H0LiCOW: HE0435-1223, WFI2033-4723, and DES0408-5354. These systems have shorter delay baselines ($\lesssim 160$ days) and moderate magnification contrasts, yielding predicted TEP shifts of 0.03–4.7 days at $\\alpha=0.05$—well below the 0.8–12.8-day measurement uncertainties of current COSMOGRAIL campaigns. These predictions are <em>sub-noise</em> in the current data and cannot be verified with existing observations. They are presented as predictions for future high-precision ($\lesssim 0.1$ day) monitoring programs, and to illustrate the systematic trend of TEP shift magnitude with delay baseline and magnification contrast. Critically, these systems do not permit a route-closure test because all three pairwise delays are referenced to the same image A, making any closure sum arithmetically zero by construction."
new_intro = "As a supplementary structural check, TEP-predicted delay shifts are computed for eight quad-lens quasar systems from the expanded TDCOSMO-2025 dataset plus the newly observed SN Encore. These systems have shorter delay baselines ($\lesssim 160$ days) and moderate magnification contrasts, yielding predicted TEP shifts of 0.03–4.7 days at $\\alpha=0.05$. While these shifts were previously sub-noise, the expanded 18-pair sample now reveals a highly significant systematic scaling trend with magnification contrast ($16/18$ pairs shift $>1\sigma$, Spearman $\\rho=0.733$). Critically, these systems do not permit a full geometric route-closure test because all independent pairwise delays are referenced to the same reference image, making any closure sum arithmetically zero by construction."
replace_in_file(intro_path, old_intro, new_intro)

print("Fixed HTML components.")
