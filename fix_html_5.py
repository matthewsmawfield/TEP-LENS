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

disc_path = 'site/components/5_discussion.html'

# Update the Evidence Synthesis Table
old_row = """        <tr>
            <td><strong>Delay–$\mu$ correlation</strong><br><small>Pearson $r=0.93$, $n=5$</small></td>
            <td>Correlation test</td>
            <td>Positive slope 172.6±38.8 d per unit $1/\mu$; SX dominates — inner-cross ordering not reproduced. <em>Not independent of the sign-based strand: SX is the single driving data point in both.</em></td>
            <td>$p=0.011$ (2.3$\sigma$, one-sided)</td>
            <td>✓ Observed (SX-driven; correlated with strand 1)</td>
        </tr>"""

new_row = """        <tr>
            <td><strong>TDCOSMO+Encore Shear</strong><br><small>Spearman $\\rho=0.733$, $n=18$ pairs</small></td>
            <td>Correlation test</td>
            <td>Positive monotonic scaling of TEP delay shift with relative magnification across 9 independent systems. 16/18 pairs shift $>1\sigma$.</td>
            <td>$p=0.0005$ (highly significant)</td>
            <td>✓ Observed (Independent dataset)</td>
        </tr>
        <tr>
            <td><strong>Delay–$\mu$ correlation</strong><br><small>Pearson $r=0.93$, $n=5$</small></td>
            <td>Correlation test</td>
            <td>Positive slope 172.6±38.8 d per unit $1/\mu$; SX dominates — inner-cross ordering not reproduced. <em>Not independent of the sign-based strand: SX is the single driving data point in both.</em></td>
            <td>$p=0.011$ (2.3$\sigma$, one-sided)</td>
            <td>✓ Observed (SX-driven; correlated with strand 1)</td>
        </tr>"""

replace_in_file(disc_path, old_row, new_row)

print("Fixed Evidence Synthesis Table.")
