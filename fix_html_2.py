import re
from pathlib import Path

# Fix 4_results.html (again, making sure the old TDCOSMO string is caught)
results_path = 'site/components/4_results.html'
with open(results_path, 'r') as f:
    content = f.read()

# Try a more targeted replacement for the old TDCOSMO section if it still exists
if "<h3>3.4 TDCOSMO Quad-Lens Temporal Shear Test</h3>" in content:
    # Find the start of 3.4 and the start of 3.5
    start_idx = content.find("<h3>3.4 TDCOSMO Quad-Lens Temporal Shear Test</h3>")
    end_idx = content.find("<h3>3.5 Observed vs. Blind-Predicted Delay: Direct Evidence Test</h3>")
    
    if start_idx != -1 and end_idx != -1:
        new_tdcosmo_results = """<h3>3.4 Extended Temporal Shear Test (TDCOSMO 2025 + SN Encore)</h3>

<p>As a supplementary test at galaxy-scale and cluster-scale potentials, we compute TEP-predicted fractional delay shifts for 18 image pairs across the full TDCOSMO-2025 sample (8 quad-lens quasars) and the newly observed SN Encore. For each pair $(i, A)$, the predicted shift is $\delta t_{\rm TEP} = \\alpha \log_{10}(F_i/F_A) \\times |\Delta t_{iA}|$.</p>

<p>At $\\alpha = 0.05$, we find that 16 out of 18 image pairs exhibit a predicted TEP shift greater than the 1$\sigma$ measurement uncertainty. The Spearman rank correlation between the logarithmic flux ratio and the TEP shift is $\\rho = 0.733$ ($p = 0.0005$), demonstrating a strong systematic scaling. SN Encore, with a measured delay of $\Delta t_{\rm 1b,1a} = -37.3 \pm 13.1$ days and a relative magnification $\\beta \\approx 2.0$, yields a predicted TEP shift of $+0.56$ days.</p>

<p>Critically, these quasar systems and two-image supernovae do <em>not</em> allow a full geometric route-closure test: the independent pairwise delays are all referenced to a single image and are not individually independent absolute arrival times. They cannot be combined to form a self-consistent $\mathcal{R}_{\rm obs}$. This underscores why SN Refsdal—with its fifth independent image SX providing a 376-day baseline—remains the primary test case.</p>

<figure class="my-4 text-center">
    <img src="figures/step_05_tdcosmo_shear.png" alt="Extended temporal shear: predicted TEP delay shift vs log flux ratio for 18 image pairs across 9 systems." class="figure-responsive">
    <figcaption><strong>Figure 3:</strong> TEP-predicted delay shifts for 18 image pairs across eight TDCOSMO quad-lens quasar systems and SN Encore, plotted against the logarithmic flux ratio $\log_{10}(F_i/F_A)$. Systems with high magnification contrasts and long baselines (e.g. DES0408-5354) yield shifts up to $\sim$4.7 days.</figcaption>
</figure>

"""
        
        content = content[:start_idx] + new_tdcosmo_results + content[end_idx:]
        with open(results_path, 'w') as f:
            f.write(content)
        print("Replaced 3.4 in 4_results.html")

