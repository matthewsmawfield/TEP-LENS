import numpy as np
from scipy import stats

np.random.seed(42)
n_sim = 10000
sigmas = np.array([66., 68., 40., 20., 53., 40., 51., 16.])

p_wilcoxon = []
p_wmean = []

for _ in range(n_sim):
    deltas = np.random.normal(0, sigmas)
    
    # Wilcoxon
    try:
        w_res = stats.wilcoxon(deltas, alternative='greater')
        pw = w_res.pvalue
    except:
        pw = 0.5
    p_wilcoxon.append(pw)
    
    # Weighted mean
    w = 1.0 / sigmas**2
    m = np.sum(w * deltas) / np.sum(w)
    s = 1.0 / np.sqrt(np.sum(w))
    z = m / s
    pm = stats.norm.sf(z)
    p_wmean.append(pm)

p_wilcoxon = np.array(p_wilcoxon)
p_wmean = np.array(p_wmean)

# Calculate Fisher's combined p-value for each simulation
fisher_stat = -2 * (np.log(p_wilcoxon) + np.log(p_wmean))
fisher_p = stats.chi2.sf(fisher_stat, 4)

print(f"False positive rate (alpha=0.05) for Wilcoxon: {np.mean(p_wilcoxon < 0.05):.4f}")
print(f"False positive rate (alpha=0.05) for WMean: {np.mean(p_wmean < 0.05):.4f}")
print(f"False positive rate (alpha=0.05) for Fisher combined: {np.mean(fisher_p < 0.05):.4f}")

# Correlation
print(f"Spearman correlation between p-values: {stats.spearmanr(p_wilcoxon, p_wmean).statistic:.4f}")
