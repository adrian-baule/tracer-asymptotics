"""
Plot weighted histogram of scattering increments A vs theory.
Usage: python scripts/plot_distribution.py --data data/output.csv
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as gammaf

def mittag_leffler_moments(n, nu=0.25):
    """Moments of bilateral Mittag-Leffler M_nu: <z^n> = (1+(-1)^n)/2 * n! / Gamma(nu*n+1)"""
    if n % 2 != 0:
        return 0.0
    return float(np.math.factorial(n)) / gammaf(nu * n + 1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default="data/output.csv")
    parser.add_argument("--nbins", type=int, default=100)
    args = parser.parse_args()

    data = np.loadtxt(args.data, delimiter=",", skiprows=1)
    b, A, w = data[:,0], data[:,1], data[:,2]

    # Normalise weights
    w = w / w.mean()

    # Weighted histogram
    bins = np.linspace(np.percentile(A, 1), np.percentile(A, 99), args.nbins)
    counts, edges = np.histogram(A, bins=bins, weights=w, density=True)
    centres = 0.5 * (edges[:-1] + edges[1:])

    # Diagnostics
    N_eff = w.sum()**2 / (w**2).sum()
    mean_A = np.average(A, weights=w)
    var_A  = np.average((A - mean_A)**2, weights=w)
    print(f"N_eff / N_traj = {N_eff / len(w):.3f}")
    print(f"Weighted mean(A) = {mean_A:.4f} (should be ~0)")
    print(f"Weighted var(A)  = {var_A:.4f}")
    print(f"Weighted std(A)  = {np.sqrt(var_A):.4f}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(centres, counts, 'k.', ms=3, label='Simulation')

    # Rescale A by std to compare shape with M_{1/4}
    sigma = np.sqrt(var_A)
    z = centres / sigma
    # M_{1/4} moments: <z^2> = 2/Gamma(3/4)
    m2 = 2.0 / gammaf(0.75)
    z_theory = np.linspace(z.min(), z.max(), 500)
    # M_{1/4} via its characteristic function (numerical inverse Fourier)
    # For now just mark the Gaussian with same variance for comparison
    from scipy.stats import norm
    ax.semilogy(centres, norm.pdf(centres, 0, sigma), 'r--', label='Gaussian (same var)', alpha=0.7)

    ax.set_xlabel('A')
    ax.set_ylabel('P(A)')
    ax.set_title('Single-scattering displacement distribution')
    ax.legend()
    plt.tight_layout()
    plt.savefig('data/distribution.png', dpi=150)
    print("Saved data/distribution.png")
    plt.show()

if __name__ == "__main__":
    main()
