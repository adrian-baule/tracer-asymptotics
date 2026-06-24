"""
Validate the non-universal constant C for the bmshort process via <A^2>(T) = 2*C*T.

Usage:
    python scripts/check_C.py --data data/mu2_bmshort.csv [--fit_frac 0.5]

Loads mu2_bmshort.csv (columns: t, mu2), plots mu2(t) vs t, fits a straight line
to the large-t regime to estimate slope = 2*C, and reports C.
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default="data/mu2_bmshort.csv")
    parser.add_argument("--fit_frac", type=float, default=0.5,
                        help="Fraction of large-t data to use for linear fit (default 0.5)")
    args = parser.parse_args()

    data = np.loadtxt(args.data, delimiter=",", skiprows=1)
    t, mu2 = data[:, 0], data[:, 1]

    # Fit straight line through origin to large-t regime
    t_thresh = t[-1] * (1.0 - args.fit_frac)
    mask = t >= t_thresh
    # Weighted least-squares: mu2 = slope * t, no intercept
    slope = np.sum(t[mask] * mu2[mask]) / np.sum(t[mask] ** 2)
    C = slope / 2.0

    print(f"Linear fit over t >= {t_thresh:.1f}:")
    print(f"  slope = 2*C = {slope:.6g}")
    print(f"  C           = {C:.6g}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(t, mu2, 'k.', ms=2, label=r'$\langle A^2 \rangle(t)$')
    ax.plot(t, slope * t, 'r--', lw=1.5, label=f'Linear fit: slope = {slope:.4g}  (C = {C:.4g})')
    ax.axvline(t_thresh, color='gray', lw=0.8, linestyle=':')
    ax.set_xlabel('t')
    ax.set_ylabel(r'$\langle A^2 \rangle(t)$')
    ax.set_title('bmshort: weighted second moment vs time\n'
                 r'Expect $\langle A^2 \rangle = 2Ct$ (linear growth)')
    ax.legend()
    plt.tight_layout()
    out = args.data.replace('.csv', '_fit.png').replace('mu2_bmshort', 'check_C')
    plt.savefig(out, dpi=150)
    print(f"Plot saved to {out}")
    plt.show()

if __name__ == "__main__":
    main()
