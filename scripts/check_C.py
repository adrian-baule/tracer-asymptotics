"""
Validate the non-universal constant C for the bmshort process via <A^2>(T) ~ C * log(T).

For a 2D BM swimmer with a short-range Gaussian force, the scattering increment variance
grows logarithmically: <A^2>(T) ~ C * log(T). This follows from 2D BM being null-recurrent:
the occupation time of the force region scales as (sigma^2/D) * log(DT/b^2).

Usage:
    python scripts/check_C.py --data data/mu2_bmshort.csv [--fit_frac 0.5]

Loads mu2_bmshort.csv (columns: t, mu2), plots mu2(t) vs log(t), fits a straight line
C * log(t) to the large-t regime, and reports C.
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default="data/var_bmshort.csv")
    parser.add_argument("--fit_frac", type=float, default=0.5,
                        help="Fraction of large-t data to use for log fit (default 0.5)")
    args = parser.parse_args()

    data = np.loadtxt(args.data, delimiter=",", skiprows=1)
    t, mu2 = data[:, 0], data[:, 1]

    # Fit mu2 = C * log(t) to large-t regime (no intercept, through origin in log space)
    t_thresh = t[-1] * (1.0 - args.fit_frac)
    mask = t >= t_thresh
    log_t = np.log(t[mask])
    # Weighted least-squares: mu2 = C * log(t), no intercept
    C = np.sum(log_t * mu2[mask]) / np.sum(log_t ** 2)

    print(f"Log fit over t >= {t_thresh:.1f}:")
    print(f"  <A^2>(t) = C * log(t)")
    print(f"  C        = {C:.6g}")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left: mu2 vs t (linear scale)
    ax = axes[0]
    ax.plot(t, mu2, 'k.', ms=2, label=r'$\langle A^2 \rangle(t)$')
    ax.plot(t, C * np.log(t), 'r--', lw=1.5, label=f'Log fit: C·log(t),  C = {C:.4g}')
    ax.axvline(t_thresh, color='gray', lw=0.8, linestyle=':', label='fit region start')
    ax.set_xlabel('t')
    ax.set_ylabel(r'$\langle A^2 \rangle(t)$')
    ax.set_title('bmshort: weighted second moment vs time')
    ax.legend(fontsize=8)

    # Right: mu2 vs log(t) — should be linear if scaling is correct
    ax = axes[1]
    ax.plot(np.log(t), mu2, 'k.', ms=2, label=r'$\langle A^2 \rangle(t)$')
    ax.plot(np.log(t), C * np.log(t), 'r--', lw=1.5, label=f'C·log(t),  C = {C:.4g}')
    ax.axvline(np.log(t_thresh), color='gray', lw=0.8, linestyle=':', label='fit region start')
    ax.set_xlabel('log(t)')
    ax.set_ylabel(r'$\langle A^2 \rangle(t)$')
    ax.set_title(r'bmshort: $\langle A^2 \rangle$ vs log(t) — expect linear')
    ax.legend(fontsize=8)

    plt.suptitle(r'bmshort: $\langle A^2 \rangle(T) \sim C\,\log T$ (2D BM, null-recurrent)', y=1.01)
    plt.tight_layout()
    out = args.data.replace('.csv', '_fit.png').replace('mu2_bmshort', 'check_C')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {out}")
    plt.show()

if __name__ == "__main__":
    main()
