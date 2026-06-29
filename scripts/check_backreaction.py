"""
Analyse the back-reaction diagnostic for the bmlong process.

Loads var_bmlong.csv (columns: t, VarA1, VarA2, meanA1, meanA2, ratio) and produces:
  1. Log-log plot of Var(A1) and Var(A2) vs t
     -- expect Var(A1) ~ t*log(t) (marginal Coulomb), Var(A2) ~ t*(log t)^a
  2. Ratio Var(A2)/Var(A1) vs log(t)
     -- flat => a=1 (adiabatic truncation benign); linear slope => a-1 (log promotion)
  3. Mean(A1) and Mean(A2) vs t (symmetry check; should be ~0)

Usage:
    python scripts/check_backreaction.py --data data/var_bmlong.csv [--fit_frac 0.5]
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default="data/var_bmlong.csv")
    parser.add_argument("--fit_frac", type=float, default=0.5,
                        help="Fraction of large-t data for ratio slope fit (default 0.5)")
    args = parser.parse_args()

    data = np.loadtxt(args.data, delimiter=",", skiprows=1)
    t      = data[:, 0]
    VarA1  = data[:, 1]
    VarA2  = data[:, 2]
    meanA1 = data[:, 3]
    meanA2 = data[:, 4]
    ratio  = data[:, 5]
    kurtA1 = data[:, 6] if data.shape[1] > 6 else None

    # Fit slope of ratio vs log(t) over large-t regime
    t_thresh = t[-1] * (1.0 - args.fit_frac)
    mask = (t >= t_thresh) & (VarA1 > 0)
    log_t = np.log(t[mask])
    r_fit = ratio[mask]
    # Linear fit: ratio = a0 + slope * log(t)
    coeffs = np.polyfit(log_t, r_fit, 1)
    slope_ratio, intercept_ratio = coeffs
    ratio_final = ratio[mask][-1]

    print(f"Fit region: t >= {t_thresh:.1f}")
    print(f"Var(A2)/Var(A1) at t=T:      {ratio_final:.4g}")
    print(f"Slope of ratio vs log(t):     {slope_ratio:.4g}")
    print(f"  => a - 1 = {slope_ratio:.4g}  (flat => a=1; positive => log promotion)")
    print(f"  Note: A2 is O(mu^2) and noisier; interpret large-t trend only.")
    print()

    # Symmetry check
    max_mA1 = np.max(np.abs(meanA1))
    max_mA2 = np.max(np.abs(meanA2))
    flag = " [WARNING: not ~0]" if max_mA1 > 0.1 * np.sqrt(VarA1[-1]) else ""
    print(f"max|mean(A1)|: {max_mA1:.3g}{flag}")
    flag2 = " [WARNING: not ~0]" if max_mA2 > 0.1 * np.sqrt(VarA2[-1]) else ""
    print(f"max|mean(A2)|: {max_mA2:.3g}{flag2}")

    fig, axes = plt.subplots(1, 4 if kurtA1 is not None else 3, figsize=(19 if kurtA1 is not None else 15, 5))

    # Panel 1: Var(A1) and Var(A2) vs t on log-log
    ax = axes[0]
    ax.loglog(t, VarA1, 'b.-', ms=3, label=r'Var($A_1$)')
    ax.loglog(t, VarA2, 'r.-', ms=3, label=r'Var($A_2$)')
    # Reference: t*log(t)
    ref = t * np.log(np.maximum(t, 1.1))
    ref *= VarA1[len(t)//2] / ref[len(t)//2]
    ax.loglog(t, ref, 'b--', lw=1, alpha=0.5, label=r'$\propto t\log t$')
    ax.set_xlabel('t')
    ax.set_ylabel('Variance')
    ax.set_title(r'Var($A_1$) and Var($A_2$) vs $t$')
    ax.legend(fontsize=8)

    # Panel 2: ratio vs log(t) with linear fit
    ax = axes[1]
    ax.plot(np.log(t), ratio, 'k.-', ms=3, label='Var(A2)/Var(A1)')
    log_t_all = np.log(t[t > 0])
    ax.plot(log_t_all, intercept_ratio + slope_ratio * log_t_all, 'r--', lw=1.5,
            label=f'Fit: slope = {slope_ratio:.3g}')
    ax.axvline(np.log(t_thresh), color='gray', lw=0.8, linestyle=':', label='fit region')
    ax.set_xlabel(r'$\log t$')
    ax.set_ylabel(r'Var($A_2$)/Var($A_1$)')
    ax.set_title(r'Back-reaction ratio vs $\log t$')
    ax.legend(fontsize=8)

    # Panel 3: means vs t (symmetry check)
    ax = axes[2]
    ax.plot(t, meanA1, 'b.-', ms=2, label=r'mean($A_1$)')
    ax.plot(t, meanA2, 'r.-', ms=2, label=r'mean($A_2$)')
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlabel('t')
    ax.set_ylabel('Mean')
    ax.set_title(r'Symmetry check: $\langle A_1\rangle$ and $\langle A_2\rangle$ (should be $\approx 0$)')
    ax.legend(fontsize=8)

    # Panel 4: kurtosis of A1 vs t (if available)
    if kurtA1 is not None:
        ax = axes[3]
        ax.semilogx(t, kurtA1, 'g.-', ms=3, label=r'Kurt($A_1$)')
        ax.axhline(3, color='k', lw=0.8, linestyle='--', label='Gaussian (=3)')
        ax.set_xlabel('t')
        ax.set_ylabel(r'Kurtosis($A_1$)')
        ax.set_title(r'Kurtosis of $A_1$ vs $t$')
        ax.legend(fontsize=8)

    plt.suptitle('bmlong back-reaction diagnostic', y=1.01)
    plt.tight_layout()
    out = args.data.replace('var_bmlong.csv', 'check_backreaction.png').replace('.csv', '_backreaction.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {out}")
    plt.show()

if __name__ == "__main__":
    main()
