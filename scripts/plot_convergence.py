"""
Check convergence of the distribution with simulation time T.
Expects data files: data/output_T100.csv, data/output_T1000.csv, data/output_T10000.csv
Run the simulation three times with --T 100, 1000, 10000 and --output data/output_T{T}.csv
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os

def weighted_hist(A, w, bins):
    w = w / w.mean()
    counts, edges = np.histogram(A, bins=bins, weights=w, density=True)
    return 0.5*(edges[:-1]+edges[1:]), counts

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", default="data")
    args = parser.parse_args()

    T_values = [100, 1000, 10000]
    fig, ax = plt.subplots(figsize=(7,5))
    bins = np.linspace(-20, 20, 80)

    for T in T_values:
        fname = os.path.join(args.data_dir, f"output_T{T}.csv")
        if not os.path.exists(fname):
            print(f"Missing: {fname}, skipping")
            continue
        data = np.loadtxt(fname, delimiter=",", skiprows=1)
        b, A, w = data[:,0], data[:,1], data[:,2]
        centres, counts = weighted_hist(A, w, bins)
        ax.semilogy(centres, counts, label=f'T={T}')

    ax.set_xlabel('A')
    ax.set_ylabel('P(A)')
    ax.set_title('Convergence with simulation time T')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(args.data_dir, 'convergence.png'), dpi=150)
    print(f"Saved {args.data_dir}/convergence.png")
    plt.show()

if __name__ == "__main__":
    main()
