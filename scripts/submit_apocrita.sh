#!/bin/bash
# Submit one SLURM job per swimmer process on Apocrita.
# Usage: bash scripts/submit_apocrita.sh
# Assumes binaries are already built in ./build/

mkdir -p logs

for PROC in aoup abp rtp levy1 levy2 bmshort bmlong bmshort_exact rtp_exact; do
    # Defaults, overridden per process below
    N_TRAJ=1000000
    T_SIM=1000
    WALLTIME=01:00:00
    EXTRA_ARGS=""

    if [ "$PROC" = "levy2" ]; then
        EXTRA_ARGS="--beta 0.5 --tau_0 1.0"
    elif [ "$PROC" = "bmshort" ]; then
        EXTRA_ARGS="--D_bm 1.0 --V_0 1.0 --sigma 1.0"
    elif [ "$PROC" = "bmlong" ]; then
        # b_max auto-computed as sqrt(4*D_bm*T); do not hardcode
        EXTRA_ARGS="--sigma 1.0 --b_min 0.5 --mu 1.0"
    elif [ "$PROC" = "bmshort_exact" ]; then
        EXTRA_ARGS="--D_bm 1.0 --V_0 1.0 --sigma 1.0"
    elif [ "$PROC" = "rtp_exact" ]; then
        # T=1e4 at dt=1e-3 is 1e7 steps per trajectory (10x the other runs),
        # so N_traj is cut and the walltime raised to match.
        # b_min=1.0 and b_max (auto-computed from D_eff) are left at defaults.
        N_TRAJ=100000
        T_SIM=10000
        WALLTIME=08:00:00
        EXTRA_ARGS="--v_A 1.0 --omega 1.0 --p 1.0"
    fi

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=tracer_${PROC}
#SBATCH --ntasks=48
#SBATCH --time=${WALLTIME}
#SBATCH --mem=16G
#SBATCH --output=logs/tracer_${PROC}_%j.out

module load gcc
export OMP_NUM_THREADS=48
mkdir -p data
./build/tracer_${PROC} --N_traj ${N_TRAJ} --T ${T_SIM} ${EXTRA_ARGS}
EOF
    echo "Submitted tracer_${PROC}  (N_traj=${N_TRAJ}, T=${T_SIM}, walltime=${WALLTIME})"
done
