#!/bin/bash
# Submit one SLURM job per swimmer process on Apocrita.
# Usage: bash scripts/submit_apocrita.sh
# Assumes binaries are already built in ./build/

mkdir -p logs

for PROC in aoup abp rtp levy1 levy2 bmshort; do
    if [ "$PROC" = "levy2" ]; then
        EXTRA_ARGS="--beta 0.5 --tau_0 1.0"
    elif [ "$PROC" = "bmshort" ]; then
        EXTRA_ARGS="--D_bm 1.0 --V_0 1.0 --sigma 1.0 --b_max 5.0"
    else
        EXTRA_ARGS=""
    fi

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=tracer_${PROC}
#SBATCH --ntasks=48
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/tracer_${PROC}_%j.out

module load gcc
export OMP_NUM_THREADS=48
mkdir -p data
./build/tracer_${PROC} --N_traj 1000000 --T 1000 ${EXTRA_ARGS}
EOF
    echo "Submitted tracer_${PROC}"
done
