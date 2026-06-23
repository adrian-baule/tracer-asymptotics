#!/bin/bash
# Submit one SLURM job per swimmer process on Apocrita.
# Usage: bash scripts/submit_apocrita.sh
# Assumes binaries are already built in ./build/

mkdir -p logs

for PROC in aoup abp rtp levy; do
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
./build/tracer_${PROC} --N_traj 1000000 --T 1000
EOF
    echo "Submitted tracer_${PROC}"
done
