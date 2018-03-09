#! /bin/bash
# This file is is used to submit the snakemake job
# Usage: sbatch --cpus-per-task=32 --mem=64g snakemake.sh

module load snakemake samtools peakranger macs gem music PePr || exit 1

sbcmd="sbatch --cpus-per-task={threads} --mem={cluster.mem}"
sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
sbcmd+=" --out={cluster.out} {cluster.extra}"
snakemake -pr --keep-going --local-cores $SLURM_CPUS_PER_TASK \
             --jobs 10 --cluster-config cluster.json --cluster "$sbcmd" \
             --latency-wait 120 all

snakemake -pr --keep-going --local-cores $SLURM_CPUS_PER_TASK \
    --jobs 10 --cluster-config cluster.json --cluster "$sbcmd" \
    --latency-wait 120 all